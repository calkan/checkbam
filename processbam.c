#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>

#include "processbam.h"

void init_queue(Queue* queue){
	if(queue == NULL){
		queue = malloc(sizeof(Queue));
	}
	queue->head = NULL;
	queue->tail = NULL;
	queue->len = 0;
}

void push(Queue* queue, job_t* job){
	Node* node = malloc(sizeof(Node));
	node->job = job;
	node->next = NULL;

	if(queue->tail == NULL){
		queue->tail = node;
		queue->head = node;
	}
	else{
		queue->tail->next = node;
		queue->tail = node;
	}
	queue->len += 1;
}

void pop(Queue* queue, job_t** job){
	// Queue is empty
	if(queue->head == NULL){
		(*job) = NULL;
		return;
	}
	// Queue has only one element
	else if(queue->head->next == NULL){
		queue->tail = NULL;
	}
	(*job) = queue->head->job;
	Node *temp = queue->head;
	queue->head = queue->head->next;
	queue->len--;

	if(temp!=NULL) {
		free(temp);
	}
}

void destroy_job(job_t* job) {
	if(job->job_type == READ_SAMPLE && job->data != NULL) {
		bam_destroy1((bam1_t*)(job->data));
	}
}

void destroy_thread_args(thread_args_t* args) {
	if(args != NULL) {
		if(args->hash_bam != NULL) {
			free(args->hash_bam);
		}
		free(args);
	}
}

void destroy_bam_info(bam_info* in_bam) {
	if(in_bam != NULL) {
		if(in_bam->bam_header != NULL) {
			bam_hdr_destroy(in_bam->bam_header);
		}
		if(in_bam->bam_index != NULL) {
			hts_idx_destroy(in_bam->bam_index);
		}
		if(in_bam->bam_file != NULL) {
			hts_close(in_bam->bam_file);
		}
		free(in_bam);	
	}
}

int load_bam( bam_info* in_bam, char* path, int limit, int samMode)
{
	/* Variables */
	htsFile* bam_file;
	bam_hdr_t* bam_header;
	hts_idx_t* bam_index;

	fprintf( stderr, "Processing BAM file %s\n", path);

	/* Open the BAM file for reading. htslib automatically detects the format
		of the file, so appending "b" after "r" in mode is redundant. */
	bam_file = safe_hts_open( path, "r");
	if(!bam_file) {
		return -1;
	}
	bam_index = sam_index_load( bam_file, path);

	/* Read in BAM header information */
	if(samMode) {
		bam_header = sam_hdr_read( bam_file);
	}
	else{
		bam_header = bam_hdr_read( ( bam_file->fp).bgzf);
	}

	if(!bam_header) {
		return -1;
	}

	if (bam_index == NULL){
		fprintf(stderr, "BAM index not found. Limit option will not work!\nSkipping stats.\n");
		in_bam->bam_index = bam_index;
	}
	else {
		int i;
		float t_mapped=0, t_unmapped=0;
		for(i=0;i<bam_header->n_targets;i++){
			uint64_t _m, _u;
			hts_idx_get_stat(bam_index, i, &_m, &_u);
			t_mapped += _m;
			t_unmapped += _u;
		}
		
		if((t_mapped/(t_mapped+t_unmapped)) < ((float)limit/100)){
			fprintf(stderr, "Mapped read percentage is under the limit: %0.f\t%0.f\t%.2f\n", t_mapped, t_unmapped, (t_mapped/(t_mapped+t_unmapped))*100);
			return -1;
		}
		in_bam->bam_index = bam_index;
	}

	get_sample_name( in_bam, bam_header->text);
	in_bam->bam_file = bam_file;
	in_bam->bam_header = bam_header;
	return 1;
}

void *read_thread(void *_args)
{
	int i,j,k;
	bam1_t* bam_alignment;
	bam1_core_t bam_alignment_core;

	char sequence[MAX_SEQ];
	char read[MAX_SEQ];
	char qual[MAX_SEQ];
	char read2[MAX_SEQ];
	BYTE hash_input_read[MAX_SEQ];

	int read_len;
	int ref_len;

	char next_char;
	const uint32_t *cigar;
	int n_cigar;

	char rand_loc[MAX_SEQ];
	int result;
	char md[MAX_SEQ];
	int chrom_id;
	int start; int end;
	char map_chr[MAX_SEQ];
	int map_tid;
	int map_loc;
	char ref_seq[MAX_SEQ];
	char ref_seq2[MAX_SEQ];
	int loc_len;
	int soft_clips[2] = {0};
	int cigar_add_len;
	int clipped;
	int reversed = 0;
	int read_counter = 0;

	char* hash_read;

	job_t* job;

	thread_args_t* args = (thread_args_t*) _args;
	while(_stop_flag==0){
		pthread_mutex_lock(&(args->buffer.mutex));

		if(args->buffer.queue.len == 0){
			pthread_cond_wait(&(args->buffer.can_consume), &(args->buffer.mutex));
		}

		pop( &(args->buffer.queue), &job);

		pthread_cond_signal(&(args->buffer.can_produce));
		pthread_mutex_unlock(&(args->buffer.mutex));

		if(job->job_type == END_SIGNAL){
			destroy_job(job);
			break;
		}
		else if(job->job_type == READ_SAMPLE){
			bam_alignment = (bam1_t*) job->data;
			if(bam_alignment!=NULL && _stop_flag==0){

				bam_alignment_core = bam_alignment->core;

				if (bam_alignment_core.flag & (BAM_FSECONDARY|BAM_FSUPPLEMENTARY)) // skip secondary, supplementary reads before hashing
				{
					args->aligned_read_count++;
					continue;
				}

				reversed = bam_alignment_core.flag & BAM_FREVERSE;

				// Pull out the cigar field from the alignment
				cigar = bam_get_cigar(bam_alignment);
				// Number of cigar operations
				n_cigar = bam_alignment_core.n_cigar;

				// Copy the whole sequence into char array
				strncpy( sequence, bam_get_seq( bam_alignment), bam_alignment_core.l_qseq);
				sequence[bam_alignment_core.l_qseq] = '\0';
				// Copy the quality string
				strncpy( qual, bam_get_qual( bam_alignment), bam_alignment_core.l_qseq);
				qual[bam_alignment_core.l_qseq] = '\0';
				qual_to_ascii(qual);

				int seq_len = strlen(sequence);
				for( i = 0; i < seq_len; i++){
					read[i] = base_as_char(bam_seqi(sequence, i));
					if(reversed){
						hash_input_read[(seq_len-1)-i] = complement_char(read[i]);
					}
					else{
						hash_input_read[i] = read[i];
					}
				}
				read[i] = '\0';
				read_len = i;

				sha256_hash(hash_input_read, &(args->hash_bam));

				strcpy(read2, read);

				if(bam_alignment_core.flag & BAM_FUNMAP){
					continue;
				}
				// Need to look after hashing
				if (bam_aux_get(bam_alignment, "MD") == NULL){
					continue;
				}

				clipped=0;
				cigar_add_len = 0;

				//		//fprintf(stdout, "\nCIGAR: ");
				for (i=0; i<n_cigar; i++){
					if (bam_cigar_opchr(cigar[i]) == 'H'){
						clipped = 1;
						break;
					}
					else if (bam_cigar_opchr(cigar[i]) == 'D')
						cigar_add_len += bam_cigar_oplen(cigar[i]);
					else if (bam_cigar_opchr(cigar[i]) == 'I')
						cigar_add_len -= bam_cigar_oplen(cigar[i]);
				}

				if (clipped){
					fprintf(stdout, "Hard Clipped\n");
					_stop_flag = args->thread_id;
					continue;
				}

				strcpy(md, bam_aux_get(bam_alignment, "MD"));

				map_tid = bam_alignment_core.tid;
				map_loc = bam_alignment_core.pos;

				//fprintf(stdout, "map_tid: %d\t map_loc: %d\n", map_tid, map_loc);

				// Get chromosome name to never use again!
				//strcpy(map_chr, bam_header->target_name[map_tid]);

				//fprintf(stdout, "map_chr: %s\n", map_chr);

				strncpy(ref_seq, args->params->chrom_seq[map_tid]+map_loc, read_len+cigar_add_len);
				ref_len = strlen(ref_seq);

				// Sometimes ref chromosome is not long enough to cover the match. This incident should be reported.
				// For now, fill the rest with N bases.
				if(ref_len < (read_len+cigar_add_len)){
					for(i=ref_len; i<read_len+cigar_add_len; i++){
						ref_seq[i] = 'N';
					}
				}
				ref_seq[read_len+cigar_add_len] = '\0';
				strcpy(ref_seq2, ref_seq);

				apply_cigar_md(ref_seq, read, md+1, n_cigar, cigar);

				if (readcmp(read, ref_seq)){
					fprintf(stdout, "%s\n", bam_get_qname(bam_alignment));
					fprintf(stdout, "%s\n", read);
					fprintf(stdout, "%s\n",	qual);

					fprintf(stdout, "n_cigar: %d\n", n_cigar);
					for (i=0; i<n_cigar; i++){
						fprintf(stdout, "%d\t%c\t%d\t", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]), bam_cigar_type(cigar[i]));
						fprintf(stdout, "%d%c\n", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]));
					}
					fprintf(stdout, "MD: %s\n", md);

					fprintf(stdout, "\npre\n%s\n%s\n", read2, ref_seq2);
					fprintf(stdout, "\npos\n%lu %s\n%lu %s\n", strlen(read), read, strlen(ref_seq), ref_seq);
					fprintf(stdout, "\ntotal aligned reads\n%d\n", args->aligned_read_count);
					_stop_flag = args->thread_id;
					continue;
				}
				else{
					args->aligned_read_count++;
				}

				/*if((args->aligned_read_count % 500000) == 0){
					fprintf(stdout, "Thread %d: number of processed reads is %d\nCurrent buffer size is %d\n", args->thread_id, args->aligned_read_count, args->buffer.queue.len);
					fflush(stdout);
				}*/
			}
			destroy_job(job);
		}
	}
	pthread_exit(NULL);
}

void send_job(job_t* job, thread_args_t* args){
	pthread_mutex_lock(&(args->buffer.mutex));
	if(args->buffer.queue.len == BUFFER_SIZE) {
		pthread_cond_wait(&(args->buffer.can_produce), &(args->buffer.mutex));
	}

	push( &(args->buffer.queue), job);

	pthread_cond_signal(&(args->buffer.can_consume));
	pthread_mutex_unlock(&(args->buffer.mutex));
}

verifybam_result_t* read_alignment( bam_info* in_bam, parameters *params)
{
	hts_idx_t* bam_index;
	bam_hdr_t* bam_header;
	htsFile *bam_file;
	int return_value;
	int i,j,k,t;
	int aligned_read_count = 0;
	pthread_t threads[params->threads];
	thread_args_t* args[params->threads];
	verifybam_result_t* result = init_verifybam_result();

	_stop_flag = 0;

	bam_file = in_bam->bam_file;
	bam_header = in_bam->bam_header;
	bam_index = in_bam->bam_index;

	bam1_t*	bam_alignment;

	if (bam_index == NULL){
		fprintf(stderr, "BAM index not found.\n");
//		result->code = EXIT_COMMON;
//		return(result);
	}

	if (bam_header == NULL){
		fprintf(stderr, "BAM header not found.\n");
		result->code = EXIT_COMMON;
		return(result);
	}

	for(t=0; t<params->threads; t++){
		args[t] = (thread_args_t*) malloc(sizeof(thread_args_t));
		args[t]->thread_id = t;
		args[t]->aligned_read_count = 0;
		args[t]->params = params;
		init_sha256_block( &(args[t]->hash_bam));
		init_queue( &(args[t]->buffer.queue));
		pthread_mutex_init(&(args[t]->buffer.mutex), NULL);
		pthread_cond_init(&(args[t]->buffer.can_produce), NULL);
		pthread_cond_init(&(args[t]->buffer.can_consume), NULL);

		int rc = pthread_create(&threads[t], NULL, read_thread, (void *)args[t]);
		if (rc){
			fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", rc);
			result->code = -1;
			return(result);
		}
	}

	struct timespec start, finish;
	struct timespec io_start, io_finish;
	double elapsed;

	clock_gettime(CLOCK_MONOTONIC, &start);
	j=0;
	bam_alignment = bam_init1();
	if(params->samMode) {
		return_value = sam_read1( bam_file, bam_header, bam_alignment);
	}
	else{
		return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);
	}

	while( return_value != -1 && _stop_flag==0 ){

		int index = j%(params->threads);
		job_t* new_job = malloc(sizeof(job_t));
		new_job->job_type = READ_SAMPLE;
		new_job->data = bam_alignment;

		send_job(new_job, args[index]);

		bam_alignment = bam_init1();
		if(params->samMode) {
			return_value = sam_read1( bam_file, bam_header, bam_alignment);
		}
		else{
			return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);
		}

		j++;

		if((j % 500000) == 0){
			fprintf(stdout, "Total IO READ COUNT %d\r", j);
			fflush(stdout);
		}
	}

	BYTE *hash_bam;
	init_sha256_block(&hash_bam);
	void* status;

	for(t=0; t<params->threads; t++){
		job_t* new_job = malloc(sizeof(job_t));
		new_job->job_type = END_SIGNAL;
		send_job(new_job, args[t]);

		pthread_join(threads[t], &status);

		aligned_read_count += args[t]->aligned_read_count;
		for( k = 0; k < SHA256_DIGEST_LENGTH; k++){
			hash_bam[k] += args[t]->hash_bam[k];
		}
		destroy_thread_args(args[t]);
	}

	clock_gettime(CLOCK_MONOTONIC, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

	fprintf(stdout, "\nAll reads are matched.\nTotal aligned read count is %d <> %d\n", aligned_read_count, j);
	fprintf(stdout, "It took %f seconds to finish\n", elapsed);

	for( k = 0; k < SHA256_DIGEST_LENGTH; k++){
		sprintf(result->hash + strlen(result->hash), "%02x", hash_bam[k]);
	}

	//free(hash_bam);
	result->code = EXIT_SUCCESS;
	return result;
}


void get_sample_name( bam_info* in_bam, char* header_text)
{
	/* Delimit the BAM header text with tabs and newlines */

	char *tmp_header = NULL;
	set_str( &( tmp_header), header_text);
	char* p = strtok( tmp_header, "\t\n");
	char sample_name_buffer[1024];

	while( p != NULL)
	{
		/* If the current token has "SM" as the first two characters,
			we have found our Sample Name */
		if( p[0] == 'S' && p[1] == 'M')
		{
			/* Get the Sample Name */
			strncpy( sample_name_buffer, p + 3, strlen( p) - 3);

			/* Add the NULL terminator */
			sample_name_buffer[strlen( p) - 3] = '\0';

			/* Exit loop */
			break;
		}
		p = strtok( NULL, "\t\n");
	}

	set_str( &( in_bam->sample_name), sample_name_buffer);
	free( tmp_header);
}

int readcmp(char* read1, char* read2){
	int len1 = strlen(read1);
	int len2 = strlen(read2);

	if(len1!=len2){
		return 1;
	}

	int i;
	for(i=0; i<len1; i++){
		if(char_as_base(read1[i]) & char_as_base(read2[i]) <= 0){
			fprintf(stdout, "Not equal: %c %c\n", read1[i], read2[i]);
			return 1;
		}
	}
	return 0;
}

verifybam_result_t* init_verifybam_result() {
	verifybam_result_t* result = malloc(sizeof(verifybam_result_t));
	result->code = -1;
	result->hash = malloc(sizeof(char)*200);
	result->hash[0] = '\0';
	return result;
}