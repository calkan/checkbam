#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>

/* tardis headers */
#include "processbam.h"

void fix_n_base( char *ref, char *read){
	int i;
	int len = strlen(ref);
	for(i=0; i<len; i++){
		if(ref[i] == 'N')
			ref[i] = read[i];
	}
}

void load_bam( bam_info* in_bam, char* path)
{
	/* Variables */
	htsFile* bam_file;
	bam_hdr_t* bam_header;
	hts_idx_t* bam_index;

	fprintf( stderr, "Processing BAM file %s.\n", path);

	/* Open the BAM file for reading. htslib automatically detects the format
		of the file, so appending "b" after "r" in mode is redundant. */
	bam_file = safe_hts_open( path, "r");
	bam_index = sam_index_load( bam_file, path);

	if (bam_index == NULL){
		fprintf(stderr, "BAM index not found.\n");
		exit (EXIT_COMMON);
	}
	/* Read in BAM header information */
	bam_header = bam_hdr_read( ( bam_file->fp).bgzf);

	get_sample_name( in_bam, bam_header->text);
	in_bam->bam_file = bam_file;
	in_bam->bam_index = bam_index;
	in_bam->bam_header = bam_header;

}

void read_alignment( bam_info* in_bam, parameters *params)
{
	bam1_core_t bam_alignment_core;
	bam1_t*	bam_alignment;
	hts_idx_t* bam_index;
	int return_value;
	int i;
	int j;

	char sequence[MAX_SEQ];
	char read[MAX_SEQ];
	char qual[MAX_SEQ];
	char read2[MAX_SEQ];

	BYTE hash[MAX_SEQ] = {0x0};
	int max_readlen = 0;
	SHA256_CTX ctx;
	BYTE buf[SHA256_BLOCK_SIZE];

	int read_len;
	int ref_len;

	char next_char;
	const uint32_t *cigar;
	int n_cigar;
	htsFile *bam_file;

	bam_hdr_t* bam_header;
	char rand_loc[MAX_SEQ];
	int result;
	char md[MAX_SEQ];
	int chrom_id;
	int start; int end;
	char map_chr[MAX_SEQ];
	int map_tid;
	int map_loc;
	/* char *ref_seq; */
	char ref_seq[MAX_SEQ];
	char ref_seq2[MAX_SEQ];
	int loc_len;
	int soft_clips[2] = {0};
	int cigar_add_len;
	int clipped;
	int aligned_read_count=0;
	int reversed = 0;

	bam_file = in_bam->bam_file;
	bam_header = in_bam->bam_header;
	bam_index = in_bam->bam_index;

	bam_alignment = bam_init1();
	j=0;

	if (bam_index == NULL){
		fprintf(stderr, "BAM index not found.\n");
		exit (EXIT_COMMON);
	}

	if (bam_header == NULL){
		fprintf(stderr, "BAM header not found.\n");
		exit (EXIT_COMMON);
	}

	return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);

	while( return_value != -1 ){

		bam_alignment_core = bam_alignment->core;

		if (bam_alignment_core.flag & (BAM_FSECONDARY|BAM_FSUPPLEMENTARY|BAM_FUNMAP)) // skip secondary, supplementary and unmapped alignments
		{
			// hts_itr_destroy(iter);
			// it is ok to demux these and write out
			return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);
			continue;
		}
		reversed = bam_alignment_core.flag & BAM_FREVERSE;

		if (bam_aux_get(bam_alignment, "MD") == NULL){
				return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);
				continue;
		}

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

		int seq_len = strlen( sequence);
		for( i = 0; i < seq_len; i++){
			read[i] = base_as_char( bam_seqi( sequence, i));
			if(reversed){
				hash[(seq_len-1)-i] += complement_char(read[i]);
			}
			else{
				hash[i] += read[i];
			}
		}
		read[i] = '\0';
		read_len = i;
		if(read_len > max_readlen) max_readlen = read_len;
		strcpy(read2, read);

		clipped=0;
		cigar_add_len = 0;

		//		//fprintf(stdout, "\nCIGAR: ");
		for (i=0; i<n_cigar; i++){
			if (bam_cigar_opchr(cigar[i]) == 'H'){
				return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);

				clipped = 1;
				break;
			}
			else if (bam_cigar_opchr(cigar[i]) == 'D')
				cigar_add_len += bam_cigar_oplen(cigar[i]);
			else if (bam_cigar_opchr(cigar[i]) == 'I')
				cigar_add_len -= bam_cigar_oplen(cigar[i]);
		}

		if (clipped){
			//fprintf(stdout, "Clipped\n");
			break;
		}


		strcpy(md, bam_aux_get(bam_alignment, "MD"));

		map_tid = bam_alignment_core.tid;
		map_loc = bam_alignment_core.pos;

		//fprintf(stdout, "map_tid: %d\t map_loc: %d\n", map_tid, map_loc);

		// Get chromosome name to never use again!
		strcpy(map_chr, bam_header->target_name[map_tid]);

		//fprintf(stdout, "map_chr: %s\n", map_chr);

		strncpy(ref_seq, params->chrom_seq[map_tid]+map_loc, read_len+cigar_add_len);
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

		// Debug fix.
		// Strange condition. Ref genome has N base, MD field suggest different.
		fix_n_base(ref_seq, read);

		//		//fprintf(stdout, "\npos\n%s\n%s\n", read, ref_seq);

		if (bam_alignment_core.flag & BAM_FREVERSE) {
			printf("Reverse strand is found\n%s\n%s\n%s\n", bam_get_qname(bam_alignment), read2, ref_seq2);
		}

		if (strcmp(read, ref_seq)){
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
			fprintf(stdout, "\npos\n%s\n%s\n", read, ref_seq);
			fprintf(stdout, "\ntotal aligned reads\n%d\n", aligned_read_count);
			return;
		}
		else{
			aligned_read_count++;
			//fprintf(stdout, "\nread aligned\n%s\n%s\n", read, ref_seq);
		}

		j++;

		return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);

	}

	sha256_init(&ctx);
	sha256_update(&ctx, hash, max_readlen);
	sha256_final(&ctx, buf);

	BYTE hash2[MAX_SEQ] = {0x0};
	int max_readlen2 = 0;
	SHA256_CTX ctx2;
	BYTE buf2[SHA256_BLOCK_SIZE];
	for(i=0; i<params->num_fastq_files;i++){
		gzFile fp = gzopen(params->fastq_files[i], "r");
		kseq_t *seq = kseq_init(fp);
		int l;
		while ((l = kseq_read(seq)) >= 0) {
			for(j=0;j<strlen(seq->seq.s);j++){
				printf("%c", seq->seq.s[j]);
				hash2[j] += seq->seq.s[j];
			}
			if(max_readlen2 < strlen(seq->seq.s)) max_readlen2 = strlen(seq->seq.s);
		}
	}

	sha256_init(&ctx2);
	sha256_update(&ctx2, hash2, max_readlen2);
	sha256_final(&ctx2, buf2);

	int hash_result = 1;
	for(i=0;i<SHA256_BLOCK_SIZE;i++){
		hash_result = hash_result & (buf[i]==buf2[i]);
	}

	fprintf(stdout, "\nAll reads are matched. Total aligned read count is %d\nBamhash result is %d\n", aligned_read_count, hash_result);
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
