#include "fqhash.h"

int main( int argc, char** argv)
{
	parameters* params;
	int return_value;
	int i, j, k;

	MD5_CTX ctx;
	BYTE buf[MD5_BLOCK_SIZE];

	/* Set program parameters */
	init_params( &params);

	/* Parse command line arguments */
	return_value = parse_command_line( argc, argv, params, EXE_FQHASH);
	if( return_value == 0)
	{
		exit( EXIT_SUCCESS);
	}
	else if( return_value != 1)
	{
		exit( return_value);
	}

	if ( DEBUG)
	{
		print_params( params);
	}

	if(params->num_fastq_files <= 0){
		fprintf(stderr, "No fastq file is present to compute hash. Check arguments");
		print_params( params);
		exit( EXIT_SUCCESS);
	}

	j = 0;
	BYTE hash_fastq[MD5_BLOCK_SIZE] = {0x0};
	for(i=0; i<params->num_fastq_files;i++){
		fprintf(stdout, "Reading %s\n", params->fastq_files[i]);
		gzFile fp = gzopen(params->fastq_files[i], "r");
		kseq_t* seq = kseq_init(fp);
		int l;
		while ((l = kseq_read(seq)) >= 0) {
			int length = strlen(seq->seq.s);

			md5_init(&ctx);
			md5_update(&ctx, seq->seq.s, length);
			md5_final(&ctx, buf);

			for( k = 0; k < MD5_BLOCK_SIZE; k++){
				hash_fastq[k] += buf[k];
			}

			j++;
			if((j % 100000) == 0){
				fprintf(stdout, "Processed %d reads\n", j);
				fflush(stdout);
			}
		}
	}

	fprintf(stdout, "Hash of fastq files:\n");
	for( k = 0; k < MD5_BLOCK_SIZE; k++){
		fprintf(stdout, "%x", hash_fastq[k]);
	}
	fprintf(stdout, "\n");

	return RETURN_SUCCESS;
}
