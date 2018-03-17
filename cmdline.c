#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "verifybam.h"
#include "cmdline.h"

int parse_command_line( int argc, char** argv, parameters* params, int exe)
{
	int index,i;
	short flag;
	int o;
	static struct option long_options[] =
	{
		{"input"	, required_argument,	0, 'i'},
		{"ref"		, required_argument,	0, 'f'},
		{"fastq"	, required_argument,	0, 'q'},
		{"fq-list"	, required_argument,	0, 'l'},
		{"output"	, required_argument,	0, 'o'},
		{"limit"	, required_argument,	0, 'c'},
		{"threads"	, required_argument,	0, 't'},
		{"help"		, no_argument,			0, 'h'},
		{"server"	, no_argument,			0, 's'},
		{"daemon"	, no_argument,			0, 'd'},
		{"sam"		, no_argument,			0, 'S'},
		{"version"	, no_argument,			0, 'v'},
		{0			, 0, 					0, 0}
	};

	if( argc == 1)
	{
		print_help(exe);
		return 0;
	}

	while( ( o = getopt_long( argc, argv, "hdv:i:f:t:o:q:l:c:", long_options, &index)) != -1)
	{
		switch(o)
		{
			case 'i':
				set_str( &( params->bam_file), optarg);
			break;

			case 'f':
				set_str( &( params->ref_genome), optarg);
			break;

			case 'q':
				set_str( &( params->fastq_files[params->num_fastq_files]), optarg);
				params->num_fastq_files++;
			break;

			case 'o':
				set_str( &( params->output_file), optarg);
			break;

			case 'l':
				set_str( &( params->fastq_list), optarg);
				parse_fastq_list(params);
			break;

			case 't':
				params->threads = atoi( optarg);
			break;

			case 'h':
				print_help(exe);
				return 0;
			break;

			case 'd':
				params->daemon = 1;
			break;

			case 's':
				params->server = 1;
			break;

			case 'S':
				params->samMode = 1;
			break;

			case 'c':
				flag = 1;
				for(i=0; i<strlen(optarg); i++){
					flag &= isdigit(optarg[i]);
					if(!flag) break;
				}
				if(flag){
					params->limit = atoi(optarg);
				}
				else{
					fprintf( stderr, "\nLimit argument only accepts integers.\n");
				}
				params->limit = atoi(optarg);
			break;

			case 'v':
				if(exe == EXE_VERIFYBAM){
					fprintf( stderr, "\nVERIFYBAM: BAM validity checking tool.\n");
					fprintf( stderr, "Version %s\n\tLast update: %s, build date: %s\n\n", VERSION, UPDATE, BUILD_DATE);
					fprintf( stderr, "It is bigger on the inside!\n\n");
				}
				else if(exe == EXE_FQHASH){
					fprintf( stderr, "\nFQHASH: Multiple FASTQ file hashing tool.\n");
					fprintf( stderr, "Version %s\n\tLast update: %s, build date: %s\n\n", VERSION, UPDATE, BUILD_DATE);
					fprintf( stderr, "It is lesser on the inside!\n\n");
				}
				return 0;
			break;
		}
	}

	/* TODO: check parameter validity */

	if(exe == EXE_VERIFYBAM && !(params->server)){
		if( params->bam_file == NULL)
		{
			fprintf( stderr, "[VERIFYBAM CMDLINE ERROR] Please enter an input BAM file using the --input option.\n");
			return EXIT_PARAM_ERROR;
		}
		/* check if --ref	 is invoked */
		if( params->ref_genome == NULL)
		{
			fprintf( stderr, "[VERIFYBAM CMDLINE ERROR] Please enter reference genome file (FASTA) using the --ref option.\n");
			return EXIT_PARAM_ERROR;
		}
	}
	else if(exe == EXE_VERIFYBAM){
		if(!(params->limit > 0 && params->limit <= 100)){
			fprintf( stderr, "Limit must be in (0,100] range.\n");
			return EXIT_PARAM_ERROR;
		}
	}

	if(exe == EXE_FQHASH){
		if( params->num_fastq_files == 0)
		{
			fprintf( stderr, "[FQHASH CMDLINE ERROR] Please enter at least one fastq file (FASTQ) using the --fastq option.\n");
			return EXIT_PARAM_ERROR;
		}
	}

	/* check if threads>0 */
	if( params->threads <= 0)
	{
		fprintf( stderr, "[CMDLINE WARNING] Invalid number of threads was entered (%d). Resetted to 1.\n", params->threads);
		params->threads = 1;
	}

	return 1;

}

void parse_fastq_list( parameters* params){
	char* file_address = malloc(sizeof(char)*200);
	size_t len;
	ssize_t read;
	FILE *f = fopen(params->fastq_list, "r");
	if(f == NULL){
		// One of fastq files could not be read. Report the incident.
	}
	while((read=getline(&file_address, &len, f)) != -1){
		// Remove trailing new line
		file_address[strlen(file_address)-1]='\0';
		set_str( &( params->fastq_files[params->num_fastq_files]), file_address);
		params->num_fastq_files++;
	}
	fclose(f);

	if(file_address)
		free(file_address);
}

void print_help( int exe)
{
	if(exe == EXE_VERIFYBAM){
		fprintf( stdout, "\nVERIFYBAM: BAM validity checking tool.\n");
		fprintf( stdout, "Version %s\n\tLast update: %s, build date: %s\n\n", VERSION, UPDATE, BUILD_DATE);
		fprintf( stdout, "\t--input  [BAM file]  : Input file in sorted and indexed BAM format.\n");
		fprintf( stdout, "\t--ref    [reference] : Reference genome in FASTA format.\n");
		fprintf( stdout, "\t--output [File]      : Write result code and BamHash to an output file.\n");
		fprintf( stdout, "\t--threads            : Number of threads to run.\n");
		fprintf( stdout, "\t--daemon             : Start or interact with daemon mode.\n");
		fprintf( stdout, "\t--version            : Print version and exit.\n");
		fprintf( stdout, "\t--help               : Print this help screen and exit.\n\n");
	}
	else if(exe == EXE_FQHASH){
		fprintf( stdout, "\nFQHASH: Multiple FASTQ file hashing tool.\n");
		fprintf( stdout, "Version %s\n\tLast update: %s, build date: %s\n\n", VERSION, UPDATE, BUILD_DATE);
		fprintf( stdout, "\t--fastq [Fastq file] : A fastq file.\n");
		fprintf( stdout, "\t--fq-list [file]  : A file that contains addresses of multiple fastq files.\n");
		fprintf( stdout, "\t--output [file]      : Output file to write final hash.\n");
		fprintf( stdout, "\t--version            : Print version and exit.\n");
		fprintf( stdout, "\t--help               : Print this help screen and exit.\n\n");
	}
}
