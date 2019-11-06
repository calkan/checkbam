#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "verifybam.h"
#include "cmdline.h"

int parse_command_line( int argc, char** argv, parameters* params)
{
	int index,i;
	short flag;
	int o;
	static struct option long_options[] =
	{
		{"input"	, required_argument,	0, 'i'},
		{"ref"		, required_argument,	0, 'f'},
		{"output"	, required_argument,	0, 'o'},
		{"limit"	, required_argument,	0, 'c'},
		{"threads"	, required_argument,	0, 't'},
		{"mode"		, required_argument,	0, 'm'},
		{"help"		, no_argument,			0, 'h'},
        {"hash"     , no_argument,          0, 'H'},
		{"daemon"	, no_argument,			0, 'd'},
		{"sam"		, no_argument,			0, 'S'},
		{"version"	, no_argument,			0, 'v'},
		{0			, 0, 					0, 0}
	};

	if( argc == 1)
	{
		print_help();
		return 0;
	}

	while( ( o = getopt_long( argc, argv, "hdv:i:f:t:m:o:q:l:c:", long_options, &index)) != -1)
	{
		switch(o)
		{
			case 'i':
				set_str( &( params->bam_file), optarg);
			break;

			case 'f':
				set_str( &( params->ref_genome), optarg);
			break;

			case 'o':
				set_str( &( params->output_file), optarg);
			break;

			case 't':
				params->threads = atoi( optarg);
			break;

			case 'h':
				print_help();
				return 0;
			break;

            case 'H':
                params->hashing_enabled = 1;
            break;

			case 'm':
				if(strcmp(optarg, "server") == 0) {
					params->mode = SERVER;
				} else if(strcmp(optarg, "client") == 0) {
					params->mode = CLIENT;
				}
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
				fprintf( stderr, "\nVERIFYBAM: BAM validity checking tool.\n");
				fprintf( stderr, "Version %s\n\tLast update: %s, build date: %s\n\n", VERSION, UPDATE, BUILD_DATE);
				fprintf( stderr, "It is bigger on the inside!\n\n");
				return 0;
			break;
		}
	}

	/* TODO: check parameter validity */

	if(params->mode == SEQUENTIAL || params->mode == SERVER){
		if( params->ref_genome == NULL)
		{
			fprintf( stderr, "[VERIFYBAM CMDLINE ERROR] Please enter reference genome file (FASTA) using the --ref option.\n");
			return EXIT_PARAM_ERROR;
		}
	}
	if(params->mode == SEQUENTIAL || params->mode == CLIENT) {
		if( params->bam_file == NULL)
		{
			fprintf( stderr, "[VERIFYBAM CMDLINE ERROR] Please enter an input BAM file using the --input option.\n");
			return EXIT_PARAM_ERROR;
		}
	}
	if(!(params->limit > 0 && params->limit <= 100)){
		fprintf( stderr, "Limit must be in (0,100] range.\n");
		return EXIT_PARAM_ERROR;
	}
	if( params->threads <= 0)
	{
		fprintf( stderr, "[CMDLINE WARNING] Invalid number of threads was entered (%d). Resetted to 1.\n", params->threads);
		params->threads = 1;
	}

	return 1;

}

void print_help()
{
	fprintf( stdout, "\nVERIFYBAM: BAM validity checking tool.\n");
	fprintf( stdout, "Version %s\n\tLast update: %s, build date: %s\n\n", VERSION, UPDATE, BUILD_DATE);
	fprintf( stdout, "\t--input  [BAM file]  : Input file in sorted and indexed BAM format.\n");
	fprintf( stdout, "\t--ref    [reference] : Reference genome in FASTA format.\n");
	fprintf( stdout, "\t--mode               : Running mode. Default is sequential\n");
	fprintf( stdout, "\t\tserver       : Start verifybam and wait for tasks to process. Reference is required\n");
	fprintf( stdout, "\t\tclient       : Send a task to running verifybam server. Bam file is required\n");
	fprintf( stdout, "\t\tsequential   : Run verifybam once. Reference and Bam file is required\n");
	fprintf( stdout, "\t--threads            : Number of threads to run while processing bam file.\n");
	fprintf( stdout, "\t--client             : Run in client mode. Only bam file is required\n");
	fprintf( stdout, "\t--version            : Print version and exit.\n");
	fprintf( stdout, "\t--help               : Print this help screen and exit.\n\n");
}
