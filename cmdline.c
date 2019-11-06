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
		{"ref"		, required_argument,	0, 'f'},
		{"job-dir"	, required_argument,	0, 'j'},
		{"limit"	, required_argument,	0, 'c'},
		{"threads"	, required_argument,	0, 't'},
		{"help"		, no_argument,			0, 'h'},
		{"server"	, no_argument,			0, 's'},
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
			case 'j':
				set_str( &( params->job_dir), optarg);
			break;

			case 'f':
				set_str( &( params->ref_genome), optarg);
			break;

			case 't':
				params->threads = atoi( optarg);
			break;

			case 'h':
				print_help(exe);
				return 0;
			break;

			case 's':
				params->server = 1;
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

	if(!(params->server) && params->job_dir == NULL){ // Client mode and job directory is not given
		fprintf( stderr, "[VERIFYBAM CMDLINE ERROR] Client mode must define a job directory to work on.\n");
		return EXIT_PARAM_ERROR;
	}
	else if(params->server && params->ref_genome == NULL){ // Server mode and reference is not given
		fprintf( stderr, "[VERIFYBAM CMDLINE ERROR] Server mode must define a reference file to work on.\n");
		return EXIT_PARAM_ERROR;
	}

	/* check if threads>0 */
	if( params->threads <= 0)
	{
		fprintf( stderr, "[CMDLINE WARNING] Invalid number of threads was entered (%d). Resetted to 1.\n", params->threads);
		params->threads = 1;
	}

	return 1;
}

void print_help( int exe)
{
	fprintf( stdout, "\nVERIFYBAM: BAM validity checking tool.\n");
	fprintf( stdout, "Version %s\n\tLast update: %s, build date: %s\n\n", VERSION, UPDATE, BUILD_DATE);
	fprintf( stdout, "\t--job-dir  [path]     : Job directory .\n");
	fprintf( stdout, "\t--ref     [reference] : Reference genome in FASTA format.\n");
	fprintf( stdout, "\t--server              : Server or client mode.\n");
	fprintf( stdout, "\t--threads             : Number of threads to run.\n");
	fprintf( stdout, "\t--version             : Print version and exit.\n");
	fprintf( stdout, "\t--help                : Print this help screen and exit.\n\n");
}
