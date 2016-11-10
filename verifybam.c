#include "verifybam.h"

int main( int argc, char** argv)
{
	bam_info* in_bam;
	parameters* params;
	int return_value;
	char username[MAX_SEQ];
	int i;
	int j;
	FILE *refFile;

	/* Set program parameters */
	init_params( &params);

	/* Seed random number generator */
	srand(time(NULL));

	/* Parse command line arguments */
	return_value = parse_command_line( argc, argv, params, EXE_VERIFYBAM);
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

	if(params->daemon){
		if(is_daemon_running()){
			// A progress is already running. Send request to there.
			printf("verfiybam server is already runnning. Initialized in client mode.\n");
			if( params->bam_file == NULL)
			{
				fprintf( stderr, "[VERIFYBAM CMDLINE ERROR] Please enter an input BAM file using the --input option.\n");
				exit(EXIT_PARAM_ERROR);
			}
			init_client(params);
		}
		else{
			// Verifybam is started for first time.
			// Go into daemon mode and then start working on current request.
			printf("verifybam server is not present. Initialized in server mode.\n");
			if( params->ref_genome == NULL)
			{
				fprintf( stderr, "[VERIFYBAM CMDLINE ERROR] Please enter reference genome file (FASTA) using the --ref option.\n");
				exit(EXIT_PARAM_ERROR);
			}
			init_server(&params);
		}
		exit(EXIT_SUCCESS);
	}
	else{
		params->ref_fai = fai_load(params->ref_genome);

		load_chrom_properties(params);

		in_bam = ( bam_info*) malloc( sizeof( bam_info));
		in_bam->sample_name = NULL;
		load_bam( in_bam, params->bam_file);

		/* Run actual verification process */
		return read_alignment(in_bam, params);
	}
	
}

void init_server(parameters **params){
	bam_info* in_bam;

	unsigned int s, s2;
	struct sockaddr_un local, remote;
	int len;

	if((s = socket(AF_UNIX, SOCK_STREAM, 0)) == -1){
		exit(EXIT_FAILURE);
	}

	local.sun_family = AF_UNIX;  /* local is declared before socket() ^ */
	strcpy(local.sun_path, SOCK_PATH);
	unlink(local.sun_path);
	len = strlen(local.sun_path) + sizeof(local.sun_family);
	bind(s, (struct sockaddr *)&local, len);


	/* The second argument, 5, is the number of incoming connections 
	   that can be queued before you call accept(), below. 
	   If there are this many connections waiting to be accepted, 
	   additional clients will generate the error ECONNREFUSED. */
	if( listen(s, 5) == -1){
		exit(EXIT_FAILURE);
	}

	printf("Started verifybam server\n");

	daemon(0,0);

	(*params)->ref_fai = fai_load((*params)->ref_genome);

	load_chrom_properties(*params);

	while(1){
		len = sizeof(struct sockaddr_un);
		s2 = accept(s, &remote, &len);

		int filename_len;
		recv(s2, &filename_len, sizeof(int), 0);
	
		char buf[filename_len];
		
		recv(s2, buf, filename_len, 0);
		set_str(&((*params)->bam_file), buf);
		printf("Read params %s %d\n", (*params)->bam_file, filename_len);
		
		in_bam = ( bam_info*) malloc( sizeof( bam_info));
		in_bam->sample_name = NULL;
		load_bam( in_bam, (*params)->bam_file);
		int result = read_alignment(in_bam, (*params));

		send(s2, &result, sizeof(int), 0);
	}
}

void init_client(parameters* params){
	unsigned int s;
	struct sockaddr_un local, remote;
	int len;

	if ((s = socket(AF_UNIX, SOCK_STREAM, 0)) == -1) {
        perror("socket");
        exit(1);
    }

    remote.sun_family = AF_UNIX;
    strcpy(remote.sun_path, SOCK_PATH);
    len = strlen(remote.sun_path) + sizeof(remote.sun_family);
    if (connect(s, (struct sockaddr *)&remote, len) == -1) {
        perror("connect");
        exit(1);
    }

    int bamfile_len = strlen(params->bam_file);

    if (send(s, &bamfile_len, sizeof(int), 0) == -1) {
        perror("send");
        exit(1);
    }

    if (send(s, params->bam_file, sizeof(char)*bamfile_len, 0) == -1) {
        perror("send");
        exit(1);
    }

    int result;
    recv(s, &result, sizeof(int), 0);

    FILE* outstream = stdout;
	if(params->output_file != NULL){
		 outstream = fopen(params->output_file, "w");
	}

	fprintf(outstream, "%d\n", result);

    close(s);
}