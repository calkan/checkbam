#include "verifybam.h"
#include <time.h>

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
			fprintf(stdout, "verifybam server is already runnning. Initialized in client mode.\n");
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
			fprintf(stdout, "verifybam server is not present. Initialized in server mode.\n");
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
		load_bam( in_bam, params->bam_file, params->limit);

		/* Run actual verification process */
		int result = -1;
		FILE* outstream = stdout;
		if(params->output_file != NULL){
			outstream = fopen(params->output_file, "w");
		}
		result = read_alignment(in_bam, params);

		fprintf(outstream, "%d\n", result);
		return EXIT_SUCCESS;
	}

}

void init_server(parameters **params){
	bam_info* in_bam;

	unsigned int s, s2;
	struct sockaddr_un local, remote;
	int len;

	if((s = socket(AF_UNIX, SOCK_STREAM, 0)) == -1){
		perror("socket");
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
		perror("listen");
		exit(EXIT_FAILURE);
	}

	fprintf(stdout, "Started verifybam\n");

	/* Redirect STDIO to a log file*/
	switch_stdio(stdout, "out_verifybam.log");
	switch_stdio(stderr, "err_verifybam.log");

	daemon(0,1);

	char * datetime = get_datetime();
	fprintf(stdout, "STARTED verifybam daemon: %s\n", datetime);
	fprintf(stdout, "Loading reference genome\n");

	free(datetime);

	(*params)->ref_fai = fai_load((*params)->ref_genome);

	load_chrom_properties(*params);

	// Should spawn thread for each client.
	while(1){
		len = sizeof(struct sockaddr_un);
		s2 = accept(s, &remote, &len);

		int filename_len;
		recv(s2, &filename_len, sizeof(int), 0);

		char buf[filename_len];

		recv(s2, buf, filename_len, 0);
		set_str(&((*params)->bam_file), buf);
		fprintf(stdout, "Read params %s %d\n", (*params)->bam_file, filename_len);

		in_bam = ( bam_info*) malloc( sizeof( bam_info));
		in_bam->sample_name = NULL;
		load_bam( in_bam, (*params)->bam_file, (*params)->limit);
		int result = read_alignment(in_bam, (*params));

		send(s2, &result, sizeof(int), 0);

		free(in_bam);
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

	int result = -1;
	recv(s, &result, sizeof(int), 0);

	FILE* outstream = stdout;
	if(params->output_file != NULL){
		outstream = fopen(params->output_file, "w");
	}

	fprintf(outstream, "%d\n", result);

	close(s);
}

int is_daemon_running(){
	int fd = open(DAEMON_LOCK,
	    O_CREAT | //create the file if it's not present.
	    O_WRONLY,//only need write access for the internal locking semantics.
	    S_IRUSR | S_IWUSR); //permissions on the file, 600 here.

	if (fd == -1) {
	    return 1;
	}
	else if(flock(fd, LOCK_EX|LOCK_NB) == -1){
		return 1;
	}
	else{
		return 0;
	}
}

void switch_stdio(FILE * stream, const char * file_path){
	fflush(stream);
	freopen(file_path, "a+", stream);
	setlinebuf(stream);
}