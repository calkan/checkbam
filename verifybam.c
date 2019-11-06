#include "verifybam.h"
#include <time.h>

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


	/* The second argument, 1, is the number of incoming connections
	   that can be queued before you call accept(), below.
	   If there are this many connections waiting to be accepted,
	   additional clients will generate the error ECONNREFUSED. */
	if( listen(s, 1) == -1){
		perror("listen");
		exit(EXIT_FAILURE);
	}

	fprintf(stdout, "Started verifybam %s\n", get_datetime());
	fprintf(stdout, "Loading reference genome\n");

	(*params)->ref_fai = fai_load((*params)->ref_genome);

	load_chrom_properties(*params);

	fprintf(stderr, "\nWaiting for incoming tasks\n");
	signal(SIGPIPE, SIG_IGN); // Ignore pipe faults.

	// Should spawn thread for each client.
	while(1){
		len = sizeof(struct sockaddr_un);
		s2 = accept(s, &remote, &len);

		int filename_len;
		int received = recv(s2, &filename_len, sizeof(int), 0);
		if(received < 0) {
			fprintf(stderr, "Receive error\n");
			continue;
		}

		char buf[filename_len+1];

		received = recv(s2, buf, filename_len, 0); 
		if(received < 0) {
			fprintf(stderr, "Receive error\n");
			continue;
		}

		buf[filename_len] = '\0';
		set_str(&((*params)->bam_file), buf);
		fprintf(stdout, "Read params %s\n", (*params)->bam_file);

		in_bam = ( bam_info*) malloc( sizeof( bam_info));
		in_bam->sample_name = NULL;
		int load_result = load_bam( in_bam, (*params)->bam_file, (*params)->limit, (*params)->samMode);
		if(load_result < 0) {
			send(s2, &load_result, sizeof(int), 0);
		}
		else {
			verifybam_result_t* result = read_alignment(in_bam, (*params));
			destroy_bam_info(in_bam);

			send(s2, &(result->code), sizeof(int), 0);
			if(result->code >= 0) {
				len = strlen(result->hash);
				send(s2, &len, sizeof(int), 0);
				send(s2, result->hash, sizeof(char)*len, 0);
			}
			free(result->hash);
			free(result);
		}
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

	FILE* outstream = stdout;
	if(params->output_file != NULL){
		outstream = fopen(params->output_file, "w");
	}
	verifybam_result_t* result = init_verifybam_result();
	recv(s, &(result->code), sizeof(int), 0);
	if(result->code >= 0) {
		recv(s, &(len), sizeof(int), 0);
		char buf[len+1];
		recv(s, buf, len, 0); buf[len] = '\0';
		set_str(&(result->hash), buf);

		fprintf(outstream, "%d\n", result->code);
		fprintf(outstream, "%s\n", result->hash);
	}
	else {
		fprintf(outstream, "%d\n", result->code);
	}

	close(s);
}

int is_server_running(){
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
	return_value = parse_command_line( argc, argv, params);
	if( return_value == 0)
	{
		exit( EXIT_SUCCESS);
	}
	else if( return_value != 1)
	{
		exit( return_value);
	}

	if(params->mode == SERVER){
		if(is_server_running()){
			// A progress is already running.
			fprintf(stderr, "A verifybam server is already running, cannot start in server mode.\n");
			exit(EXIT_PARAM_ERROR);
		}
		else{
			// Start server mode if reference is present
			if( params->ref_genome == NULL)
			{
				fprintf( stderr, "[VERIFYBAM CMDLINE ERROR] Please enter reference genome file (FASTA) using --ref option.\n");
				exit(EXIT_PARAM_ERROR);
			}
			init_server(&params);
		}
		exit(EXIT_SUCCESS);
	}
	else if(params->mode == CLIENT) {
		if(!is_server_running()){
			// Verifybam server is not present
			fprintf(stderr, "Could not find a running verifybam server. Please start a server first\n");
			exit(EXIT_PARAM_ERROR);
		}
		else {
			if( params->bam_file == NULL)
			{
				fprintf( stderr, "[VERIFYBAM CMDLINE ERROR] Please enter input bam file using --input option.\n");
				exit(EXIT_PARAM_ERROR);
			}
			init_client(params);
		}
	}
	// Sequential mode
	else{
		// Load reference genome into memory
		params->ref_fai = fai_load(params->ref_genome);

		load_chrom_properties(params);

		in_bam = ( bam_info*) malloc( sizeof( bam_info));
		in_bam->sample_name = NULL;
		load_bam( in_bam, params->bam_file, params->limit, params->samMode);

		/* Run actual verification process */
		verifybam_result_t* result;
		FILE* outstream = stdout;
		if(params->output_file != NULL){
			outstream = fopen(params->output_file, "w");
		}
		result = read_alignment(in_bam, params);
		destroy_bam_info(in_bam);	

		fprintf(outstream, "%d\n", result->code);
		fprintf(outstream, "%s\n", result->hash);
		return EXIT_SUCCESS;
	}

}