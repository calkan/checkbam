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


	/* The second argument, 5, is the number of incoming connections
	   that can be queued before you call accept(), below.
	   If there are this many connections waiting to be accepted,
	   additional clients will generate the error ECONNREFUSED. */
	if( listen(s, 5) == -1){
		perror("listen");
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "Started verifybam %s\n", get_datetime());
	fprintf(stderr, "Loading reference genome\n");

	(*params)->ref_fai = fai_load((*params)->ref_genome);

	load_chrom_properties(*params);

	fprintf(stderr, "\nWaiting for incoming tasks\n");
	signal(SIGPIPE, SIG_IGN); // Ignore pipe faults.

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
		set_str(&((*params)->job_dir), buf);
		fprintf(stderr, "Job directory: %s\n", (*params)->job_dir);

		in_bam = ( bam_info*) malloc( sizeof( bam_info));
		in_bam->sample_name = NULL;

		char bam_file_path[1000];
		sprintf(bam_file_path, "%s/upload/output.sorted.bam", (*params)->job_dir);
		int load_result = load_bam( in_bam, bam_file_path, (*params)->limit, 0);
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

	int job_dir_len = strlen(params->job_dir);

	if (send(s, &job_dir_len, sizeof(int), 0) == -1) {
		perror("send");
		exit(1);
	}

	if (send(s, params->job_dir, sizeof(char)*job_dir_len, 0) == -1) {
		perror("send");
		exit(1);
	}

	verifybam_result_t* result = init_verifybam_result();
	recv(s, &(result->code), sizeof(int), 0);
	if(result->code >= 0) {
		recv(s, &(len), sizeof(int), 0);
		char buf[len+1];
		recv(s, buf, len, 0); buf[len] = '\0';
		set_str(&(result->hash), buf);

		fprintf(stdout, "%d\n", result->code);
		fprintf(stdout, "%s\n", result->hash);
	}
	else {
		fprintf(stdout, "%d\n", result->code);
		fprintf(stdout, "INVALID_HASH\n");
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

void switch_stdio(FILE * stream, const char * file_path){
	fflush(stream);
	freopen(file_path, "a+", stream);
	setlinebuf(stream);
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
	return_value = parse_command_line( argc, argv, params, EXE_VERIFYBAM);
	if( return_value == 0)
	{
		exit( EXIT_SUCCESS);
	}
	else if( return_value != 1)
	{
		exit( return_value);
	}

	if(params->server){
		if(is_server_running()){
			// A progress is already running. Send request to there.
			fprintf(stderr, "verifybam server is already runnning. Initialized in client mode.\n");
			exit(EXIT_PARAM_ERROR);
		}
		else{
			// Verifybam is started in server mode.
			fprintf(stderr, "verifybam server is not present. Initialized in server mode.\n");
			init_server(&params);
		}
		exit(EXIT_SUCCESS);
	}
	else{
		if(!is_server_running()){
			// Server is not present
			fprintf(stderr, "verifybam server is not runnning. Client mode cannot function.\n");
			exit(EXIT_PARAM_ERROR);
		}
		init_client(params);
	}

}