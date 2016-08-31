#ifndef __PROCESSBAM
#define __PROCESSBAM

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <memory.h>
#include <string.h>
#include <zlib.h>
#include <pthread.h>

#include "kseq.h"
#include "sha256.h"
#include "common.h"
KSEQ_INIT(gzFile, gzread)

#define BUFFER_SIZE 10000
#define READ_SAMPLE 1
#define HASH_SAMPLE 2
#define END_SIGNAL 0

/* Sample this many fragments to calculate avg/median/std per library */
#define SAMPLEFRAG 1000000

/* Maximum sequence/quality length */
#define MAX_SEQ 1000

// Multithreading stop flag
int _stop_flag;

typedef struct _bam_info
{
	htsFile* bam_file; /* file pointer to the BAM file */
	hts_idx_t *bam_index;
	bam_hdr_t* bam_header;
	char* sample_name; /* name of the sample, parsed from SM in the BAM header */
	int num_libraries; /* number of libraries, counted from the RG tags in the BAM header */
	struct library_properties** libraries; /* each library_properties struct holds statistical/other info */
} bam_info;

typedef struct _job{
	int job_type;
	void* data;
} job_t;

typedef struct _buffer{
	job_t* stack[BUFFER_SIZE];
	size_t len;
	pthread_mutex_t mutex;
	pthread_cond_t can_produce;
  pthread_cond_t can_consume;
} buffer_t;

typedef struct _thread_data {
   int thread_id;
   buffer_t buffer;
	 BYTE hash_bam[SHA256_BLOCK_SIZE];
	 BYTE hash_fastq[SHA256_BLOCK_SIZE];
	 int aligned_read_count;
	 int hashed_read_count;
	 parameters *params;
} thread_args_t;


typedef struct _hash_buffer{
	char*	stack[BUFFER_SIZE];
	size_t len;
	int end_signal;
	pthread_mutex_t mutex;
	pthread_cond_t can_produce;
  pthread_cond_t can_consume;
} hash_buffer_t;


typedef struct _hash_thread_data{
	int thread_id;
	hash_buffer_t buffer;
	BYTE hash[SHA256_BLOCK_SIZE];
	int hashed_read_count;
	parameters *params;
} hash_thread_args_t;

void load_bam( bam_info* in_bam, char* path);
int read_alignment( bam_info* in_bam, parameters *params);
int readcmp(char* read1, char* read2);



/* BAM Utility functions */
void get_sample_name( bam_info* in_bam, char* header_text);



#endif
