#ifndef __COMMON
#define __COMMON

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <errno.h>
#include <openssl/sha.h>
#include "sha2.h"

#define BYTE unsigned char

/* Exit Codes */
#define EXIT_SUCCESS 0
#define EXIT_COMMON 1
#define EXIT_MAXBAMS 2
#define EXIT_PARAM_ERROR 3
#define EXIT_EXTERNAL_PROG_ERROR 4
#define EXIT_FILE_OPEN_ERROR 5
#define EXIT_READGROUP 6

/* Return Codes */
#define RETURN_SUCCESS 0
#define RETURN_ERROR 1

#define MAX_BAMS 256
#define MAX_FASTQS 1000

// Track memory usage
extern long long memUsage;

typedef enum { MALE, FEMALE}gender;
typedef enum { SERVER = 0, CLIENT = 1, SEQUENTIAL = 2 }runMode;

typedef struct _params
{
	char* ref_genome; /* path to reference genome - fasta */
	char* bam_file; /* path to input bam file */
	int threads; /* number of threads to use for parallel mrFAST, and maybe future parallelization of TARDIS */
	int num_chrom; /* number of chromosomes */
	int* chrom_lengths; /* lengths of the chromosomes */
	char** chrom_names; /* names of the chromosomes */
	char **chrom_seq; /* chromosomes */
	faidx_t* ref_fai;
	runMode mode; /* Running mode */
	short limit; /* Bottom limit for aligned read count in percentage. */
	short samMode; /* Input is in SAM format */
    short hashing_enabled; /* Whether hash calculation is enabled */

	char* output_file;
} parameters;

typedef struct operation_s {
	unsigned short int start;
	unsigned short int length;
} operation;

/* Parameter related TARDIS functions */
void init_params( parameters**);
void print_params( parameters*);
void load_chrom_properties(parameters*);

/* FILE opening and error printing functions. For opening regular and BAM/SAM
 files safely */
void print_error( char*);
FILE* safe_fopen( char* path, char* mode);
gzFile safe_fopen_gz( char* path, char* mode);
htsFile* safe_hts_open( char* path, char* mode);

/* General BAM processing functions */
int is_concordant( bam1_core_t bam_alignment_core, int min, int max);
char base_as_char( int base_as_int);
int char_as_base( char base);
char complement_char( char base);
void qual_to_ascii( char* qual);

/* String functions */
void set_str( char **target, char *source); /* Even safer than strncpy */
void reverse_string( char* str);

/* Misc. Utility */
int compare_size_int( const void* p, const void* q);
void print_quote( void);

// Memory allocation/tracking functions
void* getMem( size_t size);
void freeMem( void* ptr, size_t size);
double getMemUsage();

void del_char(char *ref, int start, int len);
void ins_char(char *ref, char *read, int start_origin, int start_dest, int len);
//void applymd(char *ref, char *md);

void init_sha256_block(BYTE**);
void sha256_hash(char*, BYTE**);

void apply_cigar_md(char *ref, char *read, char *md, int n_cigar, const uint32_t *cigar);
char* get_datetime();

#endif
