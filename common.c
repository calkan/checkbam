#include <stdarg.h>
#include <ctype.h>

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "common.h"

// Track memory usage
long long memUsage = 0;

void init_params( parameters** params)
{
	int i;
	/* initialize parameters */
	*params = ( parameters*) malloc( sizeof( parameters));
	( *params)->ref_genome = NULL;
	( *params)->bam_file = NULL;
	( *params)->num_fastq_files = 0;
	( *params)->threads = 1;
	( *params)->daemon = 0;
	( *params)->server = 0;
	( *params)->samMode = 0;
	( *params)->limit = 90;
}

void load_chrom_properties(parameters* params)
{
	FILE* fai_file;
	int ln_count=0,i,c;
	char filename[255],first_arg[255],sec_arg[255],cache_filename[255];
	int return_value;
	int loc_len;

	sprintf( filename, "%s.fai", (params)->ref_genome);
	fai_file = safe_fopen( filename, "r");

	//count the number of chromosomes by counting the non-empty lines in the .fai file
	do{
		c=fgetc(fai_file);
		if(c=='\n') ln_count++;

	}while (c!=EOF);
	params->num_chrom=ln_count;

	/* Reset the file pointer to the start of the file */
	rewind(fai_file);

	//We need the names and lengths of each chromosome
	params->chrom_lengths = ( int*) malloc( params->num_chrom * sizeof( int));
	params->chrom_names = ( char**) malloc( params->num_chrom * sizeof( char*));
	params->chrom_seq = ( char**) malloc( params->num_chrom * sizeof( char*));
	for( i = 0; i < params->num_chrom; i++)
	{
		return_value = fscanf( fai_file, "%[^\t]\t%[^\t]%*[^\n]",first_arg,sec_arg);
		params->chrom_names[i]=NULL;
		if (first_arg[0] == '\n' || first_arg[0] == '\r')
			set_str( &(params->chrom_names[i]), first_arg+1);
		else
			set_str( &(params->chrom_names[i]), first_arg);
		params->chrom_lengths[i]=atoi(sec_arg);
		fprintf(stderr, "\rLoading chromosome  %s ",params->chrom_names[i]);
		params->chrom_seq[i] = faidx_fetch_seq( params->ref_fai, params->chrom_names[i], 0, params->chrom_lengths[i]-1, &loc_len);
	}

	fclose(fai_file);
}

void print_params( parameters* params)
{
	int i;
	fprintf(stdout, "BAM input: %s\n", params->bam_file);
	fprintf(stdout, "ref_genome: %s\n", params->ref_genome);
}

void print_error( char* msg)
{
	/* print error message */
	fprintf( stderr, "\n%s\n", msg);
	fprintf( stderr, "Invoke parameter -h for help.\n");
	//exit( EXIT_COMMON);
}


FILE* safe_fopen( char* path, char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
	FILE* file;
	char err[500];

	file = fopen( path, mode);
	if( !file)
	{
		sprintf( err, "[INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);

	}
	return file;
}

gzFile safe_fopen_gz( char* path, char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
        gzFile file;
	char err[500];

	file = gzopen( path, mode);
	if( !file)
	{
		sprintf( err, "[INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);
	}
	return file;
}

htsFile* safe_hts_open( char* path, char* mode)
{
	htsFile* bam_file;
	char err[500];

	bam_file = hts_open( path, mode);
	if( !bam_file)
	{
		sprintf( err, "[INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);
	}

	return bam_file;
}

int is_concordant( bam1_core_t bam_alignment_core, int min, int max)
{
	int flag = bam_alignment_core.flag;

	if( ( flag & BAM_FPAIRED) == 0)
	{
		/* Read is single-end. Skip this by calling it concordant */
		return 1;
	}

	if( ( flag & BAM_FPROPER_PAIR) == 0)
	{
		/* Not proper pair */
		return 0;
	}

	if( ( flag & BAM_FUNMAP) != 0)  // c.a.
	{
		/* Read unmapped; Orphan or OEA */
		return 0;
	}

	if( ( flag & BAM_FMUNMAP) != 0) // c.a.
	{
		/* Mate unmapped; Orphan or OEA */
		return 0;
	}

	if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) != 0)
	{
		/* -- orientation = inversion */
		return 0;
	}

	if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) == 0)
	{
		/* ++ orientation = inversion */
		return 0;
	}

	if( bam_alignment_core.tid != bam_alignment_core.mtid)
	{
		/* On different chromosomes */
		return 0;
	}

	if( bam_alignment_core.pos <= bam_alignment_core.mpos) // c.a.
	{
		/* Read is placed BEFORE its mate */
		if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) == 0)
		{
			/* -+ orientation = tandem duplication */
			return 0;
		}
	}
	else
	{
		/* Read is placed AFTER its mate */
		if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) != 0)
		{
			/* +- orientation = tandem duplication */
			return 0;
		}
	}

	/* Passed all of the above. proper pair, both mapped, in +- orientation. Now check the isize */
	if( abs(bam_alignment_core.isize) < min || abs(bam_alignment_core.isize) > max) // c.a.
	{
		/* Deletion or Insertion */
		return 0;
	}

	/* All passed. Read is concordant */
	return 1;
}

char iupac_base[16]       = {'0','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'};
int  iupac_base_map[26]        = {1,14,2,13,0,0,4,11,0,0,12,0,3,15,0,0,0,5,6,8,0,7,9,0,10,0};
char iupac_complements_map[26] = {'T','V','G','H','X','X','C','D','X','X','M','X','K','N','X','X','X','R','S','A','X','B','W','X','Y','X'};

/* Decode 4-bit encoded bases to their corresponding characters */
char base_as_char( int base_as_int)
{
	return iupac_base[base_as_int];
}

int char_as_base( char base){
	if(base >= 'a' && base <='z')
		base += ('A' - 'a');

	base -= 'A';

	return iupac_base_map[base];
}

/* Return the complement of a base */
char complement_char( char base)
{
	if(base >= 'a' && base <='z')
		base += ('A' - 'a');

	base -= 'A';

	return iupac_complements_map[base];
}

/* Add 33 to the integer value of the qual characters to convert them to ASCII */
void qual_to_ascii( char* qual)
{
	int i;
	for( i = 0; i < strlen( qual); i++)
	{
		qual[i] = qual[i] + 33;
	}
}

/* Even safer than strncpy as it dynamically allocates space for the string if
 there hasn't been already */
void set_str( char** target, char* source)
{
	if( *target != NULL)
	{
		free( ( *target));
	}

	if (source != NULL)
	{
		( *target) = ( char*) malloc( sizeof( char) * ( strlen( source) + 1));
		strncpy( ( *target), source, ( strlen( source) + 1));
	}
	else
	{
		( *target) = NULL;
	}
}


/* Reverse a given string */
void reverse_string( char* str)
{
	int i;
	char swap;
	int len = strlen( str);

	for( i = 0; i < len / 2; i++)
	{
		swap = str[i];
		str[i] = str[len - i - 1];
		str[len - i - 1] = swap;
	}
}

int compare_size_int( const void* p, const void* q)
{
	int i = *( const int*) p;
	int j = *( const int*) q;

	if( i < j)
	{
		return -1;
	}
	else if( i == j)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

void* getMem( size_t size)
{
	void* ret;

	ret = malloc( size);
	if( ret == NULL)
	{
		fprintf( stderr, "Cannot allocate memory. Currently addressed memory = %0.2f MB, requested memory = %0.2f MB.\nCheck the available main memory.\n", getMemUsage(), ( float) ( size / 1048576.0));
		exit( 0);
	}

	memUsage = memUsage + size;
	return ret;
}

void freeMem( void* ptr, size_t size)
{
	memUsage = memUsage - size;
	free( ptr);
}

double getMemUsage()
{
	return memUsage / 1048576.0;
}

void del_char(char *ref, int start, int len){
	int ref_len = strlen(ref);
	int i;
	for (i=start;i<start+len && i<ref_len;i++)
		ref[i]=ref[i+len];
	while (i < ref_len - len){
		ref[i]=ref[i+len];
		i++;
	}
	ref[i]=0;
}

void ins_char(char *ref, char *read, int start_origin, int start_dest, int len){
	int ref_len = strlen(ref);
	int i;

	for (i=ref_len+len; i>=start_dest; i--)
		ref[i]=ref[i-len];

	ref[ref_len+len]=0;
	ref_len = ref_len+len;

	for (i=0; i<len && i+start_dest<ref_len; i++){
		//printf("%c", read[i + start_origin]);
		ref[i + start_dest] = read[i + start_origin]; //'.';
	}
	//printf("\n");
}

pid_t proc_find(const char* name)
{
    DIR* dir;
    struct dirent* ent;
    char* endptr;
    char buf[512];

    if (!(dir = opendir("/proc"))) {
        perror("can't open /proc");
        return -1;
    }

    while((ent = readdir(dir)) != NULL) {
        /* if endptr is not a null character, the directory is not
         * entirely numeric, so ignore it */
        long lpid = strtol(ent->d_name, &endptr, 10);
        if (*endptr != '\0') {
            continue;
        }

        /* try to open the cmdline file */
        snprintf(buf, sizeof(buf), "/proc/%ld/cmdline", lpid);
        FILE* fp = fopen(buf, "r");

        if (fp) {
            if (fgets(buf, sizeof(buf), fp) != NULL) {
                /* check the first token in the file, the program name */
                char* first = strtok(buf, " ");
                if (!strcmp(first, name)) {
                    fclose(fp);
                    closedir(dir);
                    return (pid_t)lpid;
                }
            }
            fclose(fp);
        }

    }

    closedir(dir);
    return -1;
}

void init_sha256_block(BYTE** block) {
	int i;
	/*if( *block != NULL)
	{
		free( ( *block));
	}*/

	(*block) = (BYTE*) malloc(sizeof(BYTE) * SHA256_DIGEST_LENGTH);

	for(i = 0; i<SHA256_DIGEST_LENGTH; i++) {
		(*block)[i] = 0x0;
	}
}

void sha256_hash(char * str, BYTE **ongoing)
{
	BYTE hash[SHA256_DIGEST_LENGTH];
	sha256(str, strlen(str), hash);
	int k;
	for( k = 0; k < SHA256_DIGEST_LENGTH; k++){
		(*ongoing)[k] += hash[k];
	}
}

void apply_cigar_md(char *ref, char *read, char *md, int n_cigar, const uint32_t* _cigar){
	/* Applies given CIGAR and MD operations on the read to match reference */
	int i,j,k;
	char buf[1000];
	int oplen;
	int refptr;
	int edit_loc;
	int carry_pos;
	int delcnt;
	int thisdel;
	int inserted;
	int last_i_cigar;
	operation soft_clips[2];
	int soft_clips_len = 0;
	int cigarlen[n_cigar];
	char cigar[n_cigar];

	for(i=0; i<n_cigar; i++){
		cigar[i] = bam_cigar_opchr(_cigar[i]);
		cigarlen[i] = bam_cigar_oplen(_cigar[i]);
	}

	edit_loc = 0;
	carry_pos = 0;
	for (i=0; i<n_cigar; i++){
		if (cigar[i] == 'M'){
			edit_loc += cigarlen[i];
			carry_pos += cigarlen[i];
		}
		else if (cigar[i] == 'D'){
			del_char(ref, carry_pos, cigarlen[i]);
		}
		else if (cigar[i] == 'I'){
			ins_char(ref, read, edit_loc, carry_pos, cigarlen[i]);
			edit_loc += cigarlen[i];
			carry_pos += cigarlen[i];
		}
		else if (cigar[i] == 'S'){
			soft_clips[soft_clips_len].start = edit_loc;
			soft_clips[soft_clips_len].length = cigarlen[i];
			soft_clips_len++;

			edit_loc += cigarlen[i];
		}
	}
	// Soft clipping, edit read
	for(i=soft_clips_len - 1; i>=0; i--){
		del_char(read, soft_clips[i].start, soft_clips[i].length);
		int ref_len = strlen(ref);
		del_char(ref, ref_len - (soft_clips[i].length), soft_clips[i].length);
	}

	j=0;
	buf[0]=0;
	i=0;
	delcnt=1;
	inserted=0;
	last_i_cigar = 0;
	refptr = 0;

	while (i<strlen(md)){
		if (isdigit(md[i])){
			buf[j++]=md[i];
		}
		else {
			buf[j]=0;
			j=0;
			oplen=atoi(buf);
			refptr+=oplen;
		}

		if (md[i]=='^'){ // del. skip
			thisdel = 0;
			for (k=0; k<n_cigar; k++){
				if (cigar[k] == 'D'){
					thisdel++;
					if (thisdel == delcnt){
						i+=cigarlen[k];
						delcnt++; break;
					}
				}
			}
		}
		else if (isalpha(md[i])){
			inserted = 0;
			edit_loc = 0;

			for (k=0; k<n_cigar; k++){
				if (cigar[k] == 'M'){
					edit_loc += cigarlen[k];
				}
				if (cigar[k] == 'I' && edit_loc <= refptr){
					inserted += cigarlen[k];
				}

			}
			read[refptr+inserted]=md[i];
			refptr++;
		}
		i++;
	}
}

char * get_datetime(){
	time_t     now;
    struct tm*  ts;
    void * buf = getMem(200);
    // Get current time
    time(&now);
    ts = localtime(&now);

    strftime(buf, 200, "%Y-%m-%d %H:%M:%S %Z", ts);
    return buf;
}