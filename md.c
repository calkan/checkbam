#include <stdio.h>
#include <string.h>


void del_char(char *ref, int start, int len){
	int ref_len = strlen(ref);
	int i;
	printf("deleting ");
	for (i=start;i<start+len && i<ref_len;i++){
		printf ("%c", ref[i]);
		ref[i]=ref[i+len];
	}
	while (i < ref_len - len){
		ref[i]=ref[i+len];
		i++;
	}
	ref[i]=0;
	printf("\n");
}

void ins_char(char *ref, char *read, int start_origin, int start_dest, int len){
	int ref_len = strlen(ref);
	int i;

	printf("inserting position %d to %d ", start_dest, start_origin);
	for (i=ref_len+len; i>=start_dest; i--)
		ref[i]=ref[i-len];

	ref[ref_len+len]=0;
	ref_len = ref_len+len;

	for (i=0; i<len && i+start_dest<ref_len; i++){	 
		printf("%c", read[i + start_origin]);
		ref[i + start_dest] = read[i + start_origin]; //'.';
	}
	printf("\n");
}

struct operation_s {
	unsigned short int start;
	unsigned short int length;
} init_operation = {0,0};

typedef struct operation_s operation;

void apply_cigar_md(char *ref, char *read, char *md, int n_cigar, char *cigar, int *cigarlen){
	int i,j,k;
	char buf[1000];
	int oplen;
	int refptr;
	int edit_loc;
	int carry_pos;
	int delcnt;
	int thisdel;
	int inserted;
	int skipk;
	operation soft_clips[2] = {init_operation};
	int soft_clips_len = 0;

	edit_loc = 0;
	carry_pos = 0;
	for (i=0; i<n_cigar; i++){
		if (cigar[i] == 'M'){
			edit_loc += cigarlen[i];
			carry_pos += cigarlen[i];
		}
		else if (cigar[i] == 'D'){
			del_char(ref, carry_pos, cigarlen[i]);
			//edit_loc -= cigarlen[i];
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
		//fprintf(stdout, "%d\t%c\t", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(c igar[i]));
		//fprintf(stdout, "%d\t%d\t%c\t%d\t", bam_cigar_op(cigar[i]), bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]), bam_cigar_type(cigar[i]));
	}
	// Soft clipping, edit read
	for(i=soft_clips_len - 1; i>=0; i--){
		printf("Deleting read from %d for %d\n", soft_clips[i].start, soft_clips[i].length);
		del_char(read, soft_clips[i].start, soft_clips[i].length);
		int ref_len = strlen(ref);
		printf("Deleting ref from %d for %d\n", ref_len - (soft_clips[i].length), soft_clips[i].length);
		del_char(ref, ref_len - (soft_clips[i].length), soft_clips[i].length);
	}

	char index[100] = "'         '         '         '         '         '         '         '         '         '         ";
	index[100] = '\0';
	fprintf(stdout, "\ndel\n%s\n%s\n%s\n", index, read, ref);

	j=0;
	buf[0]=0;
	i=0;
	delcnt=1;
	inserted=0;
	skipk = 0;
	refptr = 0;
	//printf("md: %s\n", md);
	while (i<strlen(md)){
		//printf("md[%d]: %c\n", i, md[i]);
		if (isdigit(md[i])){
			//printf("digit %c\n", md[i]);
			buf[j++]=md[i];
		}
		else {
			buf[j]=0;
			j=0;
			oplen=atoi(buf);
			refptr+=oplen;
			printf("buf: %s : %d - %d\n", buf, oplen, refptr);
		}

		if (md[i]=='^'){ // del. skip
			thisdel = 0;
			for (k=0; k<n_cigar; k++){
				if (cigar[k] == 'D'){
					thisdel++;
					printf("thisdel %d	delcnt %d\n", thisdel, delcnt);
					if (thisdel == delcnt){
						printf("md i %d ", i);
						i+=cigarlen[k]; 
						printf("now	%d.	md[i]=%c\n", i, md[i]);			
						delcnt++; break;
					}
				}
			}
		}
		else if (isalpha(md[i])){
			inserted = 0;
			edit_loc = 0;
			
			for (k=skipk; k<n_cigar; k++){
				printf("editloc %d	refptr %d\n", edit_loc, refptr);
				if (cigar[k] == 'I' && edit_loc <= refptr){
					// todo: until refptr
					inserted += cigarlen[k];
					skipk = k+1;
				}

				if (cigar[k] != 'S'){
					edit_loc += cigarlen[k];
				}
			}
			refptr += inserted;
			printf("refptr added %d\n", inserted);
			printf("changing %d:%c %d:%c\n", refptr, read[refptr], (i+1), md[i]);
			read[refptr++]=md[i];
			printf("postedit refptr %d\n", refptr);
			fprintf(stdout, "\npostedit\n%s\n%s\n%s\n", index, read, ref);
		}
		i++; 
		printf("lnow	%d.	md[i]=%c\n", i, md[i]);			
	}
}


int main(){
	char read[1000], ref[1000];

	char cigar[100]; 
	int n_cigar; 
	int cigarlen[100];
	char md[1000];
	int i;
	int testcnt=1;
	while (fscanf(stdin, "%d", &n_cigar) > 0){
		if (feof(stdin)) break;
		printf("Testing %d\n", testcnt++);
		for (i=0;i<n_cigar;i++){
			printf("%d. cigar operation", i);
			fscanf(stdin, "\t%d\t%c", &(cigarlen[i]), &(cigar[i]));
		}
		
		fscanf(stdin, "%s", md);
		fscanf(stdin, "%s", read);
		fscanf(stdin, "%s", ref);
	
		apply_cigar_md(ref, read, md+1, n_cigar, cigar, cigarlen);
		
		if (strcmp(read, ref)){
			fprintf(stdout, "\npos\n%s\n%s\n", read, ref);
			return;
		}
		else{	
			printf("ok\n") ;
		}
	}
}
