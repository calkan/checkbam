#include "fastatools.h"

int readSingleFasta(FILE *fastaFile){
  int cnt;
  char ch; 
  int i,j;
  int index;
  int seqcnt=0, seqlen=0;
  int maxlen;
  char dummy[SEQ_LENGTH];
  char str[SEQ_LENGTH];

  cnt = 0; i=0;
  fprintf(stderr, "Counting sequences.\n");
  maxlen=0;
  rewind(fastaFile);
  while (fscanf(fastaFile, "%c", &ch) > 0){
    if (ch == '>'){
      if (seqlen>maxlen)
	maxlen=seqlen;
      cnt++;
      seqlen=0;
      my_fgets(dummy, SEQ_LENGTH, fastaFile);
    }
    else if (!isspace(ch))
      seqlen++;
  }

  seqcnt = cnt;
  if (seqlen>maxlen)
    maxlen=seqlen;

  //length = (int *) malloc((seqcnt) * sizeof(int));


  cnt = 0; i=0;

  fprintf(stderr, "Allocating memory for %d sequences with max length %d.\n", seqcnt, maxlen);
  
  
  seqs = (char **) malloc((seqcnt) * sizeof(char *));
  
  for (i=0; i<seqcnt; i++)
    seqs[i] = (char *) malloc(maxlen);
  

  names = (char **) malloc((seqcnt) * sizeof(char *));

  for (i=0; i<seqcnt; i++)
    names[i] = (char *) malloc(SEQ_LENGTH);
  
  
  
  for (i=0; i<seqcnt; i++){
    seqs[i][0] = 0;
    names[i][0] = 0;
  }
  

  fprintf(stderr, "Reading sequences.\n");
  rewind(fastaFile);
  cnt = -1; 
  while (fscanf(fastaFile, "%c", &ch) > 0){
    if (ch == '>'){
      cnt++;
      fprintf(stderr, "\r%d\tof\t%d", (cnt+1), seqcnt );
      my_fgets(names[cnt], SEQ_LENGTH, fastaFile);
      names[cnt][strlen(names[cnt])-1] = 0;
    }
    i = 0;
    if (cnt != 0){
      seqs[cnt][i++] = toupper(ch);
    }
    do{
      if (!(fscanf(fastaFile, "%c", &ch) > 0))
	break;
      if (ch!='>' && ch!='\r' && ch!='\n'){
	seqs[cnt][i++] = toupper(ch);
      }
    } while (ch != '>');

    seqs[cnt][i] = 0;
    
    if (ch == '>'){
      cnt++;
      fprintf(stderr, "\r%d\tof\t%d", (cnt+1), seqcnt );
      if (cnt != seqcnt){
	  my_fgets(names[cnt], SEQ_LENGTH, fastaFile);
	  names[cnt][strlen(names[cnt])-1] = 0;
      }
    } // if
  } // while
	    

  fprintf(stderr, "\n[OK] %d sequences read from fasta file.\n",seqcnt);


  return seqcnt;

}


void my_fgets(char *str, int length, FILE *in){
  char ch;
  int i=0;
  while (i<length && fscanf(in, "%c", &ch)){
    if ((ch==' '  || ch=='\t')  && i!=0){
      if (str[i-1] == ' ' || str[i] == '\t')
	;
      else
	str[i++] = toupper(ch);
    }
    else
      str[i++] = toupper(ch);
    if (ch == '\n' || ch=='\r')
      break;
  }
  str[i] = 0;
}

