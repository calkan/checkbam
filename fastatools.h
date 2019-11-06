#ifndef __FASTATOOLS
#define __FASTATOOLS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>


#define SEQ_LENGTH 150

char **seqs;
char **names;

int readSingleFasta(FILE *);
void my_fgets(char *, int, FILE *);

#endif 
