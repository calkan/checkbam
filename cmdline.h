#ifndef __COMMANDLINE
#define __COMMANDLINE

#define EXE_VERIFYBAM 1
#define EXE_FQHASH 2

int parse_command_line( int, char**, parameters*, int);
void parse_fastq_list( parameters*);
void print_help( int);

#endif
