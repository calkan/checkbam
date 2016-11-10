#ifndef __TARDIS
#define __TARDIS


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include <signal.h>
#include <time.h>
#include <unistd.h>

// For unix sockets while using daemon
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>

#include "common.h"
#include "cmdline.h"
#include "processbam.h"

#define SOCK_PATH "verifybam.socket"

void init_server(parameters **params);
void init_client(parameters *params);

#endif
