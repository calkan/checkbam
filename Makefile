VERIFYBAM_VERSION := "0.0.3"
VERIFYBAM_UPDATE := "January 1, 2018"
VERIFYBAM_DEBUG := 1

FQHASH_VERSION := "0.0.1"
FQHASH_UPDATE := "September 2, 2016"
FQHASH_DEBUG := 0

BUILD_DATE := "$(shell date)"
CC=gcc
CFLAGS_VERIFYBAM =  -lpthread -O0 -g -I htslib -DVERSION=\"$(VERIFYBAM_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DUPDATE=\"$(VERIFYBAM_UPDATE)\" -DDEBUG=$(VERIFYBAM_DEBUG)
CFLAGS_FQHASH =  -lpthread -O0 -g -I htslib -DVERSION=\"$(FQHASH_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DUPDATE=\"$(FQHASH_UPDATE)\" -DDEBUG=$(FQHASH_DEBUG)
SOURCES_VERIFYBAM = verifybam.c cmdline.c common.c processbam.c md5.c
SOURCES_FQHASH = cmdline.c common.c md5.c fqhash.c
OBJECTS_VERIFYBAM = $(SOURCES_VERIFYBAM:.c=.o)
OBJECTS_FQHASH = $(SOURCES_FQHASH:.c=.o)
EXE_VERIFYBAM = verifybam
EXE_FQHASH = fqhash
INSTALLPATH = /usr/local/bin/

LIBRARY_PATH = /usr/local/lib/libhts.a
ifneq ("$(wildcard $(LIBRARY_PATH))","")
    LDFLAGS = /usr/local/lib/libhts.a -lz -lm -lpthread
    LIB_EXISTS = 1
else
    LDFLAGS = htslib/libhts.a -lz -lm -lpthread
    LIB_EXISTS = 0
endif

all: $(SOURCES_VERIFYBAM) $(SOURCES_FQHASH) $(EXE_VERIFYBAM) $(EXE_FQHASH)
	rm -rf *.o

$(EXE_VERIFYBAM): CFLAGS = $(CFLAGS_VERIFYBAM)
$(EXE_VERIFYBAM): $(OBJECTS_VERIFYBAM)
ifneq ($(LIB_EXISTS), 1)
	make libs
endif
	$(CC) $(OBJECTS_VERIFYBAM) -o $@ $(LDFLAGS)

$(EXE_FQHASH): CFLAGS = $(CFLAGS_FQHASH)
$(EXE_FQHASH): $(OBJECTS_FQHASH)
ifneq ($(LIB_EXISTS), 1)
	make libs
endif
	$(CC) $(OBJECTS_FQHASH) -o $@ $(LDFLAGS)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXE_VERIFYBAM) $(EXE_FQHASH) *.o *~ \#*

libs:
	make -C htslib

install:
	cp verifybam $(INSTALLPATH)
	cp fqhash $(INSTALLPATH)
