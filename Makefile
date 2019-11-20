VERIFYBAM_VERSION := "0.0.3"
VERIFYBAM_UPDATE := "January 1, 2018"
VERIFYBAM_DEBUG := 1

BUILD_DATE := "$(shell date)"
CC=gcc
CFLAGS_VERIFYBAM =  -lpthread -O3 -g -I htslib -DVERSION=\"$(VERIFYBAM_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DUPDATE=\"$(VERIFYBAM_UPDATE)\" -DDEBUG=$(VERIFYBAM_DEBUG)
SOURCES_VERIFYBAM = verifybam.c cmdline.c common.c processbam.c sha2.c
OBJECTS_VERIFYBAM = $(SOURCES_VERIFYBAM:.c=.o)
EXE_VERIFYBAM = verifybam
INSTALLPATH = /usr/local/bin/

LIBRARY_PATH = /usr/local/lib/libhts.a
ifneq ("$(wildcard $(LIBRARY_PATH))","")
    LDFLAGS = /usr/local/lib/libhts.a -lz -lm -lpthread -lcrypto
    LIB_EXISTS = 1
else
    LDFLAGS = htslib/libhts.a -lz -lm -lpthread -lcrypto
    LIB_EXISTS = 0
endif

all: $(SOURCES_VERIFYBAM) $(EXE_VERIFYBAM)
	rm -rf *.o

$(EXE_VERIFYBAM): CFLAGS = $(CFLAGS_VERIFYBAM)
$(EXE_VERIFYBAM): $(OBJECTS_VERIFYBAM)
ifneq ($(LIB_EXISTS), 1)
	make libs
endif
	$(CC) $(OBJECTS_VERIFYBAM) -o $@ $(LDFLAGS)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXE_VERIFYBAM) *.o *~ \#*

libs:
	make -C htslib

install:
	cp verifybam $(INSTALLPATH)
