VERIFYBAM_VERSION := "0.2-alpha"
VERIFYBAM_UPDATE := "August 01, 2016"
VERIFYBAM_DEBUG := 1
BUILD_DATE := "$(shell date)"
CC=gcc
CFLAGS =  -O0 -g -I htslib -DVERIFYBAM_VERSION=\"$(VERIFYBAM_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DVERIFYBAM_UPDATE=\"$(VERIFYBAM_UPDATE)\" -DVERIFYBAM_DEBUG=$(VERIFYBAM_DEBUG)
LIBRARY_PATH = /usr/local/lib/libhts.a
SOURCES = verifybam.c cmdline.c common.c processbam.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = verifybam
INSTALLPATH = /usr/local/bin/

ifneq ("$(wildcard $(LIBRARY_PATH))","")
    LDFLAGS = /usr/local/lib/libhts.a -lz -lm -lpthread
    LIB_EXISTS = 1
else
    LDFLAGS = htslib/libhts.a -lz -lm -lpthread
    LIB_EXISTS = 0
endif

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS)
ifneq ($(LIB_EXISTS), 1)
	make libs
endif
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECUTABLE) *.o *~ \#*

libs:
	make -C htslib

install:
	cp verifybam $(INSTALLPATH)
