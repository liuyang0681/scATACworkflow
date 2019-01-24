PROG= parse_alignment_readname sc_atac_parse_cellbarcode parse_Tn5me_adaptors
DEBUG_PROG= 

all: $(PROG)

debug: $(DEBUG_PROG)

HTSDIR = htslib-1.9
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a

CC       = gcc
CFLAGS   = -Wall -O2
DEBUG_CFLAGS   = -g -Wall -O0 -DDEBUG_MODE
DFLAGS   =
INCLUDES = -I src/ -I$(HTSDIR)/
LIBS = -lz

all:$(PROG)

# See htslib/Makefile
PACKAGE_VERSION := $(shell git describe --tags)

version.h:
	echo '#define BCFANNO_VERSION "$(PACKAGE_VERSION)"' > $@


.SUFFIXES:.c .o
.PHONY:all clean clean-all clean-plugins distclean install lib tags test testclean force plugins docs

force:

makedir:
	-mkdir -p bin

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

parse_alignment_readname: $(HTSLIB) makedir
	$(CC) $(INCLUDES) -pthread -o bin/$@ src/number.c src/thread_pool.c src/parse_alignment_readname.c $(HTSLIB) $(LIBS)

sc_atac_parse_cellbarcode: $(HTSLIB) makedir
	$(CC) $(INCLUDES) -pthread -o bin/$@ src/number.c src/thread_pool.c src/fastq.c src/sc_atac_parse_cellbarcode.c $(HTSLIB) $(LIBS)

parse_Tn5me_adaptors: $(HTSLIB) makedir
	$(CC) $(INCLUDES) -pthread -o bin/$@ src/number.c src/thread_pool.c src/fastq.c src/parse_Tn5ME_adaptors.c $(HTSLIB) $(LIBS)

clean: testclean
	-rm -f gmon.out *.o *~ $(PROG) version.h

testclean:
	-rm -f test/*.o test/*~ $(TEST_PROG)
	-rm -rf *.dSYM bin

distclean: clean
	-rm -f TAGS

clean-all: clean clean-htslib

tags:
	ctags -f TAGS *.[ch] plugins/*.[ch]
