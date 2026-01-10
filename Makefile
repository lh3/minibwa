CC=			gcc
CXX=		g++
CFLAGS=		-std=c99 -g -Wall -O3
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=
INCLUDES=
LOBJS=		kommon.o kalloc.o bwt.o l2bit.o options.o map-algo.o
AOBJS=		kthread.o QSufSort.o bwtgen.o libsais.o libsais64.o index.o bseq.o map-main.o
PROG=		minibwa
LIBS=		-lpthread -lz -lm

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

ifneq ($(omp),0)
	CPPFLAGS=-DLIBSAIS_OPENMP
	CFLAGS+=-fopenmp
	LIBS+=-fopenmp
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

libminibwa.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

minibwa:libminibwa.a $(AOBJS) main.o
		$(CC) $(CFLAGS) $(AOBJS) main.o -o $@ -L. -lminibwa $(LIBS)

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cpp)

# DO NOT DELETE

QSufSort.o: QSufSort.h
bseq.o: bseq.h kseq.h
bwt.o: kommon.h kalloc.h bwt.h
bwtgen.o: QSufSort.h
index.o: libsais64.h kommon.h ketopt.h mbpriv.h minibwa.h l2bit.h bwt.h
kalloc.o: kalloc.h
kommon.o: kommon.h
kthread.o: kthread.h
l2bit.o: kommon.h l2bit.h kseq.h
libsais.o: libsais.h
libsais64.o: libsais.h libsais64.h
main.o: kommon.h mbpriv.h minibwa.h l2bit.h bwt.h ketopt.h kseq.h
map-algo.o: mbpriv.h minibwa.h l2bit.h bwt.h kalloc.h kommon.h
map-main.o: kommon.h mbpriv.h minibwa.h l2bit.h bwt.h bseq.h kalloc.h
map-main.o: kthread.h ketopt.h
options.o: minibwa.h
