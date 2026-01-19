CC=			gcc
CXX=		g++
CFLAGS=		-std=c99 -g -Wall -O3
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=
INCLUDES=
LOBJS=		kommon.o kalloc.o bwt.o l2bit.o options.o seed.o map-algo.o lchain.o align.o format.o \
			ksw2_extz2_sse.o ksw2_extd2_sse.o ksw2_ll_sse.o
AOBJS=		kthread.o QSufSort.o bwtgen.o libsais.o libsais64.o index.o bseq.o map-main.o fastmap.o
PROG=		minibwa
LIBS=		-lpthread -lz -lm
ARCH=		$(shell uname -m)

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

ifneq ($(omp),0)
	CPPFLAGS=-DLIBSAIS_OPENMP
	CFLAGS+=-fopenmp
	LIBS+=-fopenmp
endif

ifeq ($(ARCH), x86_64)
	CFLAGS+=-msse4.2 -mpopcnt
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
align.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h kalloc.h ksw2.h
bseq.o: bseq.h kseq.h
bwt.o: kommon.h kalloc.h bwt.h
bwtgen.o: QSufSort.h
fastmap.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h ketopt.h kseq.h kalloc.h
format.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h
index.o: libsais64.h kommon.h ketopt.h mbpriv.h minibwa.h l2bit.h bwt.h
kalloc.o: kalloc.h
kommon.o: kommon.h
ksw2_extd2_sse.o: ksw2.h kalloc.h
ksw2_extz2_sse.o: ksw2.h kalloc.h
ksw2_ll_sse.o: ksw2.h kalloc.h
kthread.o: kthread.h
l2bit.o: kommon.h l2bit.h kseq.h
lchain.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h kalloc.h ksort.h
libsais.o: libsais.h
libsais64.o: libsais.h libsais64.h
main.o: kommon.h mbpriv.h minibwa.h l2bit.h bwt.h ketopt.h kseq.h
map-algo.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h kalloc.h ksort.h
map-main.o: kommon.h mbpriv.h minibwa.h l2bit.h bwt.h bseq.h kalloc.h
map-main.o: kthread.h ketopt.h
options.o: minibwa.h
seed.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h kalloc.h ksort.h
