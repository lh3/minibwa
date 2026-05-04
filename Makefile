CC=			gcc
CFLAGS=		-std=c99 -g -Wall -O3
CPPFLAGS=
LDFLAGS=
INCLUDES=
LOBJS=		kommon.o kalloc.o bwt.o l2bit.o options.o seed.o map-algo.o lchain.o align.o pe.o cs.o format.o meth.o \
			ksw2_extz2_sse.o ksw2_extd2_sse.o ksw2_ll_sse.o
AOBJS=		kthread.o libsais.o libsais64.o index.o bseq.o map-main.o fastmap.o
MALLOC_O=	mimalloc.o
PROG=		minibwa
LIBS=		-lpthread -lz -lm
ARCH=		$(shell uname -m)
omp=		$(shell printf '\043include <omp.h>\nint main(){return 0;}' | $(CC) -x c -fopenmp - 2>/dev/null && echo "1" || echo "0")

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LDFLAGS+=-fsanitize=address
	LIBS+=-ldl
endif

ifeq ($(omp),1)
	CPPFLAGS+=-DLIBSAIS_OPENMP
	CFLAGS+=-fopenmp
	LIBS+=-fopenmp
endif

ifneq ($(gpl),0)
	AOBJS+=QSufSort.o bwtgen.o
	CPPFLAGS+=-DUSE_GPL
endif

ifeq ($(mimalloc),0)
	MALLOC_O=
	CPPFLAGS+=-DHAVE_KALLOC
endif

ifeq ($(ARCH), x86_64)
	CFLAGS+=-msse4.2 -mpopcnt
endif

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

mimalloc.o:
		$(CC) -c -std=gnu99 -O3 -Wall -Wextra -DNDEBUG -DMI_MALLOC_OVERRIDE -DMI_OSX_INTERPOSE=1 -DMI_OSX_ZONE=1 -Imimalloc mimalloc/static.c -o $@

libminibwa.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

minibwa:libminibwa.a $(MALLOC_O) $(AOBJS) main.o
		$(CC) $(CFLAGS) $(LDFLAGS) $(MALLOC_O) $(AOBJS) main.o -o $@ -L. -lminibwa $(LIBS)

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cpp)

# DO NOT DELETE

QSufSort.o: QSufSort.h
align.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h kalloc.h ksw2.h
bseq.o: bseq.h kommon.h kseq.h
bwt.o: kommon.h kalloc.h bwt.h
bwtgen.o: QSufSort.h
cs.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h kalloc.h
fastmap.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h ketopt.h kseq.h
fastmap.o: kalloc.h
format.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h meth.h
index.o: libsais64.h kommon.h ketopt.h mbpriv.h minibwa.h l2bit.h bwt.h
index.o: bseq.h kseq.h
meth.o: meth.h l2bit.h
kalloc.o: kalloc.h
kommon.o: kommon.h
ksw2_extd2_sse.o: ksw2.h
ksw2_extz2_sse.o: ksw2.h
ksw2_ll_sse.o: ksw2.h
kthread.o: kthread.h
l2bit.o: kommon.h l2bit.h kseq.h
lchain.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h kalloc.h ksort.h
libsais.o: libsais.h
libsais64.o: libsais.h libsais64.h
main.o: kommon.h mbpriv.h minibwa.h l2bit.h bwt.h bseq.h ketopt.h kseq.h
map-algo.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h kalloc.h ksort.h
map-main.o: kommon.h mbpriv.h minibwa.h l2bit.h bwt.h bseq.h kalloc.h
map-main.o: kthread.h ketopt.h
options.o: minibwa.h
pe.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h kalloc.h ksw2.h
seed.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h kalloc.h ksort.h
