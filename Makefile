CC=			gcc
CXX=		g++
CFLAGS=		-std=c99 -g -Wall -O3
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=
INCLUDES=
OBJS=		hl.o utils.o kalloc.o sys.o bwt.o l2bit.o QSufSort.o bwtgen.o
PROG=		minibwa
LIBS=		-lpthread -lz -lm

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

libminibwa.a:$(OBJS)
		$(AR) -csru $@ $(OBJS)

minibwa:main.o libminibwa.a
		$(CC) $(CFLAGS) $< -o $@ -L. -lminibwa $(LIBS)

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cpp)

# DO NOT DELETE

QSufSort.o: QSufSort.h
bseq.o: bseq.h kseq.h
bwt.o: hl.h kalloc.h bwt.h
bwtgen.o: QSufSort.h
hl.o: hl.h
index.o: bwt.h
kalloc.o: kalloc.h
l2bit.o: hl.h l2bit.h kseq.h
main.o: bwt.h sys.h utils.h ketopt.h l2bit.h hl.h
msais.o: msais.h
sys.o: sys.h
utils.o: utils.h
