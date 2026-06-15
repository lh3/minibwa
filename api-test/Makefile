CC=			gcc
CFLAGS=		-std=gnu99 -Wall -O3
EXE=		mbmap-one mbmap-batch

# libminibwa.a may contain OpenMP code (parallel BWT construction), so the
# example programs must link the OpenMP runtime. Mirror the parent Makefile's
# detection so the setting matches the library that was built.
UNAME_S=	$(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	LIBOMP_PREFIX?=	$(shell test -d /opt/homebrew/opt/libomp && echo /opt/homebrew/opt/libomp || echo /usr/local/opt/libomp)
	omp_cflags=	-Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
	omp_libs=	-L$(LIBOMP_PREFIX)/lib -Wl,-rpath,$(LIBOMP_PREFIX)/lib -lomp
else
	omp_cflags=	-fopenmp
	omp_libs=	-fopenmp
endif
omp?=		$(shell printf '\043include <omp.h>\nint main(){return omp_get_max_threads();}' | $(CC) -x c $(omp_cflags) - $(omp_libs) -o /dev/null 2>/dev/null && echo "1" || echo "0")
ifeq ($(omp),1)
	LIBS+=$(omp_libs)
endif

.PHONY:all clean

all:$(EXE)

../libminibwa.a:
	(cd ..; make libminibwa.a)

mbmap-one:ex-one.c ../libminibwa.a
	$(CC) $(CFLAGS) -o $@ -I.. $^ -lz -lm $(LIBS)

mbmap-batch:ex-batch.c ../libminibwa.a
	$(CC) $(CFLAGS) -o $@ -I.. $^ -lz -lm $(LIBS)

clean:
	rm -f $(EXE)
