CC=			gcc
CFLAGS=		-std=gnu99 -Wall -O3
EXE=		mbmap-sgl mbmap-batch test-meth-cmap

.PHONY:all clean test

all:$(EXE)

test:test-meth-cmap
	./test-meth-cmap

../libminibwa.a:
	(cd ..; make libminibwa.a)

mbmap-sgl:ex-sgl.c ../libminibwa.a
	$(CC) $(CFLAGS) -o $@ -I.. $^ -lz -lm

mbmap-batch:ex-batch.c ../libminibwa.a
	$(CC) $(CFLAGS) -o $@ -I.. $^ -lz -lm

test-meth-cmap:ex-meth-cmap.c ../libminibwa.a
	$(CC) $(CFLAGS) -o $@ -I.. $^ -lz -lm -lpthread

clean:
	rm -f $(EXE)
