CC=			gcc
CFLAGS=		-std=gnu99 -Wall -O3
EXE=		mbmap-sgl mbmap-batch

.PHONY:all clean

all:$(EXE)

../libminibwa.a:
	(cd ..; make libminibwa.a)

mbmap-sgl:ex-sgl.c ../libminibwa.a
	$(CC) $(CFLAGS) -o $@ -I.. $^ -lz -lm

mbmap-batch:ex-batch.c ../libminibwa.a
	$(CC) $(CFLAGS) -o $@ -I.. $^ -lz -lm

clean:
	rm -f $(EXE)
