CC=gcc
CFLAGS=-I.
DEPS = quad_tree.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: quad_tree.o nbody.o
	$(CC) -o nbody nbody.o quad_tree.c -lm

clean:
	\rm -f *.o all *~ *# 