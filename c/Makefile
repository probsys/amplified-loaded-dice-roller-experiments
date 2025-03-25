FILES.h := $(wildcard *.h)
FILES.c := $(wildcard $(FILES.h:.h=.c))
FILES.o := $(FILES.c:.c=.o)

LIBS = -lgsl -lgmp -lm

# CFLAGS = '-pg'
CFLAGS ?= -O3 -DNDEBUG -Wno-unused-result -flto -march=native

# .PRECIOUS: %.o
%.o: %.c %.h
	gcc -c $(CFLAGS) -o $@ $<

%.s: %.c ${FILES.o}
	gcc $(CFLAGS) -S -o $@ $^ $(LIBS)

%.out: %.c ${FILES.o}
	gcc $(CFLAGS) -o $@ $^ $(LIBS)

%.valgrind: %.out
	valgrind --leak-check=full \
			 --show-leak-kinds=all \
			 --track-origins=yes \
			 --verbose \
			 --log-file=$@ \
			 ./$^

.PHONY: clean
clean:
	rm -rf \
		*.a *.o *.s *.out *.debug *.debug __pycache__ \
		*.gcno *.gcov *.gcda *.valgrind

all: main.out
