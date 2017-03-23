CC = gcc

SOURCES = main.c Lab4_IO.c datatrim.c serialtester.c

all: main datatrim serialtester

main: main.c 
	mpicc -g -Wall -o main main.c Lab4_IO.c -lm

datatrim: datatrim.c
	gcc -g -Wall -o datatrim datatrim.c

serialtester: serialtester.c
	gcc -g -Wall -o serialtester serialtester.c Lab4_IO.c -lm

clean:
	-rm -f *.o main datatrim serialtester

run:
	./main ./datatrim ./serialtester
