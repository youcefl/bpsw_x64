

CC = g++ -std=c++0x -g

all: asm


asm: lucas.o main.o
	$(CC) -o asm lucas.o main.o

lucas.o: lucas.s
	nasm -g -f elf64 -o lucas.o -F dwarf lucas.s

main.o: main.cpp
	$(CC) -c -o main.o main.cpp

clean:
	rm -f asm
	rm -f *.o

