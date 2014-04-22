

CC = g++ -std=c++0x -g
NASM_FLAGS = -D WIN64 -g -f win64


all: asm


asm: lucas.o main.o
	$(CC) -o asm lucas.o main.o

lucas.o: lucas.s
	nasm $(NASM_FLAGS) -o lucas.o lucas.s

main.o: main.cpp
	$(CC) -c -o main.o main.cpp

clean:
	rm -f asm
	rm -f *.o

