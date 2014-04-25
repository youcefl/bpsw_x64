

CXX = g++
CXXFLAGS = -std=c++0x -Ofast
NASM_FLAGS = -D WIN64 -f win64

BIN = ../bin

all: $(BIN)/asm

OBJS =  $(BIN)/arith.o $(BIN)/lucas.o $(BIN)/main.o

$(BIN)/asm: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(BIN)/asm $(OBJS)

$(BIN)/lucas.o: lucas.s
	nasm $(NASM_FLAGS) -o $(BIN)/lucas.o lucas.s

$(BIN)/arith.o: arith.s
	nasm $(NASM_FLAGS) -o $(BIN)/arith.o arith.s

$(BIN)/main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o $(BIN)/main.o main.cpp

clean:
	rm -f $(BIN)/asm.exe
	rm -f $(BIN)/*.o

