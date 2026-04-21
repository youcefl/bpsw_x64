# A primality test written in assembly x64

This is an implementation of the Baille-PSW primality test in x86_64 assembly.
It was written for fun, there is no associated claim regarding the performance.<br>
I put it here in case someone finds it useful.
I wrote this back in 2014. It still compiles and runs in 2026 without any
modification to the core algorithm (only added a `-c` flag for benchmarking).

## Building the code

**Platform support:** Windows only (MinGW-w64). Visual Studio is not supported.<br>
To build the code source msys2 is needed, the installer can be downloaded from [msys2.org](https://www.msys2.org/).<br>
Once msys2 is installed run MINGW64, install the needed packages:
```bash
$ pacman -S mingw-w64-x86_64-gcc
$ pacman -S mingw-w64-x86_64-nasm
$ pacman -S make
```

To run the build simply do
```bash
$ cd src
$ mkdir -p ../bin
$ make
nasm -D WIN64 -f win64 -o ../bin/arith.o arith.s
g++ -static -std=c++0x -Ofast -c -o ../bin/main.o main.cpp
g++ -static -std=c++0x -Ofast -o ../bin/asm ../bin/arith.o ../bin/main.o
```

## Messing with the primes
### Tests
Upon successfull build one -has to- can run the tests:
```bash
$ ../bin/asm -t
jacobi_symbol(0, 1)                                        OK
jacobi_symbol(0, 2)                                        OK
# ...snip...
is_slprp(17000000000000000023, 5) = 0                      OK
is_slprp(18446743979220271189, -7) = 0                     OK
is_slprp(18446744073709551583, 5) = 0                      OK
is_prime(0)                                                OK
is_prime(1)                                                OK
is_prime(2)                                                OK
# ...snip...
is_prime(18446744073709551613)                             OK
is_prime(18446744073709551615)                             OK
ALL 101 tests passed.
```
### Usage
If all the tests are OK one can dig deeper:
```bash
# To output all primes below L use 'asm -p L [-c]'
$ ../bin/asm -p 100

          2          3          5          7         11         13         17         19

         23         29         31         37         41         43         47         53

         59         61         67         71         73         79         83         89

         97

# Add '-c' for counting only
$ ../bin/asm -p 100 -c
Seconds: 0.0005111
Primes: 25

# To output the primes in range [S, S + L[ use 'asm -s S -l L [-c]'
# e.g. primes in [2^61-1, 2^61+99[
$ ../bin/asm -s 2305843009213693951 -l 100
2305843009213693951
2305843009213693967
2305843009213693973
2305843009213694009
2305843009213694017

# Add '-c' for counting only
$ ../bin/asm -s 2305843009213693951 -l 100 -c
Seconds: 0.0010093
Primes: 5

```

## Benchmark

This section was added solely in order to give an idea of the kind of performance one can expect from the test.

### Environment
**Processor:**	Intel(R) Core(TM) Ultra 7 155H, 1400 Mhz, (16 cores, 22 threads)<br>
**RAM:** 64.0 GB<br>
**OS:** Microsoft Windows 11 Pro<br>
**Compiler:** gcc (Rev13, Built by MSYS2 project) 15.2.0 (mingw-w64)<br>
**Assembler:** nasm version 3.01


### Results

|Range start|Range length|Primes in range|Time(s)|Throughput (tests/s)|
|-----------|------------|---------------|-------|--------------------|
|0|10^8|5'761'455|6.6|15'151'515|
|10^9|10^8|4'814'936|7|14'285'714|
|10^15|10^8|2'893'937|8.7|11'494'252|
|10^18|10^8|2'414'886|9.3|10'752'688|

