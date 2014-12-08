nvcc -O -o td -lcuda testdriver.c
./td -bin k6.bin -block 8 512 -thread 32 1 -mat 1024 -size 1024 -check
