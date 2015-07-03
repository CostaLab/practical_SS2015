#!/bin/bash

gcc -Wall -fPIC -c libdistint.c -std=c11 -O3
gcc -shared -Wl,-soname,libdistint.so.1 -o libdistint.so.1.0 libdistint.o 
