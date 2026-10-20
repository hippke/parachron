#!/bin/sh
# Variants from kernel files: newavx512 origavx512 newavx2 origavx2 (+ library kernels via lab_lib)
set -e
S=/home/claude/rebound/src
gcc -O3 -std=c99 -D_GNU_SOURCE -I$S -DSUFFIX=new512 -c variant_new512.c -o v_new512.o
gcc -O3 -std=c99 -D_GNU_SOURCE -I$S -DSUFFIX=orig512 -c variant_orig512.c -o v_orig512.o
gcc -O3 -std=c99 -D_GNU_SOURCE -I$S -c variant_newavx2.c -o v_newavx2.o
gcc -O3 -std=c99 -D_GNU_SOURCE -I$S -c variant_origavx2.c -o v_origavx2.o
printf "X(orig512)\nX(new512)\nX(origavx2)\nX(newavx2)\n" > variants.def
gcc -O3 -march=native -I$S lab.c v_orig512.o v_new512.o v_origavx2.o v_newavx2.o -o lab2 -L$S -lrebound -Wl,-rpath,$S -lm
