#!/bin/sh
# usage: ./build.sh "name:FLAGS" ...   (first is the reference; generic variants: name must start with g_)
set -e
S=/home/claude/rebound/src
: > variants.def; OBJS=""
for spec in "$@"; do
  name=${spec%%:*}; flags=${spec#*:}; [ "$flags" = "$spec" ] && flags=""
  src=variant_avx512.c; case $name in g_*) src=variant_generic.c;; esac
  gcc -O3 -std=c99 -D_GNU_SOURCE -I$S -DSUFFIX=$name $flags -c $src -o v_$name.o
  echo "X($name)" >> variants.def; OBJS="$OBJS v_$name.o"
done
gcc -O3 -march=native -I$S lab.c $OBJS -o lab -L$S -lrebound -Wl,-rpath,$S -lm
