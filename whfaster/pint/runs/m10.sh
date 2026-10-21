#!/bin/bash
cd /home/claude/rebound/whfaster/pint
./pint_long runs/m10 1000 10 360000 2 3 1 0 0 2 1 > runs/m10_a.txt 2>&1
cp runs/m10_pr.bin runs/m10_prproj.bin
./pint_long runs/m10s1 1000 10 360000 1 0 1 0 0 1 0 > /dev/null 2>&1 || true
