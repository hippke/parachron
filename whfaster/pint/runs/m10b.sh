#!/bin/bash
cd /home/claude/rebound/whfaster/pint
./pint_long runs/k4 1000 10 360000 2 4 0 0 0 2 1 0 > runs/k4.txt 2>&1
./pint_long runs/h3 1000 10 360000 2 3 0 0 0 2 1 1 > runs/h3.txt 2>&1
./pint_long runs/k5 1000 10 360000 2 5 0 0 0 2 1 0 > runs/k5.txt 2>&1
