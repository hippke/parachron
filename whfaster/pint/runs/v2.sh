#!/bin/bash
cd /home/claude/rebound/whfaster/pint
./pint_long runs/v2 1000 10 360000 2 3 1 0 0 2 1 1 1e-6 > runs/v2.txt 2>&1
./pint_long runs/v2a 1000 10 360000 1 0 1 0 0 1 0 1 1e-6 > /dev/null 2>&1
./ser_pert runs/v2p.bin 10000 360000 1.000000001e-6
