#!/bin/bash
./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 1 10 -s r;
./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 1 10 -s t;
./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 1 10 -s t a;
./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 1 10 -s b;
./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 1 10 -s i;

./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 1 1 -s r;
./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 1 1 -s t;
./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 1 1 -s t a;
./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 1 1 -s b;
./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 1 1 -s i;

./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 10 1 -s r;
./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 10 1 -s t;
./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 10 1 -s t a;
./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 10 1 -s b;
./runCellCountCases.py -T 2000 -t 10 -n 50 --dim 120 120 1 1 -i 180 --perm i 0 0 0 -m 10 1 -s i;
