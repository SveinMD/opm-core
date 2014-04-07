#!/bin/bash
nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 1 10 -s r;
nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 1 10 -s t;
nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 1 10 -s t a;
nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 1 10 -s b;
nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 1 10 -s i;

nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 1 1 -s r;
nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 1 1 -s t;
nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 1 1 -s t a;
nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 1 1 -s b;
nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 1 1 -s i;

nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 10 1 -s r;
nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 10 1 -s t;
nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 10 1 -s t a;
nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 10 1 -s b;
nice ./runTimeStepCases.py -T 2000 -n 50 --dim 120 120 20 20 -i 180 --perm i 0 0 0 -m 10 1 -s i;
