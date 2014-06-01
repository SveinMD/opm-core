#!/bin/bash
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 1 10 -s r;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 1 10 -s t;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 1 10 -s a;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 1 10 -s b;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 1 10 -s i;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 1 10 -s g;

nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 1 1 -s r;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 1 1 -s t;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 1 1 -s a;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 1 1 -s b;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 1 1 -s i;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 1 1 -s g;

nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 10 1 -s r;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 10 1 -s t;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 10 1 -s a;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 10 1 -s b;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 10 1 -s i;
nice ./runTimeStepCases_short.py -f ../bin/spe10 -T 30 -n 20 --dim 60 60 10 60 220 1 -i 50 --perm i 36 0 0 -m 10 1 -s g;
