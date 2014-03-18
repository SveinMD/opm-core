#!/bin/bash
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 1 -m 1 1 -s r;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 1 -m 1 1 -s t;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 1 -m 1 1 -s t a;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 1 -m 1 1 -s b;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 1 -m 1 1 -s i;

./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 10 -m 1 1 -s r;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 10 -m 1 1 -s t;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 10 -m 1 1 -s t a;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 10 -m 1 1 -s b;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 10 -m 1 1 -s i;

./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 100 -m 1 1 -s r;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 100 -m 1 1 -s t;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 100 -m 1 1 -s t a;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 100 -m 1 1 -s b;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 100 -m 1 1 -s i;

./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 1000 -m 1 1 -s r;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 1000 -m 1 1 -s t;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 1000 -m 1 1 -s t a;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 1000 -m 1 1 -s b;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 1000 -m 1 1 -s i;

./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 10000 -m 1 1 -s r;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 10000 -m 1 1 -s t;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 10000 -m 1 1 -s t a;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 10000 -m 1 1 -s b;
./runTimeStepCases.py -T 2000 --dim 120 120 20 20 --perm h 10000 -m 1 1 -s i;

