#!/bin/bash
nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 1 10 -s r;
nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 1 10 -s t;
nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 1 10 -s a;
nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 1 10 -s b;
nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 1 10 -s i;

nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 1 1 -s r;
nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 1 1 -s t;
nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 1 1 -s a;
nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 1 1 -s b;
nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 1 1 -s i;

nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 10 1 -s r;
nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 10 1 -s t;
nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 10 1 -s a;
nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 10 1 -s b;
nice ./runTimeStepCases.py -f ../bin/hom3d --timetrans -T 300 -n 30 --dim 10 10 10 10 60 10 -i 5 --perm h 10 -m 10 1 -s i;
