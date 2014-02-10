#!/usr/bin/python
import sys

filename = sys.argv[1]
file = open(filename,'r')

lines = file.readlines()
iterations = 0
for line in lines:
    words = line.split()
    iterations = iterations + int(words[1])

print(iterations)
