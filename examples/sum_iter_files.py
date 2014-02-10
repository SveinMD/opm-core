#!/usr/bin/python
import sys

if(len(sys.argv) < 3):
    print("Error: Expected file argument")
i = 1
while i < len(sys.argv):
    if(sys.argv[i] == "-f"):
        file_base = sys.argv[i+1]
        i = i + 1
    elif(sys.argv[i] == "-n"):
        n_files = int(sys.argv[i+1])
	i = i + 1
    else:
        print("Warning: Unknown flag " + sys.argv[i] + ".")
    i = i + 1

sum = 0
for i in range(0,n_files):
   file = open(file_base + str(i).zfill(3) + ".txt")
   lines = file.readlines()
   for line in lines:
	words = line.split()
	sum = sum + int(words[1])

print(sum)
