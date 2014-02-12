#!/usr/bin/python
from subprocess import call
import subprocess
import sys
from os import path

solver = "r"
filename = "../bin/testCase1"
end_time = "2"

i = 1
while(i < len(sys.argv)):
    arg = sys.argv[i]
    if(arg == "-s"):
        i = i + 1
        solver = sys.argv[i]
    elif(arg == "-f"):
        i = i + 1
        filename = sys.argv[i]
    elif(arg == "-T"):
        i = i + 1
        end_time = sys.argv[i]
    else:
        print("Warning: Invalid argument: " + sys.argv[i])
    i = i + 1
    
#v_dt = [0.01,0.05,0.1,0.5,1.0,3,5,7,10]
v_dt = [1,25,50,75,100,150,200,250,300, 400, 500, 600, 700, 800, 900, 1000]

arguments = [filename,"-s",solver,"-T",end_time,"-t"]
dataoutfilename = path.basename(filename)+"-cputvsdt-s-"+solver+"-T-"+str(end_time)+".data"
dataoutfile = open(dataoutfilename,'w')
for dt in v_dt:
    print("dt = " + str(dt) + ": ")
    arguments.append(str(dt))
    with open('tempData.data','w') as output_f:
        p = subprocess.Popen(arguments,stdout=output_f,stderr=output_f)
        p.wait()
    arguments.pop()
    
    with open('tempData.data','r') as output_f:
        lines = output_f.readlines()
    line = lines[len(lines)-1]
    time = float(line.split()[3])
    dataoutfile.write(str(dt) + "\t " + str(time) + "\n") 
    print(time)

dataoutfile.close()
