#!/usr/bin/python
from subprocess import call
import subprocess
import sys
from os import path

solver = "r"
filename = "../bin/testCase2"
end_time = "2"
dt = "1"
solver_flag = False
print_flag = False
muw = "1"
muo = "1"

i = 1
while(i < len(sys.argv)):
    arg = sys.argv[i]
    if(arg == "-s"):
        i = i + 1
        solver = sys.argv[i]
	if(i+1 < len(sys.argv) and not sys.argv[i+1].startswith("-")):
	        solver_flag = True
		i += 1
    elif(arg == "-m"):
        i += 1
        muw = sys.argv[i]
	if(i+1 < len(sys.argv) and not sys.argv[i+1].startswith("-")):
	        muo = sys.argv[i+1]
		i += 1
    elif(arg == "-f"):
        i = i + 1
        filename = sys.argv[i]
    elif(arg == "-T"):
        i = i + 1
        end_time = sys.argv[i]
    elif(arg == "-t"):
        i = i + 1
        dt = sys.argv[i]	
    elif(arg == "-p"):
        i += 1
	print_flag = True
    else:
        print("Warning: Invalid argument: " + sys.argv[i])
    i = i + 1

v_dim = [5,10,20,30,40,45,50]

arguments = [filename,"-s",solver,"-m",muw,muo,"-T",end_time,"-t",dt,"-d"]
extra_param = ""
if(solver_flag):
    arguments.insert(3,"a")
    extra_param = "a"
if(print_flag):
    arguments.insert(1,"-p")
dataoutfilename = path.basename(filename)+"-cputvscellc-s-"+solver+extra_param+"-T-"+str(end_time)+"-t-" + str(dt) + "-m-"+muw+"-"+muo+".data"
dataoutfile = open(dataoutfilename,'w')
for dim in v_dim:
    print("dimension = " + str(dim) + ": ")
    arguments.append(str(dim))
    with open('tempData.data','w') as output_f:
        p = subprocess.Popen(arguments,stdout=output_f,stderr=output_f)
        p.wait()
    arguments.pop()
    
    with open('tempData.data','r') as output_f:
        lines = output_f.readlines()
    line = lines[len(lines)-1]
    time = float(line.split()[3])
    dataoutfile.write(str(dim*dim) + "\t " + str(time) + "\n") 
    print(time)

dataoutfile.close()
