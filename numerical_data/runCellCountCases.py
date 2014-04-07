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
n_runs = 10
muw = "1"
muo = "1"
sizex = "10"
sizey = "10"
ncellx = "20"
ncelly = "20"
xpos = "0"
ypos = "0"
srcVol = "0.2"
layer = "0"
inhom = True
perm = "10"

i = 1
while(i < len(sys.argv)):
    arg = sys.argv[i]
    if(arg == "-s"):
        i += 1
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
        i += 1
        filename = sys.argv[i]
    elif(arg == "-T"):
        i += 1
        end_time = sys.argv[i]
    elif(arg == "-t"):
        i += 1
        dt = sys.argv[i]	
    elif(arg == "-p"):
        i += 1
	print_flag = True
    elif(arg == "-n"):
	i += 1
	n_runs = int(sys.argv[i])
    elif(arg == "--dim"):
        i += 1
        sizex = sys.argv[i]
        i += 1
        sizey = sys.argv[i]
        i += 1
        ncellx = sys.argv[i]
        i += 1
	ncelly = sys.argv[i]
    elif(arg == "--perm"):
        i += 1
        inhom = sys.argv[i] == "i"
        if(not inhom):
            i += 1
            perm = sys.argv[i]
        else:
            i += 1
            layer = sys.argv[i]
            i += 1
            posx = sys.argv[i]
            i += 1
            posy = sys.argv[i]
    elif(arg == "-i"):
        i += 1
	srcVol = sys.argv[i]
    else:
        print("Warning: Invalid argument: " + sys.argv[i])
    i += 1

v_dim = [5,10,20,30,40,45,50]

arguments = [filename,"-s",solver,"--dim",sizex,sizey,ncellx,ncelly,"--perm","i",layer,xpos,ypos,"-i",srcVol,"-m",muw,muo,"-T",end_time,"-t",dt,"-d"]
extra_param = ""
if(solver_flag):
    arguments.insert(3,"a")
    extra_param = "a"
if(print_flag):
    arguments.insert(1,"-p")
perm_str = "perm-i-"+layer+"-"+posx+"-"+posy
if(not inhom):
    perm_str = "perm-h-"+perm
    perm_ind = arguments.index("--perm")
    del arguments[perm_ind+1:perm_ind+5]
    arguments.insert(perm_ind+1,"h")
    arguments.insert(perm_ind+2,perm)

dataoutfilename = path.basename(filename)+"-cputvscellc-s-"+solver+extra_param+"-T-"+str(end_time)+"-t-" + str(dt) + "-m-"+muw+"-"+muo+"-dim-"+sizex+"-"+sizey+"-"+perm_str+"-i-"+srcVol+".data"
dataoutfile = open(dataoutfilename,'w')
dataoutfile.write("dt,\t cputime\n")
tot_time = 0
for dim in v_dim:
    print("dimension = " + str(dim) + ": ")
    arguments.append(str(dim))
    for i in range(n_runs):
        with open('tempData.data','w') as output_f:
            p = subprocess.Popen(arguments,stdout=output_f,stderr=output_f)
            p.wait()
    
        with open('tempData.data','r') as output_f:
            lines = output_f.readlines()
        line = lines[len(lines)-1]
        time = float(line.split()[3])
        tot_time += time
    arguments.pop()
    dataoutfile.write(str(dim*dim) + ",\t " + str(tot_time/n_runs) + "\n") 
    print(time)

dataoutfile.close()
