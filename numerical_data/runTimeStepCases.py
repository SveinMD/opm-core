#!/usr/bin/python
from subprocess import call
import subprocess
import sys
from os import path

solver = "r"
filename = "../bin/testCase2"
end_time = "2"
n_runs = 10
solver_flag = False
print_flag = False
muw = "1"
muo = "1"
sizex = "10"
sizey = "10"
nx = "20"
ny = "20"
posx = "0"
posy = "0"
srcVol = "0.2"
layer = "0"
perm = "10"
inhom = True

i = 1
while(i < len(sys.argv)):
    arg = sys.argv[i]
    #print(arg)
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
    elif(arg == "-n"):
	i = i + 1
	n_runs = int(sys.argv[i])
    elif(arg == "-p"):
	print_flag = True
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
    elif(arg == "--dim"):
        i += 1
        sizex = sys.argv[i]
        i += 1
        sizey = sys.argv[i]
        i += 1
        ncellx = sys.argv[i]
        i += 1
        ncelly = sys.argv[i]
    elif(arg == "-i"):
        i += 1
	srcVol = sys.argv[i]
    else:
        print("Warning: Invalid argument: " + sys.argv[i])
    i = i + 1

print("Runing " + path.basename(filename) + " " + str(n_runs) + " time(s) with solver " + solver + " and sim. time " + end_time)    

#v_dt = [0.1]
v_dt = [10,25,50,75,100,150,200,250,300, 400, 500, 600, 700, 800, 900, 1000]

arguments = [filename,"-s",solver,"--dim",sizex,sizey,nx,ny,"--perm","i",layer,posx,posy,"-i",srcVol,"-m",muw,muo,"-T",end_time,"-t"]
extra_param = ""
if(solver_flag):
    extra_param = "a"
    arguments.insert(3,extra_param)
if(print_flag):
    arguments.insert(1,"-p")
perm_str = "perm-i-"+layer+"-"+posx+"-"+posy
if(not inhom):
    perm_str = "perm-h-"+perm
    perm_ind = arguments.index("--perm")
    del arguments[perm_ind+1:perm_ind+5]
    arguments.insert(perm_ind+1,"h")
    arguments.insert(perm_ind+2,perm)
print arguments
dataoutfilename = path.basename(filename)+"-cputvsdt-s-"+solver+extra_param+"-T-"+str(end_time)+"-m-"+muw+"-"+muo+"-dim-"+sizex+"-"+sizey+"-"+nx+"-"+ny+"-"+perm_str+"-i-"+srcVol
dataoutfilename = dataoutfilename.replace(".","_")
dataoutfilename = dataoutfilename+".data"
dataoutfile = open(dataoutfilename,'w')
dataoutfile.write("dt,\t cputime\n")
for dt in v_dt:
    print("dt = " + str(dt) + ": ")
    arguments.append(str(dt))
    tot_time = 0
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
    dataoutfile.write(str(dt) + ",\t " + str(tot_time/n_runs) + "\n") 
    print(tot_time/n_runs)
dataoutfile.close()
