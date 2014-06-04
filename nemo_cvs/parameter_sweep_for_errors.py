import os
stepsize = .05
minimum = 0
maximum = 45
steps = int((maximum-minimum)/stepsize + 1)
x = 0
y = 0
z = 0
vx = 0
vy = 0
vz = 0
for i in range(0,steps):
	y = minimum+stepsize*i
	command = "mkorbit - x=" + str(x) + " y=" + str(y) + " z=" + str(z) + " vx="+ str(vx) + " vy=" + str(vy) + " vz=" + str(vz) + " potname=mpc | orbint - - nsteps=1000 dt=0.05 ndiag=1000 > /dev/null 2> param_sweep_for_errors_output.txt"
	os.system(command)
	finderrors = False
        if (os.system("grep nan param_sweep_for_errors_output.txt") == 0):
		finderrors = True
        if (os.system("grep inf param_sweep_for_errors_output.txt") == 0):
		finderrors = True
        if (os.system("grep err param_sweep_for_errors_output.txt") == 0):
		finderrors = True
	if (finderrors == True):
		addstring = "x=" + str(x) + " y=" + str(y) + " z=" + str(z) + " vx="+ str(vx) + " vy=" + str(vy) + " vz=" + str(vz) 
		if (os.system("grep \"" + addstring + "\" errors_found_in_sweep.txt") != 0):
			os.system("echo " + addstring + ">> errors_found_in_sweep.txt")
os.system("cat errors_found_in_sweep.txt")
