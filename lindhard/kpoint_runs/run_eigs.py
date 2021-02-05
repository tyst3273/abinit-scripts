import os
from timeit import default_timer as timer

for j in range(51,100):
	start_time = timer()
	os.system(f'cd eigs_{j+1} && mpirun -np 16 /home/ty/program_files/abinit-9.2.2/src/98_main/abinit *.abi > log 2> err && cd ../')
	end_time = timer()
	elapsed_time = (end_time-start_time)/60
	print(f'\n\tElapsed time:\t{elapsed_time:2.3f} minutes\n')
