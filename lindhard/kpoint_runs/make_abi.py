import netCDF4 as nc
import numpy as np
import os

nper = 14400 # kpoints per run
f_name = 'kpts.nc' # output with kpoints in it


with open('template.abi','r') as f_in:
    template = f_in.readlines()

dat = nc.Dataset(f_name)

print(dat)
exit()

nkpt = dat['nkpt'][0]
kpt = dat['kpt'][:].reshape((nkpt,3))
wtk = dat['wtk'][:]

nruns = nkpt//nper
for i in range(nruns):
    os.system(f'mkdir eigs_{i+1}')

    tmp = kpt[i*nper:(i+1)*nper,:]

    with open(f'eigs_{i+1}/eigs.abi','w') as f_out:
        for line in template:
            f_out.write(line)
        f_out.write(f'nkpt         {nper}\n')
        f_out.write(f'kpt          {tmp[0,0]:2.12f} {tmp[0,1]:2.12f} {tmp[0,2]:2.12f}\n')
        for j in range(1,nper):
            f_out.write(f'             {tmp[j,0]:2.12f} {tmp[j,1]:2.12f} {tmp[j,2]:2.12f}\n')

with open('run_eigs.py','w') as f_run:
    f_run.write('import os\n')
    f_run.write('from timeit import default_timer as timer\n\n')
    f_run.write(f'for j in range({nruns}):\n')
    f_run.write('\tstart_time = timer()\n')
    f_run.write('\tos.system(f\'cd eigs_{j+1} && mpirun -np 16 '
        '/home/ty/program_files/abinit-9.2.2/src/98_main/abinit *.abi > log 2> err && cd ../\')\n')
    f_run.write('\tend_time = timer()\n')
    f_run.write('\telapsed_time = (end_time-start_time)/60\n')
    f_run.write('\tprint(f\'\\n\\tElapsed time:\\t{elapsed_time:2.3f} minutes\\n\')\n')

