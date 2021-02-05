
in_file = 'FORCE_SETS'
out_file = 'out_FORCE_SETS'

fin = open(in_file,'r')
fout = open(out_file,'w')

num_atoms = int(fin.readline().strip())
num_disp = int(fin.readline().strip())

fout.write(str(num_atoms)+'\n')
fout.write(str(num_disp)+'\n')

for disp in range(num_disp):
    fout.write(fin.readline())
    fout.write(fin.readline())
    tmp = fin.readline().strip().split()
    for i in range(3):
        fout.write(f'{float(tmp[i])*0.529177:0.16f} ')
    fout.write('\n')
    for atom in range(num_atoms):
        fout.write(fin.readline())

    

fin.close()
fout.close()
