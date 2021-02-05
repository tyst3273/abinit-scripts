
nqpoints = 10
qpt_grid = [5,5,5]

gs_scf    = '1.0d-18' # vrs
pert_scf  = '1.0d-12' # vrs
nscf      = '1.0d-22' # wfr

####################

crystal_info = """

 ntypat 1
 natom  2
 typat  1 1
 znucl  14
 acell  3*5.40139978 Angstrom
 xred   0.00 0.00 0.00
        0.25 0.25 0.25
 rprim  0.0  0.5  0.5
        0.5  0.0  0.5
        0.5  0.5  0.0

"""

####################

electron_params = """

 nband   10
 nsppol  1
 nspden  1
 nspinor 1

"""

#####################

scf_params = """

 ngkpt   23 23 23 
 nshiftk 1
 shiftk  0.0 0.0 0.0

 pawecutdg 40
 ecut      20

 nstep  75

 occopt 3
 tsmear 2 K
 
"""

#
#
#
#
#
#
#
#
#

common = """

 getwfk  2
 getden  2
 nqpt    1
 kptopt  3
 rfphon  1
 rfdir   1 1 1 

"""

gs_txt = """

 getwfk1 0
 getden1 0
 prtden1 1
 prtwf1  1
 kptopt1 1
 rfphon1 0
 nqpt1   0

"""

gs_nscf = """

 getwfk2 1
 getden2 1
 prtden2 1
 prtwf2  1
 iscf2  -2
 nqpt2   0
 rfphon2 0

"""

gamma = """

 iqpt3   1
 kptopt3 2
 getden3 0

"""

qpoint_txt = """

 qptopt  1
 nshiftq 1
 shiftq  0.0 0.0 0.0

"""

natoms = crystal_info.strip().split()
natoms = natoms[natoms.index('natom')+1]


with open('abi-dfpt.in','w') as fid:
    ### common
    fid.write('#common\n\n autoparal 1 \n ndtset {}'.format(nqpoints*2+1)) 
    fid.write(common+' rfatpol 1 {}\n\n'.format(natoms))

    ### crystal
    fid.write('#crystal info'+crystal_info)
 
    ### electrons
    fid.write('#electrons '+electron_params)

    ### scf
    fid.write('#scf params'+scf_params+' tolvrs '+pert_scf+'\n\n')

    ### GS calc
    fid.write('#gs calc - scf'+gs_txt+' tolvrs1 '+gs_scf+'\n\n')

    ### GS NSCF
    fid.write('#gs calc nscf'+gs_nscf+' tolwfr2 '+nscf+'\n\n')

    ### qpoints
    fid.write('#qpoints'+qpoint_txt+' ngqpt {} {} {}'.format(qpt_grid[0],
        qpt_grid[1],qpt_grid[2]))

    ### gamma
    fid.write('\n\n#gamma point'+gamma+'\n')

    ### the rest of them ...
    for i in range(2,nqpoints+1):
        fid.write(""" 

#qpt no {} nscf - scf
 iscf{}  -2
 iqpt{}   {}
 tolwfr{} {}
 rfphon{} 0

 iqpt{}   {}
 getden{} 0
 getwfq{} {}
 """.format(i,i*2,i*2,i,i*2,nscf,i*2,i*2+1,i,i*2+1,i*2+1,i*2))
        

