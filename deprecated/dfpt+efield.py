
efield = True

nqpoints = 6
qpt_grid = [3,3,3]

gs_scf    = '1.0d-16' # vrs
pert_scf  = '1.0d-10' # vrs
nscf      = '1.0d-18' # wfr

####################

crystal_info = """

 ntypat 3
 natom  7
 typat  1 1 2 3 3 3 3
 znucl  57 29 8
 acell  1.3166995860E+01  1.3166995860E+01  1.3166995860E+01 Bohr
 rprim -2.8851771794E-01  2.8851771794E-01  9.1297045564E-01
        2.8851771794E-01 -2.8851771794E-01  9.1297045564E-01
        2.8851771794E-01  2.8851771794E-01 -9.1297045564E-01
 xred   0.649016293234      0.649016293234      0.000000000000
        0.350983706766      0.350983706766      0.000000000000
        0.000000000000      0.000000000000      0.000000000000
        0.000000000000      0.500000000000      0.500000000000
        0.500000000000      0.000000000000      0.500000000000
        0.750000000000      0.250000000000      0.500000000000
        0.250000000000      0.750000000000      0.500000000000


"""

####################

electron_params = """

"""

#####################

scf_params = """

 ngkpt   9 9 9
 nshiftk 1
 shiftk  0.0 0.0 0.0

 pawecutdg 45
 ecut      20

 nstep  75

 occopt 7
 tsmear 0.005
 
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

if efield == True:
    gamma = """
 
 iscf3  -3
 iqpt3   1
 rfphon3 0
 rfelfd3 2
 kptopt3 2
 tolwfr3 {}

 getddk4 3
 rfelfd4 3
 iqpt4   1
 kptopt4 2
 getden4 0

    """.format(nscf)
else:
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
    if efield == True:
        shift = 1
    else:
        shift = 0
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
 """.format(i,
     i*2+shift,
     i*2+shift,
     i,
     i*2+shift,
     nscf,
     i*2+shift,
     i*2+1+shift,
     i,
     i*2+1+shift,
     i*2+1+shift,
     i*2+shift))
        

