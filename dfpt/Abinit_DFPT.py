import numpy as np
import os

class workflow:

    def __init__(self,input_file,sbatch_job_file):
            
        self.input_file = input_file
        self.sbatch_job_file = sbatch_job_file
        self.extra = ''
        self.parse_input_file()

    def parse_input_file(self):

        with open(self.input_file) as fin:
            num_lines = len(fin.readlines())

        with open(self.input_file,'r') as fin:
            for l in range(num_lines):
                tmp = fin.readline().strip().split('=')

                if tmp[0].strip() == 'EFIELD':
                    if tmp[1].strip()[0] == 'T' or tmp[1].strip()[0] == 't':
                        self.efield = True
                    else:
                        self.efield = False
                if tmp[0].strip() == 'NUM_PROC':
                    self.num_proc = int(tmp[1].strip())

                if tmp[0].strip() == 'ABI_VERSION':
                    self.abi_version = tmp[1].strip()

                if tmp[0].strip() == 'NGQPT':
                    self.ngqpt = tmp[1].strip()
                if tmp[0].strip() == 'NUM_QPOINTS':
                    self.num_qpoints = int(tmp[1].strip())
                if tmp[0].strip() == 'QPOINTS':
                    self.qpoints = np.zeros((self.num_qpoints,3))
                    self.qpoints[0,:] = np.array(tmp[1].strip().split()).astype(float)
                    for i in range(self.num_qpoints-1):
                        self.qpoints[i+1,:] = np.array(fin.readline().strip().split()).astype(float)
                        l=l+1
                        
                if tmp[0].strip() == 'GS_TOLVRS':
                    self.gs_tolvrs = tmp[1].strip()
                if tmp[0].strip() == 'PERT_TOLVRS':
                    self.pert_tolvrs = tmp[1].strip()
                if tmp[0].strip() == 'NSCF_TOLWFR':
                    self.nscf_tolwfr = tmp[1].strip()
                if tmp[0].strip() == 'NLINE':
                    self.nline = tmp[1].strip()
    

                if tmp[0].strip() == 'IS_GGA':
                    if tmp[1].strip() == 'T':
                        self.is_gga = True
                    else:
                        self.is_gga = False

                if tmp[0].strip() == 'PP_DIRPATH':
                    self.pp_dirpath = tmp[1].strip()
                if tmp[0].strip() == 'PSEUDOS':
                    self.pseudos = tmp[1].strip()
                if tmp[0].strip() == 'NTYPAT':
                    self.ntypat = tmp[1].strip()
                if tmp[0].strip() == 'NATOM':
                    self.natom = int(tmp[1].strip())
                if tmp[0].strip() == 'TYPAT':
                    self.typat = tmp[1].strip()
                if tmp[0].strip() == 'ZNUCL':
                    self.znucl = tmp[1].strip()
                if tmp[0].strip() == 'ACELL':
                    self.acell = tmp[1].strip()
                if tmp[0].strip() == 'RPRIM':
                    self.rprim = np.zeros((3,3))
                    self.rprim[0,:] = np.array(tmp[1].strip().split()).astype(float)
                    for i in range(2):
                        self.rprim[i+1,:] = np.array(fin.readline().strip().split()).astype(float)
                        l=l+1
                if tmp[0].strip() == 'XRED':
                    self.xred = np.zeros((self.natom,3))
                    self.xred[0,:] = np.array(tmp[1].strip().split()).astype(float)
                    for i in range(self.natom-1):
                        self.xred[i+1,:] = np.array(fin.readline().strip().split()).astype(float)
                        l=l+1

                if tmp[0].strip() == 'NBAND':
                    self.nband = tmp[1].strip()
                if tmp[0].strip() == 'NSPPOL':
                    self.nsppol = tmp[1].strip()
                if tmp[0].strip() == 'NSPDEN':
                    self.nspden = tmp[1].strip()
                if tmp[0].strip() == 'SPINAT':
                    self.spinat = tmp[1].strip()

                if tmp[0].strip() == 'NGFFT':
                    self.ngfft = tmp[1].strip()
                if tmp[0].strip() == 'NGFFTDG':
                    self.ngfftdg = tmp[1].strip()

                if tmp[0].strip() == 'NGKPT':
                    self.ngkpt = tmp[1].strip()
                if tmp[0].strip() == 'NSHIFTK':
                    self.nshiftk = tmp[1].strip()
                if tmp[0].strip() == 'SHIFTK':
                    self.shiftk = tmp[1].strip()

                if tmp[0].strip() == 'ECUT':
                    self.ecut = tmp[1].strip()
                if tmp[0].strip() == 'PAWECUTDG':
                    self.pawecutdg = tmp[1].strip()
                if tmp[0].strip() == 'NSTEP':
                    self.nstep = tmp[1].strip()
                if tmp[0].strip() == 'OCCOPT':
                    self.occopt = tmp[1].strip()
                if tmp[0].strip() == 'TSMEAR':
                    self.tsmear = tmp[1].strip()

                # extra options
                if tmp[0].strip() == 'EXTRA':
                    while True:
                        tmp2 = fin.readline()
                        if len(tmp2.strip().split()) != 0:
                            self.extra = self.extra+tmp2
                            l = l+1
                        else:
                            break


#######################################################################################################################
### Phonons from DFPT 
########################################################################################################################

    def write_abinit_phonon_input_files(self):
        
        ### Header: crystal info, etc.
        file_header = ('pp_dirpath          {}\n'
                       'pseudos             {}\n'
                       'ntypat              {}\n'
                       'natom               {}\n'
                       'typat               {}\n'
                       'znucl               {}\n'
                       'acell               {}\n'
                       'rprim               {} {} {}\n'
                       '                    {} {} {}\n'
                       '                    {} {} {}\n'
                       'xred                {} {} {}\n'.format(self.pp_dirpath,
                        self.pseudos,self.ntypat,self.natom,self.typat,self.znucl,self.acell,self.rprim[0,0],
                        self.rprim[0,1],self.rprim[0,2],self.rprim[1,0],self.rprim[1,1],self.rprim[1,2],
                        self.rprim[2,0],self.rprim[2,1],self.rprim[2,2],self.xred[0,0],self.xred[0,1],self.xred[0,2]))

        for i in range(self.natom-1):
            file_header = file_header+'                    {} {} {}\n'.format(
                    self.xred[i+1,0],self.xred[i+1,1],self.xred[i+1,2])

        file_header = file_header+('\n'
            'nband               {}\n'
            'nsppol              {}\n'
            'nspden              {}\n'
            'spinat              {}\n\n'
            'ngkpt               {}\n'
            'nshiftk             {}\n'
            'shiftk              {}\n\n'
            'nline               {}\n'
            'ecut                {}\n'
            'pawecutdg           {}\n'
            'nstep               {}\n'
            'occopt              {}\n'
            'tsmear              {}\n\n'.format(self.nband,self.nsppol,self.nspden,self.spinat,self.ngkpt,self.nshiftk,
                self.shiftk,self.nline,self.ecut,self.pawecutdg,self.nstep,self.occopt,self.tsmear))

        file_header = file_header+self.extra+'\n'

        ### GROUND STATE
        gs_file_text = ('#### GROUND STATE SCF CALCULATION ####\n\n'
                        'ndtset              1\n'
                        'autoparal           1\n\n')+ file_header

        if self.is_gga == True:
            gs_file_text = gs_file_text+'pawxcdev            0\n\n'
        else:
            gs_file_text = gs_file_text+'pawxcdev            1\n\n'

        gs_file_text = gs_file_text + ('ngfft               {}\n'
                                       'ngfftdg             {}\n\n'
                                       '# GS on reduced grid of kpoints\n'
                                       'prtwf               1\n'
                                       'prtden              1\n'
                                       'kptopt              1\n'
                                       'tolvrs              {}\n\n'.format(self.ngfft,self.ngfftdg,self.gs_tolvrs))

        ### GAMMA
        if self.efield == True:
            print('\n\tNOTE: Preparing gamma point file *WITH* efield perturbation!\n')
            gamma_file_text = ('#### GAMMA POINT RESPFN WITH EFIELD PERT ####\n\n'
                            'ndtset              2\n\n')+ file_header

            if self.is_gga == True:
                gamma_file_text = gamma_file_text+'pawxcdev            0\n\n'
            else:
                gamma_file_text = gamma_file_text+'pawxcdev            1\n\n'

            gamma_file_text = gamma_file_text + ('getwfk_filepath     \"../ground_state/gso_DS1_WFK\"\n' #\"../ground_state/gs.o_DS1_WFK\"\n'
                                                 'getden_filepath     \"../ground_state/gso_DS1_DEN\"\n' #\"../ground_state/gs.o_DS1_DEN\"\n'
                                                 'prtwf               1\n'
                                                 'prtden              1\n\n'
                                                 'ngfft               {}\n'
                                                 'ngfftdg             {}\n\n'
                                                 'kptopt              2\n'
                                                 'nqpt                1\n'
                                                 'qpt                 0 0 0\n\n'
                                                 '# d/dk wave functions (?)\n'
                                                 'iscf1              -3\n'
                                                 'tolwfr1             {}\n'
                                                 'rfelfd1             2\n\n'
                                                 '# gamma point phonons\n'
                                                 'prtwf2              0\n'
                                                 'getddk2            -1\n'
                                                 'rfasr2              1\n'
                                                 'rfelfd2             3\n'
                                                 'rfphon2             1\n'
                                                 'rfdir2              1 1 1\n'
                                                 'rfatpol2            1 {}\n'
                                                 'tolvrs2             {}\n'.format(self.ngfft,self.ngfftdg,self.nscf_tolwfr,
                                                                                self.natom,self.pert_tolvrs))
        else:
            print('\n\tNOTE: Preparing gamma point file *WITHOUT* efield perturbation!\n')
            gamma_file_text = ('#### GAMMA POINT RESPFN WITHOUT EFIELD PERT ####\n\n'
                            'ndtset              1\n\n')+ file_header

            if self.is_gga == True:
                gamma_file_text = gamma_file_text+'pawxcdev            0\n\n'
            else:
                gamma_file_text = gamma_file_text+'pawxcdev            1\n\n'

            gamma_file_text = gamma_file_text + ('getwfk_filepath     \"../ground_state/gso_DS1_WFK\"\n'
                                                 'getden_filepath     \"../ground_state/gso_DS1_DEN\"\n'
                                                 'prtwf               0\n'
                                                 'prtden              0\n\n'
                                                 'ngfft               {}\n'
                                                 'ngfftdg             {}\n\n'
                                                 'kptopt              2\n'
                                                 'nqpt                1\n'
                                                 'qpt                 0 0 0\n\n'
                                                 '# gamma point phonons\n'
                                                 'rfasr               1\n'
                                                 'rfphon              1\n'
                                                 'rfdir               1 1 1\n'
                                                 'rfatpol             1 {}\n'
                                                 'tolvrs              {}\n'.format(self.ngfft,self.ngfftdg,self.natom,
                                                                                self.pert_tolvrs))


        ### Loop over Qpoints
        q_point_files = []
        for q in range(1,self.num_qpoints):
            tmp = ('#### PHONON PERT FOR QPOINT NUMBER {}, QPT = {} {} {}\n\n'
                  'ndtset              2\n\n'.format(q,
                    self.qpoints[q,0],self.qpoints[q,1],self.qpoints[q,2])) + file_header

            if self.is_gga == True:
                tmp = tmp+'pawxcdev            0\n\n'
            else:
                tmp = tmp+'pawxcdev            1\n\n'

            tmp = tmp+('nqpt                1\n'
                       'qpt                 {} {} {}\n'
                       'prtwf               1\n'
                       'prtden              1\n\n'
                       'ngfft               {}\n'
                       'ngfftdg             {}\n\n'
                       'kptopt              3\n\n'
                       '# nscf k+q\n'
                       'getden_filepath1    \"../ground_state/gso_DS1_DEN\"\n'
                       'iscf1              -2\n'
                       'tolwfr1             {}\n\n'
                       '# pert SCF\n'
                       'getwfk_filepath2    \"../ground_state/gso_DS1_WFK\"\n'
                       'prtwf2              0\n'
                       'getwfq2            -1\n'
                       'rfasr2              1\n'
                       'rfphon2             1\n'
                       'rfdir2              1 1 1\n'
                       'rfatpol2            1 {}\n'
                       'tolvrs2             {}\n'.format(self.qpoints[q,0],self.qpoints[q,1],self.qpoints[q,2],
                                                self.ngfft,self.ngfftdg,self.nscf_tolwfr,self.natom,self.pert_tolvrs))
            q_point_files.append(tmp)

        dirs = 'ground_state '
        for q in range(self.num_qpoints):
            dirs = dirs+' {}'.format(q)
        os.system('mkdir '+dirs)

        with open('ground_state/gs.abi','w') as gs_out:
            gs_out.write(gs_file_text)

        with open('0/run_0.abi','w') as run_0_out:
            run_0_out.write(gamma_file_text)

        for q in range(self.num_qpoints-1):
            with open('{}/run_{}.abi'.format(q+1,q+1),'w') as q_out:
                q_out.write(q_point_files[q])
        
        self.write_mrgddb_input_file()
        self.prep_anaddb_phonon_template()
        self.write_phonon_plot_script()

    def write_mrgddb_input_file(self):

        if self.efield == True:
            with open('mrgddb.in','w') as mrgddb:
                mrgddb.write('run.ddb.out\n')
                mrgddb.write('phonons on {} {} {} mesh\n'.format(self.ngqpt.split()[0],
                    self.ngqpt.split()[1],self.ngqpt.split()[2]))
                mrgddb.write('{}\n'.format(self.num_qpoints+1))
                mrgddb.write('ground_state/gso_DS1_DDB\n')
                for q in range(self.num_qpoints):
                    mrgddb.write('{}/run_{}o_DS2_DDB\n'.format(q,q))

        else:
            with open('mrgddb.in','w') as mrgddb:
                mrgddb.write('run.ddb.out\n')
                mrgddb.write('phonons on {} {} {} mesh\n'.format(self.ngqpt.split()[0],
                    self.ngqpt.split()[1],self.ngqpt.split()[2]))
                mrgddb.write('{}\n'.format(self.num_qpoints+1))
                mrgddb.write('ground_state/gso_DS1_DDB\n')
                mrgddb.write('0/run_0o_DS1_DDB\n')
                for q in range(1,self.num_qpoints):
                    mrgddb.write('{}/run_{}o_DS2_DDB\n'.format(q,q))

    def prep_anaddb_phonon_template(self):

        anaddb_text = (' ddb_filepath \"run.ddb.out\" \n ifcflag 1 \n asr 1 \n chneut 1 \n '
        'dipdip 1 \n eivec 2 \n enunit 1 \n brav 1 \n ngqpt {} {} {} \n nqshft 1 \n '
        'q1shft 0 0 0 \n'.format(self.ngqpt.split()[0],self.ngqpt.split()[1],self.ngqpt.split()[2]))
        with open('anaddb_template','w') as anaddb:
            anaddb.write(anaddb_text)
    
    def write_phonon_plot_script(self):

        with open('plot.py','w') as fout:
            fout.write('import numpy as np\n'
                       'import matplotlib.pyplot as plt\n'
                       'frq = np.loadtxt(\'run.abo_PHFRQ\')[:,1:]*27211.4253\n'
                       'nband = frq.shape[1]\n'
                       'fig, ax = plt.subplots()\n'
                       'fig.set_size_inches(8,6,forward=True)\n'
                       'fig.tight_layout(pad=3)\n'
                       'for i in range(nband):\n'
                       '\tplt.plot(frq[:,i],marker=\'o\',ms=0.0,mfc=\'w\',mew=0.0,mec=\'k\',lw=\'1.5\',color=\'k\')\n'
                       'for axis in [\'top\',\'bottom\',\'left\',\'right\']:\n'
                       '\tax.spines[axis].set_linewidth(1.1)\n'
                       'ax.set_ylim(0,frq[:,1:].max()*1.1)\n'
                       'ax.set_xlim(0,frq.shape[0])\n'
                       'ax.minorticks_on()\n'
                       'ax.tick_params(which=\'both\', width=1, labelsize=\'x-large\')\n'
                       'ax.tick_params(which=\'major\', length=5)\n'
                       'ax.tick_params(which=\'minor\', length=3, color=\'k\')\n'
                       'ax.set_ylabel(\'meV\',labelpad=3.0,fontweight=\'normal\',fontsize=\'x-large\')\n'
                       'plt.show()')



#######################################################################################################################

    def write_desktop_run_script(self):

        with open('run_jobs.py','w') as fout:
            fout.write('import os\n\n')
            fout.write('num_qpoints = {}\n\n'.format(self.num_qpoints))
            fout.write(f'os.system(\'cd ground_state && mpirun -np {self.num_proc:g} /home/ty/program_files/abinit-{self.abi_version}/src/98_main/abinit'
                    ' *.abi > log 2> err && cd ../\')\n'
                    'for i in range(num_qpoints):\n'
                    f'\tos.system(\'cd {{}} && mpirun -np {self.num_proc:g} /home/ty/program_files/abinit-{self.abi_version}/src/98_main/abinit'
                    ' *.abi > log 2> err && rm *WF* && cd ../\'.format(i))\n'
                    f'os.system(\'/home/ty/program_files/abinit-{self.abi_version}/src/98_main/mrgddb < mrgddb.in > mrgddb_log\')')


    def write_eagle_run_script(self):
	
        with open(self.sbatch_job_file,'r') as fid:
            job_text = fid.readlines()
        with open('ground_state/job.sh','w') as fout:
            for line in job_text:
                try:
                    if line.strip().split('=')[0] == '#SBATCH --job-name':
                        fout.write(line.strip()+'-gs\n')
                    else:
                        fout.write(line)
                        continue
                except:
                    fout.write(line)
                    continue
        for q in range(self.num_qpoints):
            with open('{}/job.sh'.format(q),'w') as fout:
                for line in job_text:
                    try:
                        if line.strip().split('=')[0] == '#SBATCH --job-name':
                            fout.write(line.strip()+'-{}\n'.format(q))
                        else:
                            fout.write(line)
                            continue
                    except:
                        fout.write(line)
                        continue
        with open('sub_q_point_jobs.py','w') as fout:
            fout.write('import os\n\n')
            for q in range(self.num_qpoints):
                fout.write('os.system(\'cd {} && sbatch job.sh && cd ../\')\n'.format(q))
	
		










                

