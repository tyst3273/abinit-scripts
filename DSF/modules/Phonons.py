import numpy as np
import matplotlib.pyplot as plt
import os
import Constants

constants = Constants.constants()


class eigenvectors:

    def run_phonopy(self,params,data):

        print('\n\tWARNING: This doesnt work right with PHONOPY yet. The structure factors are bad.\n')

        import phonopy
        from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

        phonon = phonopy.load(unitcell_filename=params.poscar,
                supercell_matrix=params.supercell,
                primitive_matrix=params.primitive_matrix,
                force_sets_filename=params.force_sets_filename)
 
        phonon.run_band_structure([list(params.qpoints)],with_eigenvectors=True)
        self.phonopy_output = phonon.get_band_structure_dict()

        data.num_Qpoints = params.num_Qpoints
        data.num_bands = params.num_atoms*3
        data.qpoints = {}
        data.Qpoints = {}
        data.cart_Qpoints = {}
        data.frequencies = {}
        data.eigenvectors = {}
        data.cart_eigenvectors = {}

        for i in range(data.num_Qpoints):
            data.frequencies[f'{i}'] = self.phonopy_output['frequencies'][0][i]*constants.thz2mev
            data.eigenvectors[f'{i}'] = self.phonopy_output['eigenvectors'][0][i]
            data.cart_eigenvectors[f'{i}'] = np.copy(data.eigenvectors[f'{i}'])
            data.qpoints[f'{i}'] = self.phonopy_output['qpoints'][0][i]
            for b in range(data.num_bands):
                for j in range(params.num_atoms):
                    data.cart_eigenvectors[f'{i}'][b,j] = data.cart_eigenvectors[f'{i}'][b,j]/np.sqrt(params.masses[f'{j}'])
                    data.cart_eigenvectors[f'{i}'][b,j] = data.cart_eigenvectors[f'{i}'][b,j]*np.exp(1j*
                            np.dot(data.qpoints[f'{i}'],params.cart_positions[j,:]))
            data.Qpoints[f'{i}'] = params.Qpoints[i,:]
            data.cart_Qpoints[f'{i}'] = np.matmul(params.reciprocal_lattice_vectors,data.Qpoints[f'{i}'])

    ##########################################################################

    def run_abinit(self,params,data):
        
        abi_q = params.qpoints # reduced Qpoints only ?

        ## call anaddb (ABINIT) to compute phonons
        anaddb_text = (' ddb_filepath \"{}\" \n ifcflag 1 \n asr 1 \n chneut 1 \n '
        'dipdip 1 \n eivec 2 \n enunit 1 \n brav 1 \n ngqpt {} {} {} \n nqshft 1 \n '
        'q1shft 0 0 0 \n prtphbands 2 \n nph1l {} \n'.format(
            params.ddb_filename,params.ngqpt[0],params.ngqpt[1],
            params.ngqpt[2],params.num_Qpoints))
        anaddb_text = anaddb_text+' qph1l '+'{:3.6f} {:3.6f} {:3.6f} 1.0 \n '.format(abi_q[0,0],abi_q[0,1],abi_q[0,2])
        for i in range(params.num_Qpoints-1):
            anaddb_text = anaddb_text+'      '+'{:3.6f} {:3.6f} {:3.6f} 1.0 \n '.format(abi_q[i+1,0],abi_q[i+1,1],abi_q[i+1,2])

        with open('anaddb.in','w') as fid:
            fid.write(anaddb_text)
        
        os.system(params.path_to_anaddb+' anaddb.in > eigenvectors 2> err')
        self.read_abinit_eigenvectors(params,data)
        if params.clean_directory == True:
            os.system('rm *abo* *DOS* anaddb* err eigenvectors')

    def read_abinit_eigenvectors(self,params,data):

        data.num_Qpoints = params.num_Qpoints
        data.num_bands = params.num_atoms*3
        data.qpoints = {}
        data.Qpoints = {}
        data.cart_Qpoints = {}
        data.frequencies = {}
        data.cart_eigenvectors = {}

        ### Read Eigenvectors from log file
        with open('eigenvectors','r') as fid:

            ### Qpoint list
            if params.num_Qpoints == 1:
                data.qpoints['0'] = params.qpoints[0]
                data.Qpoints['0'] = params.Qpoints
                data.cart_Qpoints['0'] = np.matmul(params.reciprocal_lattice_vectors,params.Qpoints[0])
            else:
                while True:
                    tmp = fid.readline()
                    if tmp.strip() == 'qph1l':
                        break
                for i in range(params.num_Qpoints):
                    data.qpoints[f'{i}'] = np.array(fid.readline().strip().split()[:-1]).astype(float) 
                    data.Qpoints[f'{i}'] = params.Qpoints[i,:]
                    data.cart_Qpoints[f'{i}'] = np.matmul(params.reciprocal_lattice_vectors,data.Qpoints[f'{i}']) 

            ### eigendisplacements and frequencies
            print('\n\tAbinit eigendisplacements are in cartesian coordinates and are divided by sqrt(M) in atomic units\n')
            for q in range(params.num_Qpoints):
                while True:
                    tmp = fid.readline()
                    if tmp.strip() == 'Eigendisplacements':
                        frequencies = np.zeros(data.num_bands)
                        eigenvectors = np.zeros((data.num_bands,params.num_atoms,3)).astype(complex)
                        for b in range(data.num_bands):
                            while True:
                                tmp = fid.readline()
                                if tmp.split()[0].strip() == 'Mode':
                                    frequencies[b] = float(tmp.split()[-1].strip())*constants.hartree2meV
                                    for j in range(params.num_atoms):
                                        while True:
                                            tmp = fid.readline()
                                            if tmp.split()[0].strip() == f'{j+1}':
                                                real = tmp.strip().split()[1:]
                                                imag = fid.readline().strip().split()
                                                eigenvectors[b,j,0] = float(real[0])+1j*float(imag[0])
                                                eigenvectors[b,j,1] = float(real[1])+1j*float(imag[1])
                                                eigenvectors[b,j,2] = float(real[2])+1j*float(imag[2])
                                                break
                                    break
                        data.cart_eigenvectors[f'{q}'] = eigenvectors
                        data.frequencies[f'{q}'] = frequencies
                        break


                                    

#########################################################################################################

class structure_factors:

    def compute_sqw(self,params,data):

        data.structure_factors = np.zeros((data.num_bands,data.num_Qpoints))
        
        print('\n#################### SQW ########################\n')
        if params.experiment_type == 'INS':
            print('\n\tExperiment type: INS\n')
            data.b_dict = {}
            data.b = np.zeros((1,params.num_atoms))
            for i in range(params.num_atoms):
                data.b_dict[f'{i}'] = constants.ins_scattering_lengths[f'{params.atom_types[i]}']
                data.b[0,i] = constants.ins_scattering_lengths[f'{params.atom_types[i]}']
            data.b = np.tile(data.b,reps=(data.num_bands,1))  ### Femtometers
        elif params.experiment_type == 'IXS':
            print('\n\tIXS not supported yet\n\tSorry, but I am gonna crash now!\n')
            exit()
        else:
            print(f'\n\tUnknown experiment type: {params.experiment_type}\n\tSorry, but I am gonna crash now!\n')
            exit()

        print('\n\tDebye-Waller factor not implemented yet.\n')
        print('\n\tThe structure factors produced by sum_on_d match very well with SNAXS. See the images in ./pngs\n'
                '\tI havent implemented the k\'/k or 1/w_s terms yet, so this isnt the true intensity. I also need to\n'
                '\tinclude the resolution and atomic form factors for xrays. Also, this is only the phonon emission term\n'
                '\tI need to implement the phonon absorption term too to integrate energy from -infinity to infinity.\n')

        self.compute_structure_factors(params,data)
        data.structure_factors = data.structure_factors/np.matrix(data.structure_factors).max() # normalize to max = 1

    def compute_structure_factors(self,params,data):
       
        for q in range(data.num_Qpoints):
            data.q = q
            self.compute_sum_on_d(params,data)

    def compute_sum_on_d(self,params,data):

        exp_ikd = params.cart_positions # (num_atoms,3_dir) in cartesian coords
        k_vector = np.tile(data.cart_Qpoints[f'{data.q}'].reshape(1,3),reps=(params.num_atoms,1)) # (num_atoms,3_dir)
        exp_ikd = np.exp(1j*np.multiply(k_vector,exp_ikd).sum(axis=1)) # exp(i*k.d)
        
        k_dot_eig = np.tile(k_vector,reps=(data.num_bands,1)) # (num_atoms*num_bands,3_dir)
        eigendisplacements = data.cart_eigenvectors[f'{data.q}'].reshape(
                data.num_bands*params.num_atoms,3)*np.sqrt(constants.proton2electron_mass) # (num_atoms*num_bands,3_dir)

        k_dot_eig = np.multiply(k_dot_eig,eigendisplacements).sum(axis=1) # k.eig_qv,d
        k_dot_eig = k_dot_eig.reshape(data.num_bands,params.num_atoms) # (num_bands,num_atoms)

        exp_ikd = np.tile(exp_ikd.reshape(1,params.num_atoms),reps=(data.num_bands,1)) # (num_bands,num_atoms)
        exp_ikd = np.multiply(exp_ikd,data.b) # multiply by scattering length. 1/sqrt(mass) is included in eigenvectors
        
        structure_factors = np.multiply(exp_ikd,k_dot_eig).sum(axis=1) # sum over atoms 
        data.structure_factors[:,data.q] = abs(structure_factors)**2 # norm squared

    ###################################################################################3

    def apply_gaussian_resolution(self,params,data,fwhm=1,dE=0.01,E_max=100):
        """
        widths and spectra vals are in meV
        the Gaussian convolution uses a PBC for negative energy, dont use a huge FWHM or it might be bad
        Also, I assume E_max is set sufficiently above the highest phonon energy so that fringe effects 
        on the upper boundary don't matter.
        Any intensities below E=0 are ignored
        """
        print('\n################### CONVOLUTION #####################\n')
        print(f'\n\tConvolution with Gaussian function, FWHM = {fwhm} meV\n')

        data.fwhm = fwhm
        c = fwhm/2.35482

        data.dE = dE
        data.E_max = E_max
        data.spectra_E = np.arange(0,data.E_max+data.dE,data.dE)
        data.spectra_num_E = len(data.spectra_E)
        data.spectra = np.zeros((data.spectra_num_E,params.num_Qpoints))
        data.smooth_spectra = np.zeros((data.spectra_num_E,params.num_Qpoints))
        structure_factors = []
        energies = []

        ### sum intensity of degenerate bands
        if params.sum_degenerate_bands == True:
            print('\n\tSumming degenerate bands before convolution (using convolution dE as tolerance)\n')
            for q in range(params.num_Qpoints):
                sfac = data.structure_factors[:,q]
                energy = data.frequencies[f'{q}']
                reduced_energies = []
                summed_sfac = []
                while True:
                    if len(energy) == 0:
                        break
                    test_energy = energy[0]
                    reduced_energies.append(test_energy)
                    indicies = np.intersect1d(np.argwhere(energy <= (test_energy+data.dE)),
                        np.argwhere(energy > (test_energy-data.dE)))
                    summed_sfac.append(sfac[indicies].sum())
                    sfac = np.delete(sfac,indicies)
                    energy = np.delete(energy,indicies)
                energies.append(reduced_energies)
                structure_factors.append(summed_sfac)
        else:
            print('\n\tWARNING: You should definitely sum degenerate bands!!!\n')
            for q in range(params.num_Qpoints):
                energies.append(data.frequencies[f'{q}'])
                structure_factors.append(data.structure_factors[:,q])

        ### populate array for heatmap
        ### try statement takes care of negative energies
        for q in range(params.num_Qpoints):
            for b in range(len(structure_factors[q][:])):
                try: # if there are negative modes, argwhere returns an empty vector and the slice crashes
                    data.spectra[np.argwhere(data.spectra_E <= 
                        energies[q][b]).max(),q] = structure_factors[q][b]
                except:
                    continue

        if params.bose_factor == True:
            print('\n\tWARNING: Bose factor isnt verified. Need to compare to SNAXS.\n')
            if params.temperature < 5:
                temperature = 5
            else:
                temperature = params.temperature
            inds = np.argwhere(data.spectra_E <= 0.5)
            tmp_e = np.copy(data.spectra_E)
            tmp_e[inds] = 0.5
            bose = 1+1/(np.exp(tmp_e/(constants.kb*1000*temperature))-1)
            bose = np.tile(bose.reshape((data.spectra_num_E,1)),reps=(1,params.num_Qpoints))
            data.spectra = np.multiply(data.spectra,bose)
            data.spectra = data.spectra/np.max(data.spectra)

        ### gaussian convolution using for loops, slow but very little memory utilization
        g_energy = np.append(data.spectra_E-data.spectra_E.max(),data.spectra_E[1:])
        gaussian = np.exp(-0.5*g_energy**2/c**2)/c/np.sqrt(2*np.pi)
        gaussian = np.tile(gaussian.reshape((gaussian.shape[0],1)),(1,data.num_Qpoints))
        tmp = np.append(data.spectra,data.spectra,axis=0)[1:,:]
        for e in range(data.spectra_num_E):
            if e%50 == 0:
                print(f'\t------ {e}/{data.spectra_num_E} -------')
            data.smooth_spectra[e,:] = np.trapz(tmp*np.roll(gaussian,shift=e,axis=0),g_energy,axis=0)
        print('\n\tDone convolving!\n')
        data.smooth_spectra = data.smooth_spectra/np.max(data.smooth_spectra)

#        if params.random_background == True:
#            data.smooth_spectra = data.smooth_spectra+(np.random.normal(0,1,
#                (data.smooth_spectra.shape[0],data.smooth_spectra.shape[1])))*0.001
        
        plt.imshow(data.smooth_spectra,origin='lower',aspect='auto',cmap='hot')
        plt.show()





        






















