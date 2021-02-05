import numpy as np
import Constants


class params:

    def __init__(self):

        self.sqw_prefactor = True
        self.debye_waller_factor = True
        self.debye_waller_grid = [9,9,9]
        self.experiment_type = 'INS' # INS or IXS, determines scattering lengths/form factors
        self.degeneracy_tolerance_decimals = 2 # decimal places to round
        self.sum_degenerate_bands = True
        self.poscar = 'POSCAR'
        self.temperature = 100 # Kelvin
        self.bose_factor = False
        self.Ei = 100 # meV incident energy
        self.save_figs = False
        self.clean_directory = True
        self.random_background = False

        ### PHONOPY options
        self.supercell = [2,2,2]        
        self.primitive_matrix = 'auto'

        ### ABINIT options
        self.ngqpt = [4,4,4]
        self.ddb_filename = 'si.ddb.out'

    def parse_unitcell(self):

        constants = Constants.constants()       
        self.lattice_vectors = np.zeros((3,3))  
        self.cartesian = False                  

        with open(self.poscar,'r') as fin: 
            fin.readline() 

            scale = float(fin.readline().strip().split()[0])
            self.lattice_vectors[0,:] = fin.readline().strip().split()  
            self.lattice_vectors[1,:] = fin.readline().strip().split()
            self.lattice_vectors[2,:] = fin.readline().strip().split()
            self.lattice_vectors = self.lattice_vectors*scale
            self.reciprocal_lattice()

            poscar_types = fin.readline().strip().split()   
            num_types = fin.readline().strip().split()         
            self.atom_types = []                            
            self.masses = {}

            for i in range(len(poscar_types)):
                for j in range(int(num_types[i])):
                    self.atom_types.append(poscar_types[i])     
                    self.masses[f'{j}'] = constants.atomic_masses[poscar_types[i]]

            self.position_mode = fin.readline().strip().split()[0]
            if self.position_mode == 'C':                           
                self.cartesisian = True

            self.num_atoms = len(self.atom_types)      
            self.mass_array = np.zeros(self.num_atoms)
            positions = np.zeros((self.num_atoms,3))
            self.crystal_positions = np.zeros((self.num_atoms,3))
            self.cart_positions = np.zeros((self.num_atoms,3))

            for i in range(self.num_atoms):            
                positions[i,:] = fin.readline().strip().split()    
            
            print('\n#################### INFO #######################\n')
            if self.cartesian == True:
                print('\n\tAtomic positions in cartesian coordinates\n')
                self.cart_positions = np.copy(positions)
            else:
                print('\n\tAtomic positions in crystal coordinates\n\tConverting to cartesian\n')
                self.crystal_positions = np.copy(positions)
                for i in range(self.num_atoms):
                    self.cart_positions[i,:] = np.matmul(self.lattice_vectors,self.crystal_positions[i,:])

    def reciprocal_lattice(self):
    
        self.reciprocal_lattice_vectors = np.zeros((3,3))  
        self.cell_volume = self.lattice_vectors[0,:].dot(np.cross(self.lattice_vectors[1,:],
            self.lattice_vectors[2,:]))
        self.reciprocal_lattice_vectors[0,:] = 2*np.pi*np.cross(self.lattice_vectors[1,:],
                self.lattice_vectors[2,:])/self.cell_volume
        self.reciprocal_lattice_vectors[1,:] = 2*np.pi*np.cross(self.lattice_vectors[2,:],
                self.lattice_vectors[0,:])/self.cell_volume
        self.reciprocal_lattice_vectors[2,:] = 2*np.pi*np.cross(self.lattice_vectors[0,:],
                self.lattice_vectors[1,:])/self.cell_volume





















