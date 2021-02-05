import sys
sys.path.append('modules')
### data structures
import Params
import Data
import Plot
### modules
import Phonons
import QPoints_Plugins
import Write_Files



#########################################################################
#################### Initialize class instances #########################
#########################################################################

params = Params.params()
data = Data.data()
eigenvectors = Phonons.eigenvectors()
structure_factors = Phonons.structure_factors()
qplugins = QPoints_Plugins.plugins()



##########################################################################
##################### Initialize parameters ##############################
##########################################################################
# define the calculator (abinit, phonopy+VASP), calculation properties, etc

### abinit
params.ngqpt = [4,4,2]
params.ddb_filename = 'data/hbco.ddb.out'
params.path_to_anaddb = '/home/ty/program_files/abinit-9.0.4/src/98_main/anaddb'

### phonopy
params.supercell = [5,5,5]
params.force_sets_filename = 'data/FORCE_SETS.si'

### generic parameters
params.poscar = 'data/POSCAR.hbco'
params.experiment_type = 'INS' # Either IXS or INS, determines scattering lengths/atomic form factors
params.temperature = 5
params.Ei = 100
params.save_figs = False
params.degeneracy_tolerance_decimals = 1
params.sum_degenerate_bands = True
params.bose_factor = True

### read the POSCAR file
params.parse_unitcell()



############################################################################
############################## Qpoint list #################################
############################################################################
# Can devise all sorts of methods to generate lists of Q-points. 
# Qpoints.explicit_list takes a list of Qpoints and assigns it to params.Qpoints
# Qpoints.brillouin_zone_path generates a 1d cut in the BZ specified by start/end points and number of segments
# Qpoints.read_file will read a csv file of Qpoints
# This part should be plug and play 


qplugins.brillouin_zone_slice(params,Q_path=[[0,0,1],[8,0,1]],num_Qpoints=251) # generate a 1D path

#qplugins.brillouin_zone_path(params,Q_list=[[[2,0,0],[6,0,0]],[[2,1,0],[6,1,0]],[[2,0,1],[6,0,1]]],
#                                            num_Qpoints_in_shortest_segment=251) # generate a path between points

### for Si to compare to expt
#qplugins.brillouin_zone_path(params,Q_list=[[[0.0,0.0,0.0],[0.0,0.5,0.5],[0.25,0.625,0.625]],
#                                            [[0.375,0.75,0.375],[0.0,0.0,0.0],[0.5,0.5,0.5],[0.0,0.5,0.5]]],
#                                            num_Qpoints_in_shortest_segment=25) # generate a path between points



#################################################################################
###################### Compute phonon eigenvectors ##############################
#################################################################################
# run_phonopy calls phonopy to get eigenvectors. requires FORCE_SETS and POSCAR to be present. 
# Only intended to work with  VASP at this point, FORCE_SETS for other codes have different units.
# run_abinit calls anaddb to get eigenvectors. anaddb requires a ddb file, but I also require a POSCAR. 
# The POSCAR is needed to get unit cell info. Its a lot cleaner than reading the ddb for this. 
# For (eventually) phonon absorption, use params.qpoints = params.qpoints_plus


### phonopy 
#params.qpoints = params.qpoints_minus
#eigenvectors.run_phonopy(params,data)

### abinit
params.qpoints = params.qpoints_minus # phonon emission
eigenvectors.run_abinit(params,data)

### plot dispersion
#Plot.plot_dispersions(params,data)



#################################################################################
######################### Compute structure factors #############################
#################################################################################
# computes equation 3.120 in Squire's Theory of Thermal Neutron Scattering
# Will need 2 methods, 1 for INS, 1 for IXS
# INS will use a table of scattering lengths, IXS needs atomic form factors
# Plot.plot_structure_factors will make a quick plot, comment out to avoid
# Write_Files.write_strufac will write a csv file. comment out to avoid

### should work for both abinit and (eventually) phonopy
structure_factors.compute_sqw(params,data)

### plot structure factors with no resolution
#Plot.plot_structure_factors(params,data)
#Write_Files.write_strufac(params,data)



#################################################################################
######################### Apply resolution function #############################
#################################################################################
# So far, Gaussian is only implemented. Give fwhm, dE, and E_max. 
# Will need different functions for true INS resolution (e.g. ARCS) and for 
# IXS etc. This part takes A LOT longer than I want, need to vectorize in a
# way that doesn't use too much memory.

structure_factors.apply_gaussian_resolution(params,data,fwhm=5,dE=0.1,E_max=100)



























