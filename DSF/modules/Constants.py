
class constants:

    def __init__(self):

        self.proton2electron_mass = 1836.152 # N x eMass
        self.thz2mev = 4.13567
        self.hartree2meV = 27211.382543519
        self.hbar = 6.582119569e-16 # eV.s
        self.kb = 8.6173e-5 # eV.K^-1

        self.atomic_masses = {'H':1.00797, # AMU
                'He':4.00260,
                'Li':6.941,
                'Be':9.01218,
                'B':10.81,
                'C':12.011,
                'N':14.0067,
                'O':15.9994,
                'F':18.998403,
                'Ne':20.179,
                'Na':22.98977,
                'Mg':24.305,
                'Al':26.98154,
                'Si':28.0855,
                'P':30.97376,
                'S':32.06,
                'Cl':35.453,
                'K':39.0983,
                'Ar':39.948,
                'Ca':40.08,
                'Sc':44.9559,
                'Ti':47.90,
                'V':50.9415,
                'Cr':51.996,
                'Mn':54.9380,
                'Fe':55.847,
                'Ni':58.70,
                'Co':58.9332,
                'Cu':63.546,
                'Zn':65.38,
                'Ga':69.72,
                'Ge':72.59,
                'As':74.9216,
                'Se':78.96,
                'Br':79.904,
                'Kr':83.80,
                'Rb':85.4678,
                'Sr':87.62,
                'Y':88.9059,
                'Zr':91.22,
                'Nb':92.9064,
                'Mo':95.94,
                'Ru':101.07,
                'Rh':102.9055,
                'Pd':106.4,
                'Ag':107.868,
                'Cd':112.41,
                'In':114.82,
                'Sn':118.69,
                'Sb':121.75,
                'I':126.9045,
                'Te':127.60,
                'Xe':131.30,
                'Cs':132.9054,
                'Ba':137.33,
                'La':138.9055,
                'Ce':140.12,
                'Pr':140.9077,
                'Nd':144.24,
                'Sm':150.4,
                'Eu':151.96,
                'Gd':157.25,
                'Tb':158.9254,
                'Dy':162.50,
                'Ho':164.9304,
                'Er':167.26,
                'Tm':168.9342,
                'Yb':173.04,
                'Lu':174.967,
                'Hf':178.49,
                'Ta':180.9479,
                'W':183.85,
                'Re':186.207,
                'Os':190.2,
                'Ir':192.22,
                'Pt':195.09,
                'Au':196.9665,
                'Hg':200.59,
                'Tl':204.37,
                'Pb':207.2,
                'Bi':208.9804,
                'Ra':226.0254,
                'Ac':227.0278,
                'Pa':231.0359,
                'Th':232.0381,
                'Np':237.0482,
                'U':238.029} # A.M.U.
        
        self.ins_scattering_lengths = {'H':-3.7390, # femtometers
                'B':5.30,  # this is dubious, its actually 5.3-0.213i, complex
                'He':3.263,
                'Li':-1.90,
                'Be':7.79,
                'C':6.6460,
                'N':9.36,
                'O':5.803,
                'F':5.654,
                'Ne':4.566,
                'Na':3.63,
                'Mg':5.375,
                'Al':3.449,
                'Si':4.1491,
                'P':5.13,
                'S':2.847,
                'Cl':9.5770,
                'Ar':1.909,
                'K':3.67,
                'Ca':4.70,
                'Sc':12.29,
                'Ti':-3.438,
                'V':-0.3824,
                'Cr':3.635,
                'Mn':-3.73,
                'Fe':9.45,
                'Co':2.49,
                'Ni':10.3,
                'Cu':7.718,
                'Zn':5.680,
                'Ga':7.288,
                'Ge':8.185,
                'As':6.58,
                'Se':7.970,
                'Br':6.795,
                'Kr':7.81,
                'Rb':7.09,
                'Sr':7.02,
                'Y':7.75,
                'Zr':7.16,
                'Nb':7.054,
                'Mo':6.715,
                'Ru':7.03,
                'Rh':5.88,
                'Pd':5.91,
                'Ag':5.922,
                'Sn':6.225,
                'Sb':5.57,
                'Te':5.80,
                'I':5.28,
                'Xe':4.92,
                'Cs':5.42,
                'Ba':5.07,
                'La':8.24,
                'Ce':4.84,
                'Pr':4.58,
                'Nd':7.69,
                'Tb':7.38,
                'Ho':8.01,
                'Er':7.79,
                'Tm':7.07,
                'Yb':12.43,
                'Lu':7.21,
                'Hf':7.7,
                'Ta':6.91,
                'W':4.86,
                'Re':9.2,
                'Os':10.7,
                'Ir':10.6,
                'Pt':9.60,
                'Au':7.63,
                'Hg':12.692,
                'Tl':8.776,
                'Pb':9.405,
                'Bi':8.532,
                'Th':10.31,
                'U':8.417} # femtometers
       

