using LatticeUtilities
using Random
using LinearAlgebra

# define size of the lattice
Lx = 2
Ly = 2

# define initial electron density
n̄ = 1.0

# number of particles (in the canonical ensemble)
# nup = 8
# ndn = 8

# square unit cell
unit_cell = UnitCell([[1.0,0.0], [0.0,1.0]],           # lattice vectors
                               [[0.0,0.0]])       # basis vectors 

# build square lattice
lattice = Lattice([Lx,Ly],[true,true])

# define nearest neighbor bonds
bond_x = Bond((1,1), [1,0])
bond_y = Bond((1,1), [0,1])
# define next nearest neighbor bonds
bond_xy = Bond((1,1), [1,1])
bond_yx = Bond((1,1), [1,-1])

# vector of 2D bonds
bonds = [[bond_x, bond_y], [bond_xy, bond_yx]]


# 1D chain unit cell
# unit_cell = UnitCell([[1.0]],           # lattice vectors
#                     [[0.0]])       # basis vectors 

# # build 1D chain
# lattice = Lattice([Lx],[true])

# # define nearest neighbor bonds for a 1D chain
# bond_x = Bond((1,1),[1])

# # vector of 1D bonds
# bonds = [[bond_x]]


# nearest neighbor hopping
t = 1.0

# next nearest neighbor hopping
tp = 0.0

# chemical potential
μ = 0.0
opt_μ = true

# s-wave order parameter
Δs = 0.3    # initial value
opt_s  = true

# d-wave order parameter
Δd = 0.25    # initial value                            
opt_d = false   

# anti-ferromagnetic (Neél) order parameter
Δa = 0.5  # initial value
opt_a = false  

# uniform charge order parameter
Δc = 0.5  # intial value
opt_c = false       

# charge modulation parameter
Δcs = 0.5     # intial value
opt_cs = false     

# spin modulation parameter
Δss = 0.5     # initial value
opt_ss = false     

# whether model is particle-hole transformed
pht = true

# whether to output during runtime
verbose = true

# initialize random seed
seed = abs(rand(Int))

# initialize random number generator
rng = Xoshiro(seed)






























########################################################################################################

# Define unit cell
# square lattice 
# unit_cell = UnitCell([[1.0,0.0], [0.0,1.0]],           # lattice vectors
#                                [[0.0,0.0]])       # basis vectors 

# # square lattice with basis
# # unitcell = UnitCell([[1.0,0.0], [0.0,1.0]],           # lattice vectors
# #                     [[0.0,0.0],[0.5,0.0],[0.0,0.5]])       # basis vectors, (Cu, O, O) 

# # Define bonds
# # for single orbital model
# # nearest neighbors
# bond_1 = Bond([1,1], [1,0])
# bond_2 = Bond([1,1], [0,1])
# # # next-nearest neighbors
# bond_3 = Bond([1,1], [1,1])
# bond_4 = Bond([1,1], [1,-1])

# #single_square_bonds = [bond_1,bond_2,bond_3,bond_4]

# # for three orbital model
# # nearest neighbors (Cu-O)
# bond_1 = Bond(orbitals = [1,2], displacement = [0,0]) 
# bond_2 = Bond(orbitals = [1,3], displacement = [0,0]) 
# bond_3 = Bond(orbitals = [2,1], displacement = [1,0])
# bond_4 = Bond(orbitals = [1,3], displacement = [0,1])

# # next-nearest neighbors (O-O)
# # TODO
# bond_5 = Bond(orbitals = [2,3], displacement = [0,0])
# bond_6 = Bond(orbitals = [2,3], displacement = [1,0])
# bond_7 = Bond(orbitals = [2,3], displacement = [0,1])
# bond_8 = Bond(orbitals = [2,3], displacement = [1,1])


# lieb_bonds = [bond_1,bond_2,bond_3,bond_4,bond_5,bond_6,bond_7,bond_8]


# # Define linear dimensions,
# Lx, Ly = 2,2


# "PARTICLE CONFIGURATION" 
# # Define ensemble to be used
# # Canonical Ensemble (ce): fixed particle number 
# # Grand Canonical Ensemble (gce): allow particle number fluctuation
# gce = false
# ce = true

# # number of spin-up and spin-down particles
# Nup = 2
# Ndwn = 2

# # target electron density
# density = 1.0

# # particle-hole transformation
# PHT = false


# "VARIATIONAL PARAMETERS"
# # the 'opt_parname' flag will determine whether the appropriate term 
# # is added to the initial Hamiltonian
# # hopping and chemical potential terms are automatically added 

# # TIGHT BINDING MODEL
# #"nearest neighbor hopping amplitude"
# t = 1.0             # single-band 
# tpd = 1.0           # p-d orbital hopping 

# #"next-nearest neighbor hopping amplitude"
# tp = 0.0            # single-band
# tpp = 0.0           # p-orbital hopping 

# #"chemical potential"
# μ = 3.0             # initial value, single-band 
# μd = 0.0               # initial value, d-orbital
# μp = 0.0               # initial value, p-orbital
# opt_μ = false       # if optimized        
# # opt_μp = false
# # opt_μd = false           

# # HUBBARD MODEL
# #"on-site Coulomb repulsion"
# U = 8.0             # single-band
# Up = 0.0            # p-orbital 
# Ud = 0.0            # d-orbital

# #"charge transfer energy"
# Δpd = 2.65         # initial value, p-d orbitals (assuming tpd = 1.13 eV)
# opt_pd = true

# #"s-wave"
# Δ0 = 0.3           # initial value
# opt_s = false
      
# #"d-wave"
# Δd = 0.3           # initial value
# opt_d = false   

# #"antiferromagnetism"
# Δafm = 0.5         # initial value
# opt_afm = false     

# #"charge density wave"
# Δcdw = 0.5         # intial value
# opt_cdw = false       

# #"charge stripes"
# Δcstr = 0.5           # intial value
# opt_cstr = false        

# #"spin stripes"
# Δsstr = 0.5             # initial value
# opt_sstr = false        

# # setting both stripe orders to "true" 
# # initializes vectors of 2*Lx variational parameters 
# # for each stripe order and allows for optimization

# # TODO: pair density wave order parameter(s)?

# #"HOLSTEIN MODEL"
# #"phonon frequency"
# Ω = 0.0                 # mean value 

# #"Holstein coupling constant"
# g = 1.5                 # microscopic coupling of the form gn(b† + b)
# α = g * sqrt( 2 * Ω )   # microscopic coupling of the form αnX

# #"BOND SSH MODEL"
# #"phonon frequency"
# # Ω = 0.0                 # mean value 

# # "bSSH coupling constant"
# # g = 1.5                 # microscopic coupling of the form gn(b† + b)
# # α = g * sqrt( 2 * Ω )   # microscopic coupling of the form αnX

# #"OPTICAL SSH MODEL"

# "RANDOM SEED"

# # initialize random seed
# seed = abs(rand(Int))

# # initialize random number generator
# rng = Xoshiro(seed)

# "JASTROW FACTORS"
# jpar_init = 0.5         # intial 
# den_Jastrow = true      # set to 'true' to use desnity Jastrow
# spn_Jastrow = false     # set to 'true' to use spin Jastrow