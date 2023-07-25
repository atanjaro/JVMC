using LatticeUtilities
using Random
using LinearAlgebra
using Test
using DelimitedFiles
using BenchmarkTools

# files to include
# include("Parameters.jl")
include("Hamiltonian.jl")
include("ParticleConfiguration.jl")
include("Jastrow.jl")

#############################
## DEFINE MODEL PARAMETERS ##
#############################

# define the size of the lattice
Lx = 2
Ly = 2

# define initial electron density
n̄ = 1.0

# number of particles
# nup = 8
# ndn = 8

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
Δcm = 0.5     # intial value
opt_cm = false     

# spin modulation parameter
Δsm = 0.5     # initial value
opt_sm = false     

##################################
## DEFINE SIMULATION PARAMETERS ##
##################################

# whether model is particle-hole transformed
pht = true

# initialize random seed
seed = abs(rand(Int))

# initialize random number generator
rng = Xoshiro(seed)

# number of thermalization updates
N_burnin = 1000

# number of simulation updates
N_updates = 1000

# number of bins
N_bins = 100

# bin size
bin_size = div(N_updates, N_bins)

# whether to output during runtime
verbose = true

# whether to output matrices to files
write = false

# initialize addition simulation information dictionary
additional_info = Dict(
    "time" => 0.0,
    "N_burnin" => N_burnin,
    "N_updates" => N_updates,
    "N_bins" => N_bins,
    "bin_size" => bin_size,
    "Δs" => Δs,
    "local_acceptance_rate" => 0.0,
    "seed" => seed,
    "n_bar" => n̄
)

##################
## DEFINE MODEL ##
##################

# square unit cell
unit_cell = UnitCell([[1.0,0.0], [0.0,1.0]],           # lattice vectors
                               [[0.0,0.0]])            # basis vectors 

# build square lattice
lattice = Lattice([Lx,Ly],[true,true])

# define model geometry
model_geometry = ModelGeometry(unit_cell,lattice)

# define nearest neighbor bonds
bond_x = Bond((1,1), [1,0])
bond_y = Bond((1,1), [0,1])
# define next nearest neighbor bonds
bond_xy = Bond((1,1), [1,1])
bond_yx = Bond((1,1), [1,-1])

# vector of 2D bonds
bonds = [[bond_x, bond_y], [bond_xy, bond_yx]]

# define non-interacting tight binding model
tight_binding_model = TightBindingModel([t,tp],μ)

# define variational parameters
variational_parameters = VariationalParameters(["μ","Δs"], [μ, Δs], [opt_μ, opt_s])
    
###################################
## INITIALIZE TRIAL WAVEFUNCTION ##
###################################

# define particle density
Np, Ne, nup, ndn = get_particle_numbers(n̄)
# density, Np, Ne = get_particle_density(nup, ndn)

# construct mean-field Hamiltonian
H_mf = build_mean_field_hamiltonian()
# writedlm("H_mf.csv", H_mf)  # TODO: move to this to function, make write=true with flag

# initialize Slater determinant state and initial particle configuration
(D, pconfig, ε₀, M, U) = build_slater_determinant()  

# construct Jastrow factors
#(Tvec, jpar_matrix, num_jpars) = get_jastrow_factor()






