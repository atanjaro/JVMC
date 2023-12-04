using LatticeUtilities
using Random
using LinearAlgebra
using Test
using DelimitedFiles
using BenchmarkTools

# files to include
include("Hamiltonian.jl")
include("ParticleConfiguration.jl")
include("Jastrow.jl")
include("VMC.jl")
include("Utilities.jl")

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

# Hubbard-U
U = 0.5

# onsite energy
# ε = 0.0

# chemical potential
μ = 3.0

# fugacity
# μₚₕ = 0.0
# opt_μₚₕ = true

# TODO: read in initial variational parameter set
# readin_vpars = false
# path_to_vpars = /path/to/variational/parameters/

######################################
##      VARIATIONAL PARAMTERS       ##
######################################
# μ: chemical potential
# μₚₕ: fugacity
# Δs: s-wave pairing
# Δd: d-wave pairing 
# Δa: antiferromagnetic (Neél) order
# Δc: uniform charge order
# Δcm: charge modulation stripe
# Δsm: spin modulation stripe

# Δd + Δa: uniform d-wave
# Δcm + Δsm: stripe order

# parameters to be optimized
parameters_to_optimize = ["μ", "Δs"]
parameter_values = [μ, 0.3]

# TODO: read in initial variational parameter set
# readin_jpars = false
# path_to_jpars = /path/to/jastrow/parameters/

# initial electron density Jastrow parameter
vᵢⱼ = 0.5

# initial electron spin Jastrow parameter
wᵢⱼ = 0.5

# initial phonon density Jastrow parameter
# pᵢⱼ = 0.5

# initial electron-phonon density Jastrow parameter
# uᵢⱼ = 0.5


##################################
## DEFINE SIMULATION PARAMETERS ##
##################################

# whether model is particle-hole transformed (if so, automatically GCE)
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

# optimization rate
η = 0.0001

# whether to output to terminal during runtime
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

# initialize variational parameters
variational_parameters = initialize_variational_parameters(parameters_to_optimize, parameter_values)

# get particle numbers (use if initial density is specified)
(Np, Ne, nup, ndn) = get_particle_numbers(n̄)

# get particle density (use if initial particle number if specified)
# density, Np, Ne = get_particle_density(nup, ndn)
    
###########################
## SET-UP VMC SIMULATION ##
###########################

# construct mean-field Hamiltonian and return variational operators
(H_mf, V) = build_mean_field_hamiltonian()

# initialize Slater determinant state and initial particle configuration
(D, pconfig, ε, ε₀, M, U) = build_slater_determinant()  

# initialize uncorrelated phonon state and initial particle configuration
# (P, phconfig) = build_phonon_state()

# initialize variational parameter matrices
As = get_Ak_matrices(V[1], U)
Aμ = get_Ak_matrices(V[2], U)

# initialize equal-time Green's function (W matrix)
W = get_equal_greens(M, D)

# construct density Jastrow factor
(init_dTvec, init_dpar_matrix, num_dpars) = build_jastrow_factor("density")

# construct spin Jastrow factor 
(init_sTvec, init_spar_matrix, num_spars) = build_jastrow_factor("spin")

# construct electron-phonon density Jastrow factor 
# (init_pTvec, init_ppar_matrix, num_ppars) = build_jastrow_factor("electron-phonon")


#############################
## INITIALIZE MEASUREMENTS ##
#############################

# Initialize standard tight binding model VMC measurements
initialize_measurements!(model_geometry, "tight binding")

# Initialze measurements related to the Hubbard model
initialize_measurements!(model_geometry, "hubbard")

# Initialze measurements related to electron-phonon models
# initialize_measurements!(model_geometry, "electron-phonon")

# Initialize density correlation measurements
initialize_correlation_measurements!(model_geometry, "density")

# Initialize spin correlation measurements
initialize_correlation_measurements!(model_geometry, "spin")


###################################
## BURNIN/THERMALIZATION UPDATES ##
###################################

# Iterate over burnin/thermalization updates.
for n in 1:N_burnin
    # perform local updates to electron dofs
    # electron_local_update!()

    # perform local updates to phonon dofs
    # phonon_local_updates
end

# recompute W and Tvec(s) for numerical stabilization

# Iterate over the number of bin, i.e. the number of time measurements will be dumped to file.
for bin in 1:N_bins

    # Iterate over the number of updates and measurements performed in the current bin.
    for n in 1:bin_size
        # perform local updates to electron dofs
        # electron_local_update!()

        # perform local updates to phonon dofs
        # phonon_local_update!()
    end

    # Write the average measurements for the current bin to file.
    # write_measurements!()
end




