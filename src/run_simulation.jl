using LatticeUtilities
using Random
using LinearAlgebra
using DelimitedFiles
using BenchmarkTools
using Profile

# files to include
include("Hamiltonian.jl")
include("ParticleConfiguration.jl")
include("Jastrow.jl")
include("Markov.jl")
include("Utilities.jl")
include("Greens.jl")
include("StochasticReconfiguration.jl")
include("Measurements.jl")

#############################
## DEFINE MODEL PARAMETERS ##
#############################

# define the size of the lattice
Lx = 4
Ly = 4

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

# (BCS) chemical potential
μ_BCS = 3.0

# phonon fugacity
# μₚₕ = 0.0

#######################################
##      VARIATIONAL PARAMETERS       ##
#######################################

# TBD: how this will be done is a file will be generated with the headers 
#      being the parameter name followed by values. For non-uniform parameters,
#      it's indices will also be reported. 

# whether to read-in initial determinantal parameters
readin_detpars = false

# whether to read-in initial Jastrow parameters
readin_jpars = false    

# filepath
path_to_vpars = "/path/to/variational/parameters/"

# parameters to be optimized and initial value
parameters_to_optimize = ["Δs", "μ_BCS"]        # BCS wavefunction
parameter_values = [0.3, μ_BCS]                 # pht = true

# parameters to be optimized and initial value
# parameters_to_optimize = ["Δcs", "Δss"]       # stripe wavefunction
# parameter_values = [0.3, 0.3]                 # pht = false

# parameters to be optimized and initial value
# parameters_to_optimize = ["Δd", "Δa"]         # uniform d-wave wavefunction
# parameter_values = [0.3, 0.3]                 # pht = true

# parameters to be optimized and initial value
# parameters_to_optimize = ["Δa"]               # AFM (Neél) wavefunction
# parameter_values = [0.3]                      # pht = false

# parameters to be optimized and initial value
# parameters_to_optimize = ["Δc"]               # CDW wavefunction
# parameter_values = [0.3]                      # pht = false

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

# Maximum allowed error in the equal-time Green's function
# which is corrected by numerical stabilization
δW = 1e-3

# Maximum allowed error in the T vector
# which is corrected by numerical stabilization
δT = 1e-3

# SR stabilization factor
η = 1e-4      # 10⁻⁴ is probably good good for the Hubbard model

# initial SR optimization rate
dt = 0.1        # dt must be of sufficient size such that convergence is rapid and the algorithm remains stable

# whether to output to terminal during runtime
verbose = true

# debugging (this will be removed later)
debug = false

# whether to output matrices to file
write = false

# initialize addition simulation information dictionary
additional_info = Dict(
    "δW" => δW,
    "δT" => δT,
    "time" => 0.0,
    "N_burnin" => N_burnin,
    "N_updates" => N_updates,
    "N_bins" => N_bins,
    "bin_size" => bin_size,
    "fermionic_local_acceptance_rate" => 0.0,
    "initial_dt" => dt,
    "final_dt" => 0.0,
    "seed" => seed,
    "n_bar" => n̄,
    "global_energy" => 0.0,
    "μ_BCS" => 0.0,
    "Δs" => 0.0
)

##################
## DEFINE MODEL ##
##################

# square unit cell
unit_cell = UnitCell([[1.0,0.0], [0.0,1.0]],           # lattice vectors
                               [[0.0,0.0]])            # basis vectors 

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

# define model geometry
model_geometry = ModelGeometry(unit_cell,lattice, bonds)

# define non-interacting tight binding model
tight_binding_model = TightBindingModel([t,tp],μ_BCS)

# initialize determinantal parameters
determinantal_parameters = initialize_determinantal_parameters(parameters_to_optimize, parameter_values)

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
(D, pconfig, ε, ε₀, M, P) = build_determinantal_state()  

# initialize uncorrelated phonon state and initial particle configuration
# (P, phconfig) = build_phonon_state()

# initialize variational parameter matrices
A = get_Ak_matrices(V, P, ε, model_geometry)

# initialize equal-time Green's function (W matrix)
W = get_equal_greens(M, D)

# construct electron density-density Jastrow factor
e_den_den_jastrow = build_jastrow_factor("e-den-den")

# construct electron spin-spin Jastrow factor 
# e_spn_spn_jastrow = build_jastrow_factor("e-spn-spn")

# construct electron-phonon density-density Jastrow factor 
# eph_den_den_jastrow = build_jastrow_factor("eph-den-den")

# construct electron-phonon density-displacement Jastrow factor
# eph_den_dsp_jastrow = build_jastrow_factor("eph-den-dsp")

# construct phonon displacement-displacement Jastrow factor
# ph_dps_dsp_jastrow = build_jastrow_factor("eph-dsp-dsp")

# initialize variational parameters
variational_parameters = cat_vpars(determinantal_parameters, e_den_den_jastrow)


#############################
## INITIALIZE MEASUREMENTS ##
#############################

# initialize measurement container for VMC measurements
measurement_container = initialize_measurement_container(model_geometry, N_burnin, N_updates, variational_parameters)

# initialize energy measurements
initialize_measurements!(measurement_container, "energy")

# initialize correlation measurements
# initialize_correlation_measurements!(measurement_container, "density")


###########################################
## PERFORM BURNIN/THERMALIZATION UPDATES ##
###########################################

# start time for simulation
t_start = time()
if verbose
    println("|| START OF VMC SIMULATION ||")
end

# Iterate over burnin/thermalization updates.
for n in 1:N_burnin
    # perform local update to fermionic dofs
    (acceptance_rate, pconfig, jastrow, W) = local_fermion_update!(Ne, model_geometry, tight_binding_model, e_den_den_jastrow, pconfig, rng)

    # record acceptance rate
    # additional_info["fermionic_local_acceptance_rate"] += acceptance_rate

    # perform local updates to phonon dofs
    # local_phonon_update!(model_geometry, electron_phonon_model, jastrow, phconfig)

    # additional_info["phononic_local_acceptance_rate"] += acceptance_rate
end

# recompute W and Tvec(s) for numerical stabilization
# TODO: this may be moved within the updating scheme
(W, ΔW) = recalc_equal_greens(W, δW)
# recalc_Tvec(Tᵤ, δT, model_geometry)       # TODO: need to update recalc_Tvec() method


##################################################################
## PERFORM SIMULATION/MEASUREMENT UPDATES AND MAKE MEASUREMENTS ##
##################################################################

# Iterate over the number of bins, i.e. the number of measurements will be dumped to file.
for bin in 1:N_bins

    # Iterate over the number of updates and measurements performed in the current bin.
    for n in 1:bin_size
        # perform local update to fermionic dofs
        (acceptance_rate, pconfig, jastrow, W) = local_fermion_update!(Ne, model_geometry, tight_binding_model, e_den_den_jastrow, pconfig, rng)

        # record acceptance rate
        additional_info["fermionic_local_acceptance_rate"] += acceptance_rate

        # perform local updates to phonon dofs
        # local_phonon_update!(model_geometry, electron_phonon_model, jastrow, phconfig)

        # additional_info["phononic_local_acceptance_rate"] += acceptance_rate

        #TODO: add additional numerical stabilization?
    end

    # Write the average measurements for the current bin to file.
    write_measurements!(
            measurement_container = measurement_container,
            model_geometry = model_geometry,
            bin = bin,
            bin_size = bin_size
    )
end

# end time for simulation
t_end = time()
if verbose
    println("|| END OF VMC SIMULATION ||")
end

# record simulation runtime
additional_info["time"] += t_end - t_start

# normalize acceptance rate measurements
additional_info["fermionic_local_acceptance_rate"] /= (N_updates + N_burnin)
# additional_info["phonon_acceptance_rate"] /= (N_updates + N_burnin)

# write simulation information to file
# save_simulation_info(simulation_info, additional_info)

# process measurements
# process_measurements(simulation_info.datafolder, 20)



