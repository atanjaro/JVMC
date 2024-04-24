using Test

include("../src/Hamiltonian.jl")
include("../src/ParticleConfiguration.jl")
include("../src/Jastrow.jl")
include("../src/VMC.jl")
include("../src/Utilities.jl")
include("../src/Greens.jl")
include("../src/StochasticReconfiguration.jl")
include("../src/Measurements.jl")

#  begin unit test for a 2x2 square lattice Hubbard model with s-wave pairing
U = 1.0
unit_cell = UnitCell([[1.0,0.0], [0.0,1.0]],[[0.0,0.0]])
lattice = Lattice([2,2],[true,true])  
model_geometry = ModelGeometry(unit_cell,lattice)        
bonds = [[Bond((1,1), [1,0]), Bond((1,1), [0,1])]]
pht = true
tight_binding_model = TightBindingModel([1.0,0],3.0)
determinantal_parameters = initialize_determinantal_parameters(["Δs", "μ"], [0.3, 3.0])
n = 1.0
nup = 2
ndn = 2

# check that particles in the canonical ensemble are correctly returned
function test_canonical_ensemble()
    @test get_particle_numbers(n) == (4, 4, 2, 2)
    @test get_particle_density(nup,ndn) == (1.0, 4, 4)
    @test_throws AssertionError get_particle_numbers(0.8)  
    @test_throws AssertionError get_particle_density(2,3)
end

# checks mean-field Hamiltonian and determinantal state
function test_hamiltonian_build()
    bmf = build_mean_field_hamiltonian()
    ϕ = build_determinantal_state() 
    
    @test bmf[1] == bmf[1]'
    @test diagonalize(bmf[1])[1] == ϕ[3]
    @test length(ϕ[1])%length(ϕ[4]) == 0
end
