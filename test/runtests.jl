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

# check that particle numbers are correctly returned
function test_particle_numbers()
    @test get_particle_numbers(1.0) == (4, 4, 2, 2)
    @test get_particle_numbers(0.5) == (4, 2, 1, 1)
    @test_throws AssertionError get_particle_numbers(0.8)  
end

# check that particle densities are correctly returned
function test_particle_density()
    @test get_particle_density(2,2) == (1.0, 4, 4)
    @test get_particle_density(1,1) == (0.5, 4, 2)
    @test_throws AssertionError get_particle_density(2,3)
end


function test_hamiltonian_build()
    @test build_mean_field_hamiltonian() == [-3.0  -1.0  -1.0   0.0  0.3  0.0  0.0  0.0; -1.0  -3.0   0.0  -1.0  0.0  0.3  0.0  0.0; -1.0   0.0  -3.0  -1.0  0.0  0.0  0.3  0.0; 0.0  -1.0  -1.0  -3.0  0.0  0.0  0.0  0.3;  0.3   0.0   0.0   0.0  3.0  1.0  1.0  0.0]
end


function test_diagonalization()

end