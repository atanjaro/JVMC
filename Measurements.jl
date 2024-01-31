using LatticeUtilities
using LinearAlgebra
using Test


function initialize_measurement_container(model_geometry)
    lattice = model_geometry.lattice
    unit_cell = model_geometry.unit_cell

    # number of orbitals per unit cell
    n_orbitals = unit_cell.n

    # extent of the lattice in each direction
    L = lattice.L

    # initialize global measurements
    global_measurements = Dict{String, AbstractFloat}(
        "density" => zero(AbstractFloat),       # average total density ⟨n⟩
        "double_occ" => zero(AbstractFloat),    # double occupancy   
        "energy" => zero(AbstractFloat),        # global energy
        "sgn" => zero(AbstractFloat)            # fermion sign
    )

    # initialize local measurement
    local_measurements = Dict{String, Vector{AbstractFloat}}(
        "density" => zeros(AsbtractFloat, norbitals),       # average density for each orbital species
        "double_occ" => zeros(AbstractFloat, norbitals)     # average double occupancy for each orbital
    )

    # initialize measurement container
    measurement_container = (
        global_measurements         = global_measurements,
        local_measurements          = local_measurements,
        equaltime_correlations      = equaltime_correlations,
        L                           = L,
    )

    return measurement_container
end


function initialize_measurements!(measurement_container::NamedTuple,
    tight_binding_model::TightBindingModel{T,E}) where {T<:Number, E<:AbstractFloat}

(; local_measurements, global_measurements) = measurement_container

# number of orbitals per unit cell
norbital = length(tight_binding_model.ϵ_mean)

# number of types of hoppings
nhopping = length(tight_binding_model.t_bond_ids)

# initialize chemical potential as global measurement
global_measurements["chemical_potential"] = zero(Complex{E})

# initialize on-site energy measurement
local_measurements["onsite_energy"] = zeros(Complex{E}, norbital)

# initialize hopping energy measurement
if nhopping > 0
local_measurements["hopping_energy"] = zeros(Complex{E}, nhopping)
end

return nothing
end



"""
    initialize_correlation_measurments!(model_geometry::ModelGeometry, model::AbstractString )

Initialize either density correlation or spin correlation measurements.

"""
function initialize_correlation_measurements!(model_geometry, correlation_type)
    if correlation_type == "density"

    elseif correlation_type == "spin"

    else
    end
end


"""
    measure_double_occ( model_geometry::ModelGeometry, pconfig::Vector{Int} )

Measure the average double occupancy ⟨D⟩ = N⁻¹ ∑ᵢ ⟨nᵢ↑nᵢ↓⟩.

"""
function measure_double_occ(pconfig, model_geometry)
    nup_ndn = 0.0

    for i in 1:model_geometry.lattice.N
        nup_ndn += number_operator(i,pconfig)[1]*(1-number_operator(i,pconfig)[2])
    end
    
    return nup_ndn / model_geometry.lattice.N
end




"""
    measure_n( site::Int, pconfig::Vector{Int} )

Measure the local particle density ⟨n⟩.

"""
function measure_n(site)
    loc_den = number_operator(site,pconfig)[1] + 1 - number_operator(i,pconfig)[2]

    return loc_den
end


"""
    measure_ρ( site::int )

Measure the local excess hole density ⟨ρ⟩.

"""
function measure_ρ(site)
    ρ = 1 - measure_n(site)

    return ρ
end


"""
    measure_s( site::Int, pconfig::Vector{Int} )

Measure the local spin.

"""
function measure_s()
    loc_spn = number_operator(site,pconfig)[1] - 1 + number_operator(i,pconfig)[2]

    return loc_spn
end


"""
    measure_local_energy(model_geometry::ModelGeometry, tight_binding_model::TightBindingModel, jastrow::Jastrow, particle_positions:: )

Measure the local variational energy. Returns the total local energy,
local kinetic energy, and local Hubbard energy.

"""
function measure_local_energy(model_geometry, tight_binding_model, jastrow, particle_positions)
    nbr_table = build_neighbor_table(bonds[1],
                                    model_geometry.unit_cell,
                                    model_geometry.lattice)
    E_loc_kinetic = 0.0

    # loop over different electrons k
    for β in 1:Np
        k = particle_positions[β][2]
        # TBA: loop over different neighbor orders (i.e. next nearest neighbors)
        # loop over nearest neighbors
        sum_nn = 0.0
        for (i,j) in eachcol(nbr_table) # TBA: only loop over the known neighbors of β, l
            # reverse sign if system is particle-hole transformed
            if pht == true
                Tₗ = jastrow.Tvec[l]
                Tₖ = jastrow.Tvec[k]
                Rⱼ = exp(-get_jastrow_ratio(l, k, Tₗ, Tₖ))
            else
                Tₗ = jastrow.Tvec[l]
                Tₖ = jastrow.Tvec[k]
                Rⱼ = exp(get_jastrow_ratio(l, k, Tₗ, Tₖ))
            end
            sum_nn += Rⱼ * W[l, β]
        end
        # reverse sign if system is particle-hole transformed
        if pht == true
            E_loc_kinetic += tight_binding_model.t[1] * sum_nn          
        else
            E_loc_kinetic += - tight_binding_model.t[1] * sum_nn
        end
    end

    # calculate Hubbard energy
    E_loc_hubb = U * number_operator(site,pconfig)[1] *(1 - number_operator(i,pconfig)[2])

    # resultant local energy
    E_loc = E_loc_kinetic + E_loc_hubb

    return E_loc, E_loc_kinetic, E_loc_hubb
end



"""
    measure_global_energy( model_geometry::ModelGeometry )

Measure the global variational energy ⟨E⟩.

"""
function measure_global_energy(model_geometry, N_bins, bin_size)

    # binned energies
    # E_binned

    # average over all sites
    E_global = E_binned / model_geometry.lattice.N

    return E_global
end


"""
    measure_density_corr( )

Measure the density-density correlation function.

"""
function measure_density_corr()
    return nothing
end


"""
    measure_spin_corr( )

Measure the spin-spin correlation function.

"""
function measure_spin_corr()
    return nothing
end


# """
#     measure_average_X( )

# Measure average phonon displacement ⟨X⟩, in the SSH model (TBD).

# """
# function measure_average_X()
#     return nothing
# end


