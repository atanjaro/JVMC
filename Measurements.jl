using LatticeUtilities
using LinearAlgebra
using Test


"""
    initialize_measurments!(model_geometry::ModelGeometry, model:AbstractString )

Initialize measurements related to different models. Standard VMC measurements (global
energy, hopping energy, onsite energy, double occupancy, and local density) average
initialized using the "tight binding" model.

"""
function initialize_measurements!(model_geometry, model)
    lattice = model_geometry.lattice
    unit_cell = model_geometry.unit_cell
    if model == "tight binding"     # contains all the standard VMC measurements
        global_measurements = Dict{String, Complex{T}}(
            "energy" => zero(Complex{T}),
            "density" => zero(Complex{T}), # average total density ⟨n⟩
            "double_occ" => zero(Complex{T}),
            "Nsqrd" => zero(Complex{T}), # total particle number square ⟨N²⟩
            "sgn" => zero(Complex{T}), # sign(det(Gup))⋅sign(det(Gdn))
        )

        local_measurements = Dict{String, Vector{Complex{T}}}(
        "density" => zeros(Complex{T}, norbitals), # average density for each orbital species
        "double_occ" => zeros(Complex{T}, norbitals), # average double occupancy for each orbital
    )
        # measure global energy
        E_global = []
        # measure hopping energy
        E_hop = []
        # measure onsite energy
        E_onsite = []
        # measure double occupancy
        dbocc = []
        # measure local density
        ldens = []
    elseif model == "hubbard"
        # measure Hubbard energy
        E_U = []
        # measure local hole density
        hdens = []
        # measure local spin
        lspin = []
    elseif model == "electron-phonon"
        # measure average phonon position
        dispX = []
        # measure Holstein energy
        E_holst = []
        # measure SSH energy
        E_ssh = []
        # measure SSH sign switching
        sgn_sw = []
    else      
    end
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
    measure_double_occ( pconfig::Vector{Int} )

Measure the average double occupancy ⟨D⟩ = N⁻¹ ∑ᵢ ⟨nᵢ↑nᵢ↓⟩.

"""
function measure_double_occ(pconfig)
    nup_ndn = 0.0

    for i in 1:model_geometry.N
        nup_ndn += number_operator(i,pconfig)[1]*(1-number_operator(i,pconfig)[2])
    end
    
    return nup_ndn
end




"""
    measure_n( site::Int )

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
    return 1 - measure_n(site)
end


"""
    measure_s( site::Int )

Measure the local spin.

"""
function measure_s()
    loc_spn = number_operator(site,pconfig)[1] - 1 + number_operator(i,pconfig)[2]

    return loc_spn
end


"""
    measure_local_energy( )

Measure the local variational energy. Returns the total local energy,
local kinetic energy, and local Hubbard energy.

"""
function measure_local_energy()

    # calculate expectation value of kinetic energy
    E_loc_kinetic = 0.0

    # loop over different electrons k

    # loop over nearest neighbors

    # reverse sign if pht == true

    # calculate Hubbard energy
    E_loc_hubb = U * number_operator(site,pconfig)[1] *(1-number_operator(i,pconfig)[2])

    # resultant local energy
    E_loc = E_loc_kinetic + E_loc_hubb

    return E_loc, E_loc_kinetic, E_loc_hubb
endy



"""
    measure_global_energy( )

Measure the global variational energy ⟨E⟩.

"""
function measure_global_energy()

    # account for bins

    # average over all sites

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


