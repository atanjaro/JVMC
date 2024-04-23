using LatticeUtilities
using LinearAlgebra
using Test



"""
    measure_local_jpar_derivative( jastrow::Jastrow, pconfig::Vector{Int} )

Performs local logarithmic derivative Δₖ(x) = ∂lnΨ(x)/∂vₗₘ, with respect to the kth Jastrow parameter vₗₘ. Returns 
a vector of derivatives.

"""
function measure_local_jpar_derivative(jastrow, pconfig)

    # number of Jastrow parameters
    num_jpars = jastrow.num_jpars

    # map of Jastrow parameters
    jpar_map = jastrow.jpar_map

    # vector to store derivatives
    derivatives = zeros(AbstractFloat, num_jpars)
                
    # for density Jastrow
    if jastrow.jastrow_type == "density"
        for num in 1:jastrow.num_jpars
            for (i, r1) in jpar_map
                for (j, r2) in jpar_map
                    if r1[1] == r2[1]
                        if pht == true
                            derivatives[num] += -(number_operator(j[1],pconfig)[1] - number_operator(j[1],pconfig)[2]) * (
                                                  number_operator(j[2],pconfig)[1] - number_operator(j[2],pconfig)[2])
                        elseif pht == false
                            derivatives[num] += -(number_operator(j[1],pconfig)[1] + number_operator(j[1],pconfig)[2]) * (
                                                  number_operator(j[2],pconfig)[1] + number_operator(j[2],pconfig)[2])
                        else
                        end
                    else
                    end
                end
            end
        end

        return derivatives

    # for spin Jastrow
    elseif jastrow.jastrow_type == "spin"   
        for num in 1:jastrow.num_jpars
            for (i, r1) in jpar_map
                for (j, r2) in jpar_map
                    if r1[1] == r2[1]
                        if pht == true
                            derivatives[num] += -(number_operator(j[1],pconfig)[1] - number_operator(j[1], pconfig)[2]) * (
                                                  number_operator(j[2], pconfig)[1] - number_operator(j[2], pconfig)[2])
                        elseif pht == false
                            derivatives[num] += -(number_operator(j[1],pconfig)[1] + number_operator(j[1], pconfig)[2]) * (
                                                  number_operator(j[2], pconfig)[1] + number_operator(j[2], pconfig)[2])
                        else
                        end
                    end
                end
            end
        end

        return derivatives

    # for electron-phonon Jastrow
    elseif jastrow.jastrow_type == "electron-phonon"   
        return derivatives
    else
    end
end


"""
    measure_local_detpar_derivative( determinantal_parameters::DeterminantalParameters, model_geometry::ModelGeometry
                                     pconfig::Vector{Int}, W::Matrix{AbstractFloat}, A::Matrix{AbstractFloat}  )

Performs local logarithmic derivative Δₖ(x) = ∂lnΨ(x)/∂αₖ, with respect to the kth variational parameter αₖ,
in the determinantal part of the wavefunction. Returns a vector of derivatives.

"""
function measure_local_detpar_derivative(determinantal_parameters, model_geometry, pconfig, Np, W, A)  

    # dimensions
    dims = model_geometry.unit_cell.n * model_geometry.lattice.N

    # number of determinantal parameters
    num_detpars = determinantal_parameters.num_detpars
    
    # particle positions
    particle_positions = get_particle_positions(pconfig)

    # vector to store derivatives
    derivatives = zeros(AbstractFloat, num_detpars)
    

    # loop over Nₚ particles # TBA: need to fix method
    G = zeros(AbstractFloat, 2*dims, 2*dims)
    for β in 1:Np
        for j in 1:2*dims
            for (spindex, iᵦ) in particle_positions
                G[iᵦ,j] = W[j,β]
            # G[iᵦ,:] = W[:,β]
            end
        end
    end

    # loop over the number of determinantal parameters
    for num in 1:num_detpars
        derivatives[num] += sum(A[num] * G)
    end

    return derivatives
end


"""
    measure_double_occ( model_geometry::ModelGeometry, pconfig::Vector{Int} )

Measure the average double occupancy ⟨D⟩ = N⁻¹ ∑ᵢ ⟨nᵢ↑nᵢ↓⟩.

"""
function measure_double_occ(pconfig, model_geometry)
    nup_ndn = 0.0

    for i in 1:model_geometry.lattice.N
        nup_ndn += number_operator(i, pconfig)[1] * (1 - number_operator(i, pconfig)[2])
    end
    
    return nup_ndn / model_geometry.lattice.N
end




"""
    measure_n( site::Int, pconfig::Vector{Int} )

Measure the local particle density ⟨n⟩.

"""
function measure_n(site)
    loc_den = number_operator(site, pconfig)[1] + 1 - number_operator(i, pconfig)[2]

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

    # generate neighbor table
    nbr_table = build_neighbor_table(bonds[1],
                                    model_geometry.unit_cell,
                                    model_geometry.lattice)

    # gnerate neighbor map
    nbr_map = map_neighbor_table(nbr_table)

    # loop over different electrons k
    E_loc_kinetic = 0.0
    for β in 1:Np
        k = particle_positions[β][2] 
        # loop over nearest neighbors. TBA: loop over different neighbor orders (i.e. nearest and next nearest neighbors)
        sum_nn = 0.0
        for l in nbr_map[k][2]
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
    E_loc_hubb = U * number_operator(site, pconfig)[1] * (1 - number_operator(i, pconfig)[2])

    # calculate total local energy
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


