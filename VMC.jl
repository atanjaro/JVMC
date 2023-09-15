using LinearAlgebra

include("ParticleConfiguration.jl")

"""
    LocalAcceptance( acceptance::Bool, particle::Int, )

A type defining quantities related to local MC update acceptance

"""
struct LocalAcceptance
    # whether move is accepted
    acceptance::Bool
    # selected particle to be moved
    particle::Int
    # the selected particle's spin
    spin::Int
    # initial site
    isite::Int
    # final site
    fsite::Int
end


"""
    local_acceptance()

Constructor for the local acceptance type.

"""
function local_acceptance()
   acceptance, particle, spin, isite, fsite =  propose_random_hop(particle_positions)

   return LocalAcceptance(acceptance, particle, spin, isite, fsite)
end


"""
    local_jastrow_derivative(jpar_indices::CartesianIndex{2}, pconfig::Vector{AbstractFloat, jastrow_type::AbstractString)

Performs the kth local logarithmic derivative Δₖ(x) = ∂lnΨ(x)/∂vₗₘ, with respect
to the kth Jastrow parameter vₗₘ.

"""
function local_jastrow_derivative(jpar_indices, pconfig, jastrow_type)
    if jastrow_type == "density"
        if pht == true
            Δₖ = -(number_operator(jpar_indices[1],pconfig)[1] - number_operator(jpar_indices[1],pconfig)[2])*(number_operator(jpar_indices[2],pconfig)[1] - number_operator(jpar_indices[2],pconfig)[2])
        else
            Δₖ = -(number_operator(jpar_indices[1],pconfig)[1] + number_operator(jpar_indices[1],pconfig)[2])*(number_operator(jpar_indices[2],pconfig)[1] + number_operator(jpar_indices[2],pconfig)[2])
        end
    elseif jastrow_type == "spin"
        if pht == true
            Δₖ = -0.5 * (number_operator(jpar_indices[1],pconfig)[1] + number_operator(jpar_indices[1],pconfig)[2])*(number_operator(jpar_indices[2],pconfig)[1] + number_operator(jpar_indices[2],pconfig)[2])
        else
            Δₖ = -0.5 * (number_operator(jpar_indices[1],pconfig)[1] - number_operator(jpar_indices[1],pconfig)[2])*(number_operator(jpar_indices[2],pconfig)[1] - number_operator(jpar_indices[2],pconfig)[2])
        end
    elseif jastrow_type == "electron-phonon"
        # derivative of electron-phonon Jastrow factor
    else
    end
    return Δₖ
end



"""
    local_slater_derivative(vpar_indices::CartesianIndex{2}, A::Matrix{AbstractFloat}, W::Matrix{AbstractFloat})

Performs the kth local logarithmic derivative Δₖ(x) = ∂lnΨ(x)/∂αₖ, with respect
to the kth variational parameter αₖ, in the determinantal part of the wavefunction.

"""
function local_slater_derivative(vpar_indices, A, W)
    # sum_over_number_of_particles(i = 1 to Nₚ)sum_over_spindices(j = 1 to 2L)
    # Δₖ = ∑ᵢ₌₁ᴺ ∑ⱼ₌₁²ᴸ (Aₖ)ᵢᵦⱼWⱼᵦ
    return nothing
end



"""
    propose_random_hop( particle_positions::Vector{Dict{Any,Any}} )

Randomly selects a particle 'β' at site 'k' to hop to a neighboring site 'l'.

"""
function propose_random_hop(particle_positions)
    nbr_table = build_neighbor_table(bonds[1],
                                    model_geometry.unit_cell,
                                    model_geometry.lattice)
    nbrs = map_neighbor_table(nbr_table)
    beta = rand(rng, 1:trunc(Int,Np))                   # randomly select some particle in the lattice
    k = particle_positions[beta][2]                # real position 'k' of particle 'β' 
    l = rand(rng, 1:nbrs[k][2][2])                        # random neighboring site 'l'
    beta_spin = get_spindex_type(particle_positions[beta][1])
    
    # checks occupation against spin species of particle 'β'
    # if site is unoccupied by same spin species, hop is accepted
    if number_operator(l,pconfig)[beta_spin] == 1
        if verbose == true
            println("HOP REJECTED!")  
        end
        return false, beta, beta_spin, k, l # acceptance, particle number, particle spin, initial site, final site
    else
        if verbose == true
            println("HOP ACCEPTED!")
        end
        return true, beta, beta_spin, k, l  
    end
end

"""
    do_particle_hop!( pconfig::Matrix{Int}, proposed_hop:: )

If proposed particle hop is accepted, perform the particle hop.

"""
function do_particle_hop!(pconfig, proposed_hop)
    if proposed_hop.acceptance == true
        # HOP!
        k_spindex = get_spindices_from_index(proposed_hop.isite)[proposed_hop.spin]
        l_spindex = get_spindices_from_index(proposed_hop.fsite)[proposed_hop.spin]

        @assert pconfig[k_spindex] == 1
        pconfig[k_spindex] = 0
        @assert pconfig[l_spindex] == 0
        pconfig[l_spindex] = 1
        
        return pconfig 
    else 
        # DO NOTHING
        return nothing
    end
end


"""
    local_update!()

Perform a local MC update.

"""
function local_update!()
    particle_positions = get_particle_positions(pconfig)
    proposed_hop = local_acceptance()
    accepted = proposed_hop.acceptance

    pconfig = do_particle_hop(pconfig, proposed_hop)
    particle_positions = update_particle_position!(particle_positions, proposed_hop)

    # perform stochastic reconfiguration

    # measure local energy
    
    return accepted, pconfig, particle_positions #, local_energy
end


"""
    parameter_gradient(vpar)

Perform gradient descent on variational parameters for Stochastic 
Reconfiguration.

"""
function parameter_gradient(vpar_matrices, jpar_matrix, pconfig)
    # perform gradient descent on Jastrow parameters
    num_jpars = get_num_jpars(jpar_matrix)
    jpar_indices = get_parameter_indices(jpar_matrix)
    for k in 1:2*num_jpars
        Δₖ = local_jastrow_derivative(jpar_indices[k],pconfig)
    end
    # perform gradient descent on variational parameters

    return nothing
end


"""
    parameter_update!()

Update variational parameters.

"""
function parameter_update!()
    return nothing
end



"""
    parameter_indices()

Get indices of variational parameters from its respective matrix.

"""
function get_parameter_indices(par_matrix)
    nonzero_indices = findall(x -> x != 0, par_matrix)

    parameter_indices = sort(nonzero_indices, by=x->(x[1], x[2]))

    return parameter_indices
end ## TODO: move to Jastrow.jl?

