using LinearAlgebra
using Distributions

include("ParticleConfiguration.jl")

"""
    LocalAcceptance( acceptance::Bool, particle::Int, spin::Int, isite::Int, fsite::Int )

A type defining quantities related to local MC update acceptance

"""
struct LocalAcceptance
    # whether move is possible
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

Performs local logarithmic derivative Δₖ(x) = ∂lnΨ(x)/∂vₗₘ, with respect
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

Performs local logarithmic derivative Δₖ(x) = ∂lnΨ(x)/∂αₖ, with respect
to the kth variational parameter αₖ, in the determinantal part of the wavefunction.

"""
function local_slater_derivative(Ak, W, acceptance)
    iᵦ = acceptance.isite
    # j = acceptance.fsite

    Δₖ = 0.0

    # sum over number of particles
    for β in 1:Np
        # sum over lattice sites for the up and down sectors
        for j in 1:2*model_geometry.lattice.N
            Δₖ += Ak[iᵦ,j]*W[j,β]
        end
    end

    return Δₖ
end


## DEPRECATED
# """
#     propose_random_hop( particle_positions::Vector{Dict{Any,Any}} )

# Randomly selects a particle 'β' at site 'k' to hop to a neighboring site 'l'.

# """
# function propose_random_hop(particle_positions)
#     nbr_table = build_neighbor_table(bonds[1],
#                                     model_geometry.unit_cell,
#                                     model_geometry.lattice)
#     nbrs = map_neighbor_table(nbr_table)
#     beta = rand(rng, 1:trunc(Int,Np))                   # randomly select some particle in the lattice
#     k = particle_positions[beta][2]                # real position 'k' of particle 'β' 
#     l = rand(rng, 1:nbrs[k][2][2])                        # random neighboring site 'l'
#     beta_spin = get_spindex_type(particle_positions[beta][1])
    
#     # checks occupation against spin species of particle 'β'
#     # if site is unoccupied by same spin species, hop is possible
#     if number_operator(l,pconfig)[beta_spin] == 1
#         if verbose == true
#             println("HOP NOT POSSIBLE!")  
#         end
#         return false 
#     else
#         if verbose == true
#             println("HOP POSSIBLE!")
#         end
#         return true, beta, beta_spin, k, l  # acceptance, particle number, particle spin, initial site, final site
#     end
# end


"""
    metropolis( particle_positions::Vector{Dict{Any,Any}})

Perform accept/reject step of proposed hop using the Metropolis algorithm. If move is
accepted, returns acceptance, particle β and it's spindex, initial position, and final
position.

"""

function metropolis(particle_positions)
    nbr_table = build_neighbor_table(bonds[1],
                                    model_geometry.unit_cell,
                                    model_geometry.lattice)
    nbrs = map_neighbor_table(nbr_table)
    beta = rand(rng, 1:trunc(Int,Np))                   # randomly select some particle in the lattice
    k = particle_positions[beta][2]                # real position 'k' of particle 'β' 
    l = rand(rng, 1:nbrs[k][2][2])                        # random neighboring site 'l'
    beta_spin = get_spindex_type(particle_positions[beta][1])
    
    # checks occupation against spin species of particle 'β'
    # if site is unoccupied by same spin species, hop is possible
    if number_operator(l,pconfig)[beta_spin] == 1
        if verbose == true
            println("HOP NOT POSSIBLE!")  
        end
        return false 
    else
        if verbose == true
            println("HOP POSSIBLE!")
        end

        # Metropolis algorithm

        # get Jastrow ratio (element of T vector)
        Rⱼ = T[k,l]        # TODO: use views to obtain element?

        # get wavefunction ratio (correpsonding element of Green's function)
        Rₛ = W[k,l]        # TODO: use views to obtain element?

        acceptance_prob = Rⱼ * Rⱼ * Rₛ * Rₛ     

        if acceptance_prob > 1 || rand(rng,dist,1) < acceptance_prob
            if verbose == true
                println("HOP ACCEPTED!")
            end

            # do particle hop
            
            return true, beta, beta_spin, k, l  # acceptance, particle number, particle spin, initial site, final site
        else
            if verbose == false
               println("HOP REJECTED")
            end
        end
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

Perform a local MC update. Proposes move and accept rejects via Metropolis algorithm,
if accepted, updates particle positions, T vector, W matrix, and variational parameters.

"""
function local_update!()
    pconfig = do_particle_hop(pconfig, proposed_hop)
    particle_positions = update_particle_position!(particle_positions, proposed_hop)

    # accept/reject (Metropolis)
    proposed_hop = metropolis(particle_positions)

    # update paticle positions
    do_particle_hop!(pconfig, proposed_hop)

    # update T vector
    #Tvec[] = 

    # update Green's function
    #W[] = 

    # update parameters
    #vpars = []
    
    return nothing
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
        Δₖ = local_jastrow_derivative(jpar_indices[k],pconfig,type)
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

