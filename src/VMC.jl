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
   acceptance, particle, spin, isite, fsite =  metropolis(particle_positions)

   return LocalAcceptance(acceptance, particle, spin, isite, fsite)
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

        # begin Metropolis algorithm

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
    do_particle_hop!( local_acceptance::LocalAcceptance, pconfig::Matrix{Int})

If proposed particle hop is accepted, perform the particle hop.

"""
function do_particle_hop!(proposed_hop, pconfig)
    if proposed_hop.acceptance == true
        k = proposed_hop.isite
        l = proposed_hop.fsite

        # HOP!
        pconfig[k] = 0
        pconfig[l] = 1

        return pconfig
    else
        return pconfig
    end
end


"""
    local_update!()

Perform a local MC update. Proposes move and accept rejects via Metropolis algorithm,
if accepted, updates particle positions, T vector, W matrix, and variational parameters.

"""
function local_update!()

    # accept/reject (Metropolis)
    proposed_hop = metropolis(particle_positions)

    # perform hopping
    pconfig = do_particle_hop!(proposed_hop, pconfig)

    # update particle positions
    update_particle_position!(particle_positions, proposed_hop)

    # update T vector
    T_vec = update_Tvec!(LocalAcceptance.fsite, LocalAcceptance.isite, Tvec)

    # update Jastrow pseudopotentials
     # need to account for number of parameters being optimized
    # need to update value of parameter, as well as keep track of it's history over the N_updates

    # update Green's function
    W = update_equal_greens!(local_acceptance)

    # update variational parameters
    vpar = parameter_update!()
    # need to account for number of parameters being optimized
    # need to update value of parameter, as well as keep track of it's history over the N_updates
    # perhaps store the histories in a num_vaprs by N_updates matrix?
    
    return acceptance_rate, W, Tvec, vpar
end


