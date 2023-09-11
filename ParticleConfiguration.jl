using LatticeUtilities
using Random


"""
    get_particle_numbers( density::{AbstractFloat} ) 

Returns total number of particles, total number of electrons, and number of 
spin-up and spin-down electrons for lattice with N sites.

"""
function get_particle_numbers(density)
    # @assert ce = true
    Ne = density * model_geometry.lattice.N
    @assert Ne % 2 == 0
    if pht == true
        # number of up, down electrons
        nup  = Ne / 2
        ndn = model_geometry.lattice.N - nup
        # total number of particles (electrons and holes)
        Np = nup + model_geometry.lattice.N - ndn
    else
        # number of up, down electrons
        nup = Ne / 2
        ndn = Ne / 2
        # total number of particles
        Np = Ne
        @assert Np % 2 == 0
    end
    return trunc(Int,Np), trunc(Int,Ne), trunc(Int,nup), trunc(Int,ndn)
end


"""
    get_particle_density( nup::Int, ndn::Int ) 

Returns particle density given the number of spin-up and spin-down electrons 
on a lattice with N sites, as well as total particle number and total number of electrons.

"""
function get_particle_density(nup, ndn)
    Ne = nup + ndn
    @assert Ne % 2 == 0
    Np = nup + model_geometry.lattice.N - ndn
    density = Ne / model_geometry.lattice.N
    return density, Np, Ne  
end


"""
    generate_initial_configuration() 

Returns a randomly generated initial configuration of electrons.

"""
function generate_initial_configuration()
    nbr_table = build_neighbor_table(bonds[1],
                                    model_geometry.unit_cell,
                                    model_geometry.lattice)
    init_pconfig = zeros(Int, length(nbr_table[1,:]))                                                      
    while sum(init_pconfig) < nup
        init_pconfig[rand(rng, 1:model_geometry.lattice.N)] = 1
    end
    while sum(init_pconfig) < Ne
        init_pconfig[rand(rng, model_geometry.lattice.N+1:2*model_geometry.lattice.N)] = 1
    end
    return init_pconfig
end


"""
    get_spindex_type( spindex::Int ) 

Returns the spin species at a given spindex.

"""
function get_spindex_type(spindex)
    @assert spindex < (2*model_geometry.lattice.N)+1
    if spindex < model_geometry.lattice.N+1
        return 1    # spin-up
    else
        return 2    # spin-down
    end
end


"""
    get_index_from_spindex( spindex::Int ) 

Returns the lattice site i for a given spindex.

"""
function get_index_from_spindex(spindex)
    @assert spindex < (2*model_geometry.lattice.N)+1
    if get_spindex_type(spindex) == 2
        return spindex - model_geometry.lattice.N
    else
        return spindex
    end
end


"""
    number_operator( site::Int, pconfig::Vector{Int} )

Returns the number of spin-up and spin-down electrons 
occupying a real lattice site i.  

"""
function number_operator(site, pconfig)
    return pconfig[site],                                           # number of spin-down electrons
         pconfig[site+model_geometry.lattice.N],                    # number of spin-down electrons
         pconfig[site] + pconfig[site+model_geometry.lattice.N]     # total number of electrons on a site
end


"""
    spin_operator( site::Int, pconfig::Vector{Int} )

Returns the total spin of particles occupying a real lattice site i. 
"""


"""
    get_particle_positions( pconfig::Vector{Int} )

Returns a dictionary of particle positions with keys and values,
"spindex" -> "lattice site".

"""
function get_particle_positions(pconfig)
    particle_positions = Dict()
    for i in eachindex(pconfig)
        if pconfig[i] == 1
            particle_positions[i] = get_index_from_spindex(i)
        else
        end
    end

    return sort(collect(particle_positions), by = x->x[1])
end


"""
    update_particlepos!( paritcle_positions, proposed_hop )

If a particle 'β' at site 'k' successfully hops to a neighboring site 'l', update its
position in 'particle_positions' as well as 'pconfig.

"""
function update_particle_position!(particle_positions, proposed_hop)
    particle_positions[proposed_hop[2]][1] = get_spindices_from_index(proposed_hop[5])[proposed_hop[3]]
    particle_positions[proposed_hop[2]][2] = proposed_hop[5]
    return particle_positions
end













#################################################### FUNCTION GRAVEYARD #########################################################

# "CREATION/DESTRUCTION OPERATORS"
# # TODO: 'site' actually refers to the 'spindex' here. Need to fix that...

# """

#     create_boson( site::Int, state::Matrix{Int} )

# Creates a boson on some lattice site.

# """
# function create_boson(site, state)
#     amp = sqrt(state[site] + 1)
#     state[site] += 1

#     return amp * state
# end

# """

#     destroy_boson( site::Int, state::Matrix{Int} )

# Annihilates a boson on some lattice site.

# """
# function destroy_boson(site, state)
#     if state[site] == 0
#         print("vacuum state!!")

#         return 0
#     else
#         amp = sqrt(state[site] - 1)
#         state[site] -= 1

#         return amp * state
#     end
# end



# """

#     create_fermion( site::Int, state::Matrix{Int}, PHT::Bool )

# Creates a fermion on some lattice site. If PHT is 'true', creates a spin down hole.

# """
# # if PHT == true, creates a spin-down hole at some site
# function create_fermion(site, state, PHT)
#     if PHT == true
#         if state[site] == 0
#             print("state occupied!!")

#             return 0
#         else
#             amp = state[site]
#             sgn = (-1)^(sum(state[1:site]))
#             state[site] -=1

#             return sgn * amp * state
#         end
#     else
#         if state[site] == 1
#             print("state occupied!!")
    
#             return 0
#         else
#             amp = 1 - state[state]
#             sgn = (-1)^(sum(state[1:site]))
#             state[site] += 1
    
#             return sgn * amp * state
#         end
#     end
# end

# """

#     destroy_fermion( site::Int, state::Matrix{Int}, PHT::Bool )

# Annihilates a fermion on some lattice site. If PHT is 'true', annihilates a spin down hole.

# """
# function destroy_fermion(site, state, PHT)
#     if PHT == true
#         if state[site] == 1
#             print("vacuum state!!")

#             return 0
#         else
#             amp = 1 - state[state]
#             sgn = (-1)^(sum(state[1:site]))
#             state[site] += 1

#             return sgn * amp * state
#         end
#     else
#         if state[site] == 0
#             print("vacuum state!!")
    
#             return 0
#         else
#             amp = state[site]
#             sgn = (-1)^(sum(state[1:site]))
#             state[site] -=1
    
#             return sgn * amp * state
#         end
#     end
# end