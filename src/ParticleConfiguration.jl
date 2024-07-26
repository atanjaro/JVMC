using LatticeUtilities
using Random

"""
    get_particle_numbers( density::{AbstractFloat} ) 

Returns total number of particles, total number of electrons, and number of 
spin-up and spin-down electrons for lattice with N sites.

"""
function get_particle_numbers(density)
    Ne = density * model_geometry.lattice.N
    @assert Ne % 2 == 0 "Ne must be even"
    
    nup = Ne ÷ 2
    ndn = Ne ÷ 2
    
    if !pht
        Np = Ne
        @assert Np % 2 == 0 "Np must be even"
    else
        Np = nup + model_geometry.lattice.N - ndn
    end
    
    return Int(Np), Int(Ne), Int(nup), Int(ndn)
end

# function get_particle_numbers(density)
#     # @assert ce = true
#     Ne = density * model_geometry.lattice.N
#     @assert Ne % 2 == 0
#     if pht == true
#         # number of up, down electrons
#         nup  = Ne / 2
#         ndn = Ne - nup
#         # total number of particles (electrons and holes)
#         Np = nup + model_geometry.lattice.N - ndn
#     else
#         # number of up, down electrons
#         nup = Ne / 2
#         ndn = Ne / 2
#         # total number of particles
#         Np = Ne
#         @assert Np % 2 == 0
#     end
#     return trunc(Int,Np), trunc(Int,Ne), trunc(Int,nup), trunc(Int,ndn)
# end



"""
    get_particle_density( nup::Int, ndn::Int ) 

Returns particle density given the number of spin-up and spin-down electrons 
on a lattice with N sites, as well as total particle number and total number of electrons.

"""
function get_particle_density(nup, ndn)
    Ne = nup + ndn
    @assert Ne % 2 == 0 "Ne must be even"
    
    Np = nup + model_geometry.lattice.N - ndn
    density = Ne / model_geometry.lattice.N
    
    return density, Int(Np), Int(Ne)
end
# function get_particle_density(nup, ndn)
#     Ne = nup + ndn
#     @assert Ne % 2 == 0
#     Np = nup + model_geometry.lattice.N - ndn
#     density = Ne / model_geometry.lattice.N
#     return density, Np, Ne  
# end


"""
    generate_initial_electron_configuration() 

Returns a randomly generated initial configuration of electrons.

"""
function generate_initial_electron_configuration()
    init_pconfig = zeros(Int, 2*model_geometry.lattice.N)                                                      
    while sum(init_pconfig) < nup
        init_pconfig[rand(rng, 1:model_geometry.lattice.N)] = 1
    end
    while sum(init_pconfig) < Ne
        init_pconfig[rand(rng, model_geometry.lattice.N+1:2*model_geometry.lattice.N)] = 1
    end
    return init_pconfig
end


"""
    generate_initial_onsite_phonon_configuration() 

Returns a randomly generated initial configuration of onsite optical phonons.
This would apply to Holstein and optical-SSH models.

"""
function generate_initial_onsite_phonon_configuration()
    init_phconfig = zeros(Int, model_geometry.lattice.N)  

    for i in 1:model_geometry.lattice.N
        upper_bound = rand(rng, 1:typemax(Int))
        init_phconfig[i] = rand(rng, 0:upper_bound) 
    end

    return init_phconfig
end


"""
    generate_initial_bond_phonon_configuration() 

Returns a randomly generated initial configuration of bond optical phonons.
This would apply only to the bond-SSH model.

"""
function generate_initial_bond_phonon_configuration()
    init_phconfig = zeros(Int, 2 * model_geometry.lattice.N)  

    for i in 1:2 * model_geometry.lattice.N
        upper_bound = rand(rng, 1:typemax(Int))
        init_phconfig[i] = rand(rng, 0:upper_bound) 
    end

    return init_phconfig
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
    get_spindices_from_index( index::Int ) 

Returns spin-up and spin-down indices from a given site index.

"""
function get_spindices_from_index(index)
    @assert index <= model_geometry.lattice.N
    spindex1 = index
    spindex2 = index + model_geometry.lattice.N
    return spindex1, spindex2
end





"""
    number_operator( site::Int, pconfig::Vector{Int} )

Returns the number of spin-up and spin-down electrons 
occupying a real lattice site i.  

"""
function number_operator(site, pconfig)
    nup = pconfig[site]
    ndn = pconfig[site+model_geometry.lattice.N]
    Ne = pconfig[site] + pconfig[site+model_geometry.lattice.N]

    return nup, ndn, Ne
end

"""
    do_particle_hop!( pconfig::Matrix{Int})

If proposed particle hop is accepted, perform the particle hop.

"""
function do_particle_hop!(proposed_hop, pconfig)
    # spin of the particle
    spin = proposed_hop.spin

    # initial site
    k = proposed_hop.isite

    # final site
    l = proposed_hop.fsite

    # account for spin-up and spin-down sectors
    if spin == 2
        k_dn = get_spindices_from_index(k)[2]
        l_dn = get_spindices_from_index(l)[2]

        if debug
            @info "Hopping particle from site $k to site $l"
        end

        pconfig[k_dn] = 0
        pconfig[l_dn] = 1
    else
        pconfig[k] = 0
        pconfig[l] = 1 
    end
    
    return nothing
end



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

    sorted_positions = sort(collect(particle_positions), by = x->x[1])

    if length(sorted_positions) != Np
        @error "Mismatch in particle_positions length: expected $Np, got $(length(sorted_positions))"
    end
    
    return sorted_positions
end


"""
    update_particle_position!(proposed_hop , paritcle_positions)

If a particle 'β' at site 'k' successfully hops to a neighboring site 'l', update its
position in 'particle_positions' as well as 'pconfig.

"""
function update_particle_position!(proposed_hop, particle_positions)
    if proposed_hop.acceptance == 1
        # Update the pair within the vector
        particle_positions[proposed_hop.particle] = Pair(get_spindices_from_index(proposed_hop.fsite)[proposed_hop.spin], proposed_hop.fsite)
        return nothing
    else
        # DO NOTHING
        return nothing
    end
end
















