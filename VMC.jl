using LinearAlgebra
using ForwardDiff

# Define your matrix A and the parameter x
A = [1.0 2.0; 3.0 4.0]
x = 2.0

# Define a function that computes the matrix logarithm of A
function logA_function(x)
    return log(A) 
end

# Compute the derivative of the logA_function with respect to x using ForwardDiff
dlogA_dx = ForwardDiff.derivative(logA_function, x)

println("Derivative of logA with respect to x:")
println(dlogA_dx)


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
    if proposed_hop[1] == true
        # HOP!
        k_spindex = get_spindices_from_index(proposed_hop[4])[proposed_hop[3]]
        l_spindex = get_spindices_from_index(proposed_hop[5])[proposed_hop[3]]

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

Perform a local update.

"""
function local_update!()
    particle_positions = get_particle_positions(pconfig)
    proposed_hop = propose_random_hop(particle_positions)
    pconfig = do_particle_hop(pconfig, proposed_hop)
    particle_positions = update_particle_position!(particle_positions, proposed_hop)
    
    return nothing
end