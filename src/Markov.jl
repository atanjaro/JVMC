using LinearAlgebra
using Distributions

include("ParticleConfiguration.jl")

"""
    LocalAcceptance( acceptance::Bool, particle::Int, spin::Int, isite::Int, fsite::Int )

A type defining quantities related to local MC update acceptance

"""
struct LocalAcceptance
    # whether move is possible
    acceptance::Int
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




"""
    metropolis( W, jastrow, particle_positions, rng )

Perform accept/reject step of proposed hop using the Metropolis algorithm. If move is
accepted, returns acceptance, particle β and it's spindex, initial position, and final
position.

"""

function metropolis(W, jastrow, particle_positions, rng)    
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
            println("Hop impossible!")  
        end
        return LocalAcceptance(0, beta, beta_spin, k, l)
    else
        if verbose == true
            println("Hop possible!")
        end

        # begin Metropolis algorithm

        # get Jastrow ratio (element of T vector)
        Rⱼ = get_jastrow_ratio(k, l, jastrow)    

        # get wavefunction ratio (correpsonding element of Green's function)
        Rₛ = W[l, beta]  
                          

        acceptance_prob = Rⱼ * Rⱼ * Rₛ * Rₛ     

        if acceptance_prob > 1 || rand(rng, Uniform(0, 1), 1)[1] < acceptance_prob
            if verbose 
                println("Hop accepted!")
                println("Rⱼ = $Rⱼ")
                println("Rₛ = $Rₛ")
                println("accept prob. = $acceptance_prob")
            end
            
            return LocalAcceptance(1, beta, beta_spin, k, l)  # acceptance, particle number, particle spin, initial site, final site
        else
            if verbose 
               println("Hop rejected!")
            end

            return LocalAcceptance(0, beta, beta_spin, k, l)
        end
    end
end


"""
    metropolis( W, jastrow1, jastrow2, particle_positions, rng )

Perform accept/reject step of proposed hop using the Metropolis algorithm. If move is
accepted, returns acceptance, particle β and it's spindex, initial position, and final
position.

"""

function metropolis(W, jastrow1, jastrow2, particle_positions, rng)    
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
        if verbose 
            println("Hop impossible!")  
        end
        return LocalAcceptance(0, beta, beta_spin, k, l)
    else
        if verbose 
            println("Hop possible!")
        end

        # begin Metropolis algorithm

        # get Jastrow ratios (element of T vector)
        Rⱼ₁ = get_jastrow_ratio(k, l, jastrow1)    
        Rⱼ₂ = get_jastrow_ratio(k, l, jastrow2)

        # get wavefunction ratio (correpsonding element of Green's function)
        Rₛ = W[l, beta]        
                           

        acceptance_prob = Rⱼ₁ * Rⱼ₁ * Rⱼ₂ * Rⱼ₂ * Rₛ * Rₛ     

        if acceptance_prob > 1 || rand(rng, Uniform(0, 1), 1)[1] < acceptance_prob
            if verbose 
                println("Hop accepted!")
                println("Rⱼ₁ = $Rⱼ₁")
                println("Rⱼ₂ = $Rⱼ₂")
                println("Rₛ = $Rₛ")
                println("accept prob. = $acceptance_prob")
            end

            return LocalAcceptance(1, beta, beta_spin, k, l)  # acceptance, particle number, particle spin, initial site, final site
        else
            if verbose 
               println("Hop rejected!")
            end

            return LocalAcceptance(0, beta, beta_spin, k, l)
        end
    end
end


"""
    metropolis( W, jastrow1, jastrow2, jastrow3, particle_positions, rng )

Perform accept/reject step of proposed hop using the Metropolis algorithm. If move is
accepted, returns acceptance, particle β and it's spindex, initial position, and final
position.

"""

function metropolis(W, jastrow1, jastrow2, jastrow3, particle_positions, rng)    
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
        if verbose
            println("Hop impossible!")  
        end
        return LocalAcceptance(0, beta, beta_spin, k, l)
    else
        if verbose
            println("Hop possible!")
        end

        # begin Metropolis algorithm

        # get Jastrow ratios (element of T vector)
        Rⱼ₁ = get_jastrow_ratio(k, l, jastrow1)    
        Rⱼ₂ = get_jastrow_ratio(k, l, jastrow2)
        Rⱼ₃ = get_jastrow_ratio(k, l, jastrow3)

        # get wavefunction ratio (correpsonding element of Green's function)
        Rₛ = W[l, beta]       
                           

        acceptance_prob = Rⱼ₁ * Rⱼ₁ * Rⱼ₂ * Rⱼ₂ * Rⱼ₃ * Rⱼ₃ * Rₛ * Rₛ     

        if acceptance_prob > 1 || rand(rng, Uniform(0, 1), 1)[1] < acceptance_prob
            if verbose
                println("Hop accepted!")
                println("Rⱼ₁ = $Rⱼ₁")
                println("Rⱼ₂ = $Rⱼ₂")
                println("Rⱼ₃ = $Rⱼ₃")
                println("Rₛ = $Rₛ")
                println("accept prob. = $acceptance_prob")
            end

            return LocalAcceptance(1, beta, beta_spin, k, l)  # acceptance, particle number, particle spin, initial site, final site
        else
            if verbose 
               println("Hop rejected!")
            end

            return LocalAcceptance(0, beta, beta_spin, k, l)
        end
    end
end



"""
    do_particle_hop!( local_acceptance::LocalAcceptance, pconfig::Matrix{Int})

If proposed particle hop is accepted, perform the particle hop.

"""
function do_particle_hop!(proposed_hop, pconfig)
    if proposed_hop.acceptance == 1
        k = proposed_hop.isite
        l = proposed_hop.fsite

        # HOP!
        pconfig[k] = 0
        pconfig[l] = 1

        return nothing
    else
        # DO NOTHING
        return nothing
    end
end


"""
    local_fermion_update!(Ne::Int, model_geometry::ModelGeometry, tight_binding_model::TightBindingModel, 
                        jastrow::Jastrow, pconfig::Vector{Int64}, rng::Xoshiro)

Perform a local MC update. Proposes moves and accept/rejects via Metropolis algorithm,
if accepted, updates particle positions, T vector, W matrix, and variational parameters.

"""
function local_fermion_update!(Np, model_geometry, tight_binding_model, jastrow, pconfig, rng)
    if verbose
        println("Starting new Monte Carlo cycle...")
    end

    # counts number of proposed hops
    proposed_hops = 0
    # counts number of accepted hops
    accepted_hops = 0

    # perform number of metropolis steps equal to the number of particles
    for s in 1:Np
        if verbose
            println("Metropolis step = $s")
        end

        # increment number of proposed hops
        proposed_hops += 1

        # get particle positions
        particle_positions = get_particle_positions(pconfig)    

        # accept/reject (Metropolis) step
        hop_step = metropolis(W, jastrow, particle_positions, rng)     

        # whether hop was accepted
        acceptance = hop_step.acceptance

        # if hop is accepted 
        if acceptance == 1
            accepted_hops += 1

            # perform hop   
            do_particle_hop!(hop_step, pconfig)                         # possible source of bug

            # update particle positions
            update_particle_position!(hop_step, particle_positions)     # possible source of bug?

            # update Green's function
            update_equal_greens!(hop_step, W)   

            # update Jastrow factors
            update_Tvec!(hop_step, jastrow, model_geometry)         

            # update variational parameters
            sr_update!(measurement_container, determinantal_parameters, jastrow, model_geometry, tight_binding_model, pconfig, Np, W, A, η, dt)
        end
    end

    # compute acceptance rate
    acceptance_rate = accepted_hops / proposed_hops     

    return acceptance_rate, pconfig, jastrow, W
end


"""
    local_fermion_update!(Ne::Int, model_geometry::ModelGeometry, tight_binding_model::TightBindingModel, 
                        jastrow1::Jastrow, jastrow2::Jastrow, pconfig::Vector{Int64}, rng::Xoshiro)

Perform a local MC update. Proposes moves and accept/rejects via Metropolis algorithm,
if accepted, updates particle positions, T vector, W matrix, and variational parameters.

"""
function local_fermion_update!(Ne, model_geometry, tight_binding_model, jastrow1, jastrow2, pconfig, rng)
    if verbose
        println("Starting new Monte Carlo cycle...")
    end
    # counts number of proposed hops
    proposed_hops = 0
    # counts number of accepted hops
    accepted_hops = 0

    # perform number of metropolis steps equal to the number of electrons
    for s in 1:Ne
        if verbose
            println("Metropolis step = $s")
        end

        # increment number of proposed hops
        proposed_hops += 1

        # get particle positions
        particle_positions = get_particle_positions(pconfig)    

        # accept/reject (Metropolis) step
        hop_step = metropolis(W, jastrow1, jastrow2, particle_positions, rng)     # TODO: test whether configuration is reverted after proposal

        # whether hop was accepted
        acceptance = hop_step.acceptance

        # if hop is accepted 
        if acceptance == 1
            accepted_hops += 1

            # perform hop
            do_particle_hop!(hop_step, pconfig)   

            # update particle positions
            update_particle_position!(hop_step, particle_positions)     

            # update Green's function
            update_equal_greens!(hop_step, W)   

            # update Jastrow factors
            update_Tvec!(hop_step, jastrow1, model_geometry)         
            update_Tvec!(hop_step, jastrow2, model_geometry)  

            # update variational parameters
            sr_update!(measurement_container, determinantal_parameters, jastrow1, jastrow2, model_geometry, tight_binding_model, pconfig, Np, W, A, η, dt)
        end
    end

    # compute acceptance rate
    acceptance_rate = accepted_hops / proposed_hops

    return acceptance_rate, pconfig, jastrow1, jastrow2, W
end


