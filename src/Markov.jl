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
    # create neighbor table
    nbr_table = build_neighbor_table(bonds[1],
                                    model_geometry.unit_cell,
                                    model_geometry.lattice)

    # maps neighbor table to dictionary of bonds and neighbors                                
    nbrs = map_neighbor_table(nbr_table)

    # randomly select some particle in the lattice
    beta = rand(rng, 1:trunc(Int,Np))                   

    # real position 'k' of particle 'β' 
    k = particle_positions[beta][2]     
    
    # randomly selected neighboring site 'l'
    k_nbrs = nbrs[k][2]
    nbr_rand = rand(rng, 1:length(k_nbrs))
    l = k_nbrs[nbr_rand]          

    # spin of particle particle 'β' 
    beta_spin = get_spindex_type(particle_positions[beta][1],model_geometry)
    
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
        particle_positions = get_particle_positions(pconfig, model_geometry)    

        # accept/reject (Metropolis) step
        met_step = metropolis(W, jastrow, particle_positions, rng)    

        # whether hop was accepted
        acceptance = met_step.acceptance

        # DEBUG
        if debug
            prop_particle = met_step.particle
            prop_spin = met_step.spin
            prop_isite = met_step.isite
            prop_fsite = met_step.fsite

            println("Particle: $prop_particle")
            println("Spin: $prop_spin")
            println("initial site: $prop_isite")
            println("final site: $prop_fsite")

            @info "Before update:"
            @info "particle_positions: $particle_positions"
            @info "pconfig: $pconfig"
        end

        # if hop is accepted 
        if acceptance == 1
            accepted_hops += 1

            # perform hop   
            do_particle_hop!(met_step, pconfig, model_geometry)                 

            # update particle positions
            update_particle_position!(met_step, particle_positions)     

            # # update particle positions 
            # particle_positions = get_particle_positions(pconfig)  

            # update Green's function
            update_equal_greens!(met_step, W)   

            # update Jastrow factors
            update_Tvec!(met_step, jastrow, model_geometry)         

            # update variational parameters
            sr_update!(measurement_container, determinantal_parameters, jastrow, model_geometry, tight_binding_model, pconfig, particle_positions, Np, W, A, η, dt)
        end
        if debug
            @info "After update:"
            @info "particle_positions: $particle_positions"
            @info "pconfig: $pconfig"

            println("Length of particle_positions: ", length(particle_positions))
        end
    end

    # compute acceptance rate
    local_acceptance_rate = accepted_hops / proposed_hops     

    return local_acceptance_rate, pconfig, jastrow, W
end





