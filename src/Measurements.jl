using LatticeUtilities
using LinearAlgebra


"""
initialize_measurement_container(model_geometry::ModelGeometry)

Creates generic arrays for storage of associated measurment quantities. 

"""
function initialize_measurement_container(model_geometry, N_burnin, N_updates, variational_parameters)
    unit_cell = model_geometry.unit_cell
    lattice = model_geometry.lattice

    L = lattice.L
    N = lattice.N
    norbs = unit_cell.n
    N_iterations = N_burnin + N_updates

    local_measurements = 0.0 

    scalar_measurements = Dict{String, Any}([
        ("density", zeros(AbstractFloat, norbs)),        # average density per orbital species
        ("double_occ", zeros(AbstractFloat, norbs))      # average double occupancy per orbital species
    ])                                                           

    derivative_measurements = Dict{String, Any}([
        ("Δk", (zeros(AbstractFloat, norbs), local_measurements, Any[])),        
        ("ΔkΔkp", (zeros(AbstractFloat, norbs), local_measurements, Any[])),
        ("ΔkE", (zeros(AbstractFloat, norbs), local_measurements, Any[]))                  
    ])      

    parameter_measurements = Dict{String, Any}([
        ("parameters", Any[])    
    ])
    # initialize parameters values
    parameter_measurements["parameters"] = [variational_parameters]  

    correlation_measurements = Dict()

    measurement_container = (
        scalar_measurements       = scalar_measurements,
        derivative_measurements   = derivative_measurements,     
        parameter_measurements    = parameter_measurements,     
        correlation_measurements  = correlation_measurements,                       
        L                         = L,
        N                         = N,
        norbs                     = norbs,
        N_iterations              = N_iterations                           
    )

    return measurement_container
end



"""
initialize_measurements!(measurement_container::, observable::AbstractString )

For a certain type of measurment (scalar or correlation), initializes the arrays
necessary to store measurements in respective bins.

"""
function initialize_measurements!(measurement_container, observable)
    (; scalar_measurements,N,norbs) = measurement_container

    local_measurements = 0.0 
    expectation = []

    if observable == "energy"
        # For energy, store tuple (zeros(AbstractFloat, norbs), local_measurements, expectation)
        scalar_measurements["energy"] = (zeros(AbstractFloat, norbs),local_measurements,expectation)
        # scalar_measurements["global_energy"] = local_measurements, expectation
    elseif observable == "stripe"
        # For stripe, store tuple (zeros(AbstractFloat, norbs*N), local_measurements, expectation)
        scalar_measurements["site_dependent_density"] = (zeros(AbstractFloat, norbs*N),local_measurements,expectation)
        scalar_measurements["site_dependent_spin"] = (zeros(AbstractFloat, norbs*N),local_measurements,expectation)
    elseif observable == "phonon_number"
        # For phonon number
    elseif observable == "displacement"
        # For displacement
    end

    return nothing
end


# TODO: implement correlation measurmenets
function initialize_correlation_measurements!(measurement_container,  correlation)
    # type of measurements
    #   - density-density correlation 
    #   - spin-spin correlation
    #   - pair correlation

    (; correlation_measurements) = measurement_container

    if correlation == "den-den"
        # scalar_measurements["density-density correlation"] = (zeros(AbstractFloat, norbs),local_measurements = 0.0,expectation = [])
    elseif correlation == "spn-spn"
        # scalar_measurements["spin-spin correlation"] = (zeros(AbstractFloat, norbs),local_measurements = 0.0,expectation = [])
    elseif correlation == "pair"
        # scalar_measurements["pair correlation"] = (zeros(AbstractFloat, norbs),local_measurements = 0.0,expectation = [])
    end

    return nothing
end

# accumulator for measurements
# this is reflected in the key entries for each measurement dictionary
    # # through the course of the simulation, we estimate the expectation value of observable O using ⟨O⟩ ≈ N⁻¹ ∑ₓ Oₗ(x)
    # by 'local' here, I mean the argument of the sums used to obtain the expectation value
    # local_measurements = 0.0         # Oₗ(x). This is added to during the course of the simulation: local_val += measured
    # iterations = 0.0        # N. This is incremented during the simulation, up to N_burnin + N_updates: iterations +=1
    # by 'global' here, I mean the sum of local measurements to obtain the expectation value
    # expectation = []   # record the expectation value for each iteration: push!(total_local_vals, local_val/iterations)


"""
    get_local_jpar_derivative( jastrow::Jastrow, pconfig::Vector{Int} )

Calculates the local logarithmic derivative Δₖ(x) = ∂lnΨ(x)/∂vₗₘ, with respect to the kth Jastrow parameter vₗₘ. Returns 
a vector of derivatives.

"""
function get_local_jpar_derivative(jastrow, pconfig)
    # jastrow type
    jastrow_type = jastrow.jastrow_type;

    # number of Jastrow parameters
    num_jpars = jastrow.num_jpars;

    # map of Jastrow parametersß
    jpar_map = jastrow.jpar_map;

    # vector to store derivatives
    derivatives = zeros(AbstractFloat, num_jpars)
                
    # for density Jastrow
    if jastrow_type == "e-den-den"
        for num in 1:num_jpars
            for (r1i, jpars1) in jpar_map         #former: indices => (dist, jpar) : (1, 2) => (1.0, jpar)
                for (r2j, jpars2) in jpar_map
                    if r1i[1] == r2j[1]
                        if pht == true
                            derivatives[num] += -(number_operator(r1i[2][1],pconfig)[1] - number_operator(r1i[2][1],pconfig)[2]) * (
                                                  number_operator(r1i[2][2],pconfig)[1] - number_operator(r1i[2][2],pconfig)[2])
                        elseif pht == false
                            derivatives[num] += -(number_operator(r1i[2][1],pconfig)[1] + number_operator(r1i[2][1],pconfig)[2]) * (
                                                  number_operator(r1i[2][2],pconfig)[1] + number_operator(r1i[2][2],pconfig)[2])
                        else
                        end
                    else
                    end
                end
            end
        end

        return derivatives

    # for spin Jastrow
    elseif jastrow_type == "e-spin-spin"   
        for num in 1:num_jpars
            for (r1i, jpars1) in jpar_map         
                for (r2j, jpars2) in jpar_map
                    if r1i[1] == r2j[1]
                        if pht == true
                            derivatives[num] += -(number_operator(r1i[2][1],pconfig)[1] - number_operator(r1i[2][1], pconfig)[2]) * (
                                                  number_operator(r1i[2][2], pconfig)[1] - number_operator(r1i[2][2], pconfig)[2])
                        elseif pht == false
                            derivatives[num] += -(number_operator(r1i[2][1],pconfig)[1] + number_operator(r1i[2][1], pconfig)[2]) * (
                                                  number_operator(r1i[2][2], pconfig)[1] + number_operator(r1i[2][2], pconfig)[2])
                        else
                        end
                    end
                end
            end
        end

        return derivatives

    # for electron-phonon Jastrow
    elseif jastrow_type == "eph-den-den"   
        return derivatives
    else
    end
end


"""
    get_local_detpar_derivative( determinantal_parameters::DeterminantalParameters, model_geometry::ModelGeometry
                                     pconfig::Vector{Int}, W::Matrix{AbstractFloat}, A::Matrix{AbstractFloat}  )

Calculates the local logarithmic derivative Δₖ(x) = ∂lnΨ(x)/∂αₖ, with respect to the kth variational parameter αₖ,
in the determinantal part of the wavefunction. Returns a vector of derivatives.

"""
function get_local_detpar_derivative(determinantal_parameters, model_geometry, particle_positions, Np, W, A)  

    # dimensions
    dims = model_geometry.unit_cell.n * model_geometry.lattice.N

    # number of determinantal parameters
    num_detpars = determinantal_parameters.num_detpars
    
    # particle positions
    # particle_positions = get_particle_positions(pconfig)

    # vector to store derivatives
    derivatives = zeros(AbstractFloat, num_detpars)
    

    # loop over Nₚ particles 
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
    measure_Δk!(measurement_container, determinantal_parameters, jastrow, model_geometry, pconfig, Np, W, A)

Measures logarithmic derivatives for all variational parameters. The first 'p' are derivatives of 
determinantal parameters and the rest are derivatives of Jastrow parameters. Measurments are then written
to the measurement container.

"""
# PASSED
function measure_Δk!(measurement_container, determinantal_parameters, jastrow, model_geometry, pconfig, particle_positions, Np, W, A)
    # perform derivatives
    detpar_derivatives = get_local_detpar_derivative(determinantal_parameters, model_geometry, particle_positions, Np, W, A)
    jpar_derivatives = get_local_jpar_derivative(jastrow, pconfig)
    Δk = vcat(detpar_derivatives,jpar_derivatives)

    # record current expectation values
    local_measurement = measurement_container.derivative_measurements["Δk"][2] .+ Δk   
    current_expectation = local_measurement / measurement_container.N_iterations

    # write to measurement container
    push!(measurement_container.derivative_measurements["Δk"][3], current_expectation)
    updated_container = (measurement_container.derivative_measurements["Δk"][1], local_measurement, measurement_container.derivative_measurements["Δk"][3])
    measurement_container.derivative_measurements["Δk"] = updated_container

    return nothing
end


"""
    measure_ΔkE( determinantal_parameters, jastrow, model_geometry, tight_binding_model, pconfig, Np, W, A )

Measures the product of variational derivatives with the local energy. Measurments are then written
to the measurement container.

"""
# PASSED
function measure_ΔkE!(measurement_container, determinantal_parameters, jastrow, model_geometry, tight_binding_model, pconfig, particle_positions,  Np, W, A)
    # perform derivatives
    detpar_derivatives = get_local_detpar_derivative(determinantal_parameters, model_geometry, particle_positions, Np, W, A)
    jpar_derivatives = get_local_jpar_derivative(jastrow,pconfig)
    Δk = vcat(detpar_derivatives,jpar_derivatives)

    # compute local energy
    E = get_local_energy(model_geometry, tight_binding_model,jastrow,pconfig, particle_positions) 

    # compute product of local derivatives with the local energy
    ΔkE = Δk * E

    # record current expectation values
    local_measurement = measurement_container.derivative_measurements["ΔkE"][2] .+ ΔkE
    current_expectation = local_measurement / measurement_container.N_iterations

    # write to measurement container
    push!(measurement_container.derivative_measurements["ΔkE"][3], current_expectation)
    updated_container = (measurement_container.derivative_measurements["ΔkE"][1], local_measurement, measurement_container.derivative_measurements["ΔkE"][3])
    measurement_container.derivative_measurements["ΔkE"] = updated_container

    return nothing
end


"""
    measure_ΔkΔkp(  )

Measures the product of variational derivatives with other variational derivatives. Measurments are then written
to the measurement container.

"""
# PASSED
function measure_ΔkΔkp!(measurement_container, determinantal_parameters, jastrow, model_geometry, pconfig, particle_positions, Np, W, A)
    # perform derivatives
    detpar_derivatives = get_local_detpar_derivative(determinantal_parameters, model_geometry, particle_positions, Np, W, A)
    jpar_derivatives = get_local_jpar_derivative(jastrow,pconfig)
    Δk = vcat(detpar_derivatives,jpar_derivatives)

    # inner product of Δk and Δk′
    ΔkΔkp = Δk * Δk'

    # record current expectation values
    local_measurement = measurement_container.derivative_measurements["Δk"][2] .+ ΔkΔkp
    current_expectation = local_measurement / measurement_container.N_iterations

    # write to measurement container
    push!(measurement_container.derivative_measurements["ΔkΔkp"][3], current_expectation)
    updated_container = (measurement_container.derivative_measurements["ΔkΔkp"][1], local_measurement, measurement_container.derivative_measurements["ΔkΔkp"][3])
    measurement_container.derivative_measurements["ΔkΔkp"] = updated_container

    return nothing
end 


"""
    get_local_energy(model_geometry::ModelGeometry, tight_binding_model::TightBindingModel, jastrow::Jastrow, particle_positions:: )

Calculates the local variational energy. Returns the total local energy and writes to the measurement container.

"""
# PASSED
function get_local_energy(model_geometry, tight_binding_model, jastrow, pconfig, particle_positions)
    # number of sites
    N = model_geometry.lattice.N

    # generate neighbor table
    nbr_table = build_neighbor_table(bonds[1],
                                    model_geometry.unit_cell,
                                    model_geometry.lattice)

    # gnerate neighbor map
    nbr_map = map_neighbor_table(nbr_table)

    # particle positions
    # particle_positions = get_particle_positions(pconfig)

    E_loc_kinetic = 0.0
    E_loc_hubbard = 0.0

    # calculate electron kinetic energy
    for β in 1:Np
        # if β > length(particle_positions)
        #     error("Index β ($β) is out of bounds for particle_positions of length $(length(particle_positions))")
        # end
        # if length(particle_positions[β]) < 2
        #     error("particle_positions[β] does not have at least 2 elements: ", particle_positions[β])
        # end
        k = particle_positions[β][2]
        # println("β: ", β, ", k: ", k)

        # loop over nearest neighbors. TBA: loop over different neighbor orders (i.e. nearest and next nearest neighbors)
        sum_nn = 0.0
        for l in nbr_map[k][2]
            # reverse sign if system is particle-hole transformed
            if pht == true
                Rⱼ = exp(-get_jastrow_ratio(l, k, jastrow))
            else
                Rⱼ = exp(get_jastrow_ratio(l, k, jastrow))
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
    for i in 1:N
        E_loc_hubbard += U * number_operator(i, pconfig)[1] * (1 - number_operator(i, pconfig)[2])
    end
    

    # calculate total local energy
    E_loc = E_loc_kinetic + E_loc_hubbard

    return E_loc
end


"""
    get_local_energy(model_geometry::ModelGeometry, tight_binding_model::TightBindingModel, 
                    jastrow1::Jastrow, jastrow2::Jastrow, particle_positions:: )

Calculates the local variational energy. Returns the total local energy and writes to the measurement container.

"""
function get_local_energy(model_geometry, tight_binding_model, jastrow1, jastrow2, pconfig, particle_positions)
    # number of sites
    N = model_geometry.lattice.N

    # generate neighbor table
    nbr_table = build_neighbor_table(bonds[1],
                                    model_geometry.unit_cell,
                                    model_geometry.lattice)

    # gnerate neighbor map
    nbr_map = map_neighbor_table(nbr_table)

    # particle positions
    # particle_positions = get_particle_positions(pconfig)

    E_loc_kinetic = 0.0
    E_loc_hubbard = 0.0

    # calculate electron kinetic energy
    for β in 1:Np
        # loop over different electrons k
        k = particle_positions[β][2] 
        # loop over nearest neighbors. TBA: loop over different neighbor orders (i.e. nearest and next nearest neighbors)
        sum_nn = 0.0
        for l in nbr_map[k][2]
            # reverse sign if system is particle-hole transformed
            if pht == true
                Rⱼ₁ = exp(-get_jastrow_ratio(l, k, jastrow1))
                Rⱼ₂ = exp(-get_jastrow_ratio(l, k, jastrow2))
            else
                Rⱼ₁ = exp(get_jastrow_ratio(l, k, jastrow1))
                Rⱼ₂ = exp(get_jastrow_ratio(l, k, jastrow2))
            end
            sum_nn += Rⱼ₁ * Rⱼ₂ * W[l, β]
        end

        # reverse sign if system is particle-hole transformed
        if pht == true
            E_loc_kinetic += tight_binding_model.t[1] * sum_nn          
        else
            E_loc_kinetic += - tight_binding_model.t[1] * sum_nn
        end
    end

    # calculate Hubbard energy
    for i in 1:N
        E_loc_hubbard += U * number_operator(i, pconfig)[1] * (1 - number_operator(i, pconfig)[2])
    end
    

    # calculate total local energy
    E_loc = E_loc_kinetic + E_loc_hubbard

    return E_loc
end


"""
    get_local_energy(model_geometry::ModelGeometry, tight_binding_model::TightBindingModel, 
                    jastrow1::Jastrow, jastrow2::Jastrow, particle_positions:: )

Calculates the local variational energy. Returns the total local energy and writes to the measurement container.

"""
function get_local_energy(model_geometry, tight_binding_model, jastrow1, jastrow2, jastrow3, pconfig, particle_positions)
    # number of sites
    N = model_geometry.lattice.N

    # generate neighbor table
    nbr_table = build_neighbor_table(bonds[1],
                                    model_geometry.unit_cell,
                                    model_geometry.lattice)

    # gnerate neighbor map
    nbr_map = map_neighbor_table(nbr_table)

    # particle positions
    # particle_positions = get_particle_positions(pconfig)

    E_loc_kinetic = 0.0
    E_loc_hubbard = 0.0

    # calculate electron kinetic energy
    for β in 1:Np
        # loop over different electrons k
        k = particle_positions[β][2] 
        # loop over nearest neighbors. TBA: loop over different neighbor orders (i.e. nearest and next nearest neighbors)
        sum_nn = 0.0
        for l in nbr_map[k][2]
            # reverse sign if system is particle-hole transformed
            if pht == true
                Rⱼ₁ = exp(-get_jastrow_ratio(l, k, jastrow1))
                Rⱼ₂ = exp(-get_jastrow_ratio(l, k, jastrow2))
                Rⱼ₃ = exp(-get_jastrow_ratio(l, k, jastrow3))
            else
                Rⱼ₁ = exp(get_jastrow_ratio(l, k, jastrow1))
                Rⱼ₂ = exp(get_jastrow_ratio(l, k, jastrow2))
                Rⱼ₃ = exp(get_jastrow_ratio(l, k, jastrow3))
            end
            sum_nn += Rⱼ₁ * Rⱼ₂ * Rⱼ₃ * W[l, β]
        end

        # reverse sign if system is particle-hole transformed
        if pht == true
            E_loc_kinetic += tight_binding_model.t[1] * sum_nn          
        else
            E_loc_kinetic += - tight_binding_model.t[1] * sum_nn
        end
    end

    # calculate Hubbard energy
    for i in 1:N
        E_loc_hubbard += U * number_operator(i, pconfig)[1] * (1 - number_operator(i, pconfig)[2])
    end
    

    # calculate total local energy
    E_loc = E_loc_kinetic + E_loc_hubbard

    return E_loc
end



"""
    measure_local_energy!( measurement_container, model_geometry::ModelGeometry, tight_binding_model::TightBindingModel, 
                jastrow::Jastrow, particle_positions:: )

Measures the total local energy and writes to the measurement container.

"""
# PASSED
function measure_local_energy!(measurement_container, model_geometry, tight_binding_model, jastrow, pconfig, particle_positions)
    # calculate the current local energy
    E_loc = get_local_energy(model_geometry, tight_binding_model, jastrow, pconfig, particle_positions)

    # record current expectation values
    local_measurement = measurement_container.scalar_measurements["energy"][2] .+ E_loc
    current_expectation = local_measurement / measurement_container.N_iterations

    # write to measurement container
    push!(measurement_container.scalar_measurements["energy"][3], current_expectation)
    updated_container = (measurement_container.scalar_measurements["energy"][1], local_measurement, measurement_container.scalar_measurements["energy"][3])
    measurement_container.scalar_measurements["energy"] = updated_container

    return nothing
end



"""
    measure_global_energy( model_geometry::ModelGeometry )

Measure the global variational energy ⟨E⟩. This is intended to used at the end of the simulation.

"""
function measure_global_energy(model_geometry, N_bins, bin_size)

    # binned energies
    # E_binned

    # average over all sites
    E_global = E_binned / model_geometry.lattice.N

    return E_global
end


"""
    measure_double_occ( model_geometry::ModelGeometry, pconfig::Vector{Int} )

Measure the average double occupancy ⟨D⟩ = N⁻¹ ∑ᵢ ⟨nᵢ↑nᵢ↓⟩.

"""
function measure_double_occ!(measurement_container, pconfig, model_geometry)
    nup_ndn = 0.0

    for i in 1:model_geometry.lattice.N
        nup_ndn += number_operator(i, pconfig)[1] * (1 - number_operator(i, pconfig)[2])
    end

    dblocc = nup_ndn / model_geometry.lattice.N

    # write to measurment container
    
    
    return nothing
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


