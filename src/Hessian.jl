"""
    get_hessian_matrix(measurement_container, determinantal_parameters, 
                    jastrow, model_geometry, tight_binding_model, pconfig, 
                    Np, W, A )

Generates the covariance (Hessian) matrix S, for Stochastic Reconfiguration

The matrix S has elements S_kk' = <Δ_kΔk'> - <Δ_k><Δ_k'>

"""
function get_hessian_matrix(measurement_container, determinantal_parameters, jastrow, model_geometry, pconfig, particle_positions,Np, W, A )

    # measure local parameters derivatives ⟨Δₖ⟩ (also ⟨Δₖ'⟩), for this configuration
    measure_Δk!(measurement_container, determinantal_parameters, jastrow, model_geometry, pconfig, particle_positions, Np, W, A)
    Δk = measurement_container.derivative_measurements["Δk"][3][end]        

    # measure products of local derivatives ⟨ΔₖΔₖ'⟩, for this configuration 
    measure_ΔkΔkp!(measurement_container, determinantal_parameters, jastrow, model_geometry, pconfig, particle_positions, Np, W, A)
    ΔkΔkp = measurement_container.derivative_measurements["ΔkΔkp"][3][end]   
    
    # calculate the product of local derivatives ⟨Δk⟩⟨Δkp⟩
    ΔkΔk = Δk * Δk'  

    # generate covariance matrix
    S = ΔkΔkp - ΔkΔk
    
    return S
end


"""
    get_force_vector(measurement_container, determinantal_parameters, 
                jastrow, model_geometry, tight_binding_model, pconfig,
                 Np, W, A )

Generates the force vector f, for Stochastic Reconfiguration.

The vector f has elements f_k = <Δ_k><H> - <Δ_kH>

"""
# PASSED
function get_force_vector(measurement_container, determinantal_parameters, jastrow, model_geometry, tight_binding_model, pconfig, particle_positions, Np, W, A )
    
    # initialize force vector
    f = [] 

    # measure local parameters derivatives ⟨Δₖ⟩, for this configuration
    measure_Δk!(measurement_container, determinantal_parameters, jastrow, model_geometry, pconfig, particle_positions,Np, W, A)
    Δk = measurement_container.derivative_measurements["Δk"][3][end]         

    # measure local energy E = ⟨H⟩, for this configuration
    measure_local_energy!(measurement_container, model_geometry, tight_binding_model, jastrow, pconfig, particle_positions)
    E = measurement_container.scalar_measurements["energy"][3][end]      

    # measure product of local derivatives with energy ⟨ΔkE⟩, for this configuration
    measure_ΔkE!(measurement_container, determinantal_parameters, jastrow, model_geometry, tight_binding_model, pconfig, particle_positions, Np, W, A)
    ΔkE = measurement_container.derivative_measurements["ΔkE"][3][end]       

    # product of local derivative with the local energy ⟨Δk⟩⟨H⟩
    ΔktE = Δk * E
    
    for (i,j) in zip(ΔktE,ΔkE)
        fk = i - j
        push!(f, fk)
    end

    return f  # the length of f == number of vpars where the first p are the determinantal parameters and the rest are Jastrow parameters
end


"""
    parameter_gradient( S::Matrix{AbstractFloat}, η::Float64 )

Perform gradient descent on variational parameters for Stochastic 
Reconfiguration.

"""
function parameter_gradient(S, f, η)

    # Convert f to a vector of Float64
    f = convert(Vector{Float64}, f)

    # add small variation to diagonal of S for numerical stabilization
    S += η * I

    # solve for δα using LU decomposition
    δvpars = S \ f
    
    return δvpars   # the length of δvpars == number of vpars where the first p are the determinantal parameters and the rest are Jastrow parameters
end


"""
    sr_update!(measurement_container, determinantal_parameters, 
                jastrow, model_geometry, tight_binding_model, 
                pconfig, Np, W, A, η, dt)

Update variational parameters through stochastic optimization.

"""
function sr_update!(measurement_container, determinantal_parameters, jastrow, model_geometry, tight_binding_model, pconfig, particle_positions,Np, W, A, η, dt)
    # get covariance (Hessian) matrix
    S = get_hessian_matrix(measurement_container, determinantal_parameters, jastrow, model_geometry, pconfig, particle_positions, Np, W, A )
    # get force vector
    f = get_force_vector(measurement_container, determinantal_parameters, jastrow, model_geometry, tight_binding_model, pconfig, particle_positions, Np, W, A)

    # perform gradient descent
    δvpars = parameter_gradient(S,f,η)     

    # update ALL parameters
    vpars = cat_vpars(determinantal_parameters, jastrow)
    vpars += dt * δvpars    # TODO: start with a large dt and reduce as energy is minimized
    push!(measurement_container.parameter_measurements["parameters"], vpars)

    # push back Jastrow parameters
    update_jastrow!(jastrow, vpars)

    return nothing
end

