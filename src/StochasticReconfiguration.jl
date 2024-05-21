# module StochasticReconfiguration

using LinearAlgebra

# # for the covariance matrix...
# # initial sum for a single element f
# ΔkΔp = 0.0

# # stuff happens here...

# ΔkΔp += derivative


# # for the force vector...
# # initial sum for a single element of f
# ΔkE = 0.0

# # stuff happens here...

# ΔkE += derivative


"""
    get_sr_comatrix()

Generates the covariance matrix S, for Stochastic Reconfiguration

The matrix S has elements S_kk' = <Δ_kΔk'> - <Δ_k><Δ_k'>

"""
function get_sr_comatrix()

    # total number of variational parameters
    num_detpars = determinantal_parameters.num_detpars
    num_jpars = density_jastrow.num_jpars 
    num_vpars = num_detpars + num_jpars 

    # initialize covariance matrix
    S = zeros(AbstractFloat, num_vpars, num_vpars)

    # measure local parameters derivatives for this configuration
    Δk = measure_Δk(determinantal_parameters, jastrow, model_geometry, pconfig, Np, W, A)

    for (i,j) in zip(num_vpars, num_vpars)
        # populate SR matrix
    end

   

    # S = Δk_Δkprime - Δk * transpose(Δk)
    
    return S
end


"""
    get_sr_forces( determinantal_parameters::DeterminantalParameters)

Generates the force vector f, for Stochastic Reconfiguration.

The vector f has elements f_k = <Δ_k><H> - <Δ_kH>

"""
function get_sr_forces(determinantal_parameters, jastrow, model_geometry, tight_binding_model, pconfig, Np, W, A )

    # particle positions
    particle_positions = get_particle_positions(pconfig)

    # total number of variational parameters
    num_detpars = determinantal_parameters.num_detpars
    num_jpars = jastrow.num_jpars 
    num_vpars = num_detpars + num_jpars 
    
    # initialize force vector
    f = zeros(AbstractFloat, num_vpars)

    # measure local parameters derivatives for this configuration
    j_derivatives = measure_local_jpar_derivative(jastrow, pconfig)
    d_derivatives = measure_local_detpar_derivative(determinantal_parameters, model_geometry, pconfig, Np, W, A)

    # concatenate derivatives
    derivatives = vcat(d_derivatives,j_derivatives)

    # measure local energy for this configuration
    (E_loc, E_loc_kinetic, E_loc_hubb) = measure_local_energy(model_geometry, tight_binding_model, jastrow, particle_positions)

    # for i in num_vpars


    return f
end


"""
    parameter_gradient( S::Matrix{AbstractFloat}, η::Float64 )

Perform gradient descent on variational parameters for Stochastic 
Reconfiguration.

"""
function parameter_gradient(S, η)
    # add small variation to diagonal of S
    S += η + I

    # solve for δα using LU decomposition
    δvpar = S \ f
    
    return δvpar
end


"""
    sr_update!()

Update variational parameters.

"""
function sr_update!()
    # get covariance matrix
    S = get_SR_matrix()
    # get force vector
    f = get_SR_forces()

    # perform gradient descent
    δvpar = parameter_gradient(S,η)     
    # update parameter
    vpar += dt * δvpar      # start with a large dt and reduce as energy is minimized

    return vpar
end



"""
    parameter_indices()

Get indices of variational parameters from its respective matrix.

"""
function get_parameter_indices(par_matrix)
    nonzero_indices = findall(x -> x != 0, par_matrix)

    parameter_indices = sort(nonzero_indices, by=x->(x[1], x[2]))

    return parameter_indices
end 