# module StochasticReconfiguration

using LinearAlgebra

"""
    local_jastrow_derivative(jastrow::Jastrow, pconfig::Vector{AbstractFloat}, jastrow_type::AbstractString)

Performs local logarithmic derivative Δₖ(x) = ∂lnΨ(x)/∂vₗₘ, with respect
to the kth Jastrow parameter vₗₘ.

"""
function local_jastrow_derivative(jastrow, pconfig, model_geometry)
    if jastrow.jastrow_type == "density"
        result = zeros(AbstractFloat, jastrow.num_jpars)
        for jpar in 1:jastrow.num_jpars
            for i in 1:model_geometry.lattice.N
                for j in 1:model_geometry.lattice.N
                    if pht == true
                        result[jpar] = -(number_operator(i,pconfig)[1] - number_operator(i,pconfig)[2])*(number_operator(j,pconfig)[1] - number_operator(j,pconfig)[2])
                    else
                        result[jpar] = -(number_operator(i,pconfig)[1] + number_operator(i,pconfig)[2])*(number_operator(j,pconfig)[1] + number_operator(j,pconfig)[2])
                    end
                end
            end
        end
    elseif jastrow_type == "spin"
        result = zeros(AbstractFloat, jastrow.num_jpars)
        for jpar in 1:jastrow.num_jpars
            for i in 1:model_geometry.lattice.N
                for j in 1:model_geometry.lattice.N
                    if pht == true
                        result[jpar] = -0.5 * (number_operator(i,pconfig)[1] + number_operator(i,pconfig)[2])*(number_operator(j,pconfig)[1] + number_operator(j,pconfig)[2])
                    else
                        result[jpar] = -0.5 * (number_operator(i,pconfig)[1] - number_operator(i,pconfig)[2])*(number_operator(k,pconfig)[1] - number_operator(j,pconfig)[2])
                    end
                end
            end
        end
    elseif jastrow_type == "electron-phonon"
        result = zeros(AbstractFloat, jastrow.num_jpars)
        # TBA
    else
    end
    return result
end



"""
    local_determinantal_derivative( A::Matrix{AbstractFloat}, W::Matrix{AbstractFloat}, acceptance::LocalAcceptance, 
                                    model_geometry::ModelGeometry, parameters_to_optimize::Vector{AbstractString})

Performs local logarithmic derivative Δₖ(x) = ∂lnΨ(x)/∂αₖ, with respect
to the kth variational parameter αₖ, in the determinantal part of the wavefunction.

"""
function local_determinantal_derivative(A, W, acceptance, model_geometry, parameters_to_optimize)
    dims = model_geometry.unit_cell.n * model_geometry.lattice.N
    num_vpars = length(parameters_to_optimize)
    result = zeros(AbstractFloat, num_vpars)
    iᵦ = acceptance.isite

    # j = acceptance.fsite
    
    G = zeros(AbstractFloat, 2*dims, 2*dims)

    for β in 1:Np
        G[iᵦ,:] = W[:,β]
    end

    for vpar in 1:num_vpars
        result[vpar] = sum(A[vpar] * G)
    end

    # # sum over number of particles
    # for β in 1:Np
    #     # sum over lattice sites for the up and down sectors
    #     for j in 1:2*model_geometry.lattice.N
    #         Δₖ += Ak[iᵦ,j]*W[j,β]
    #     end
    # end

    return result
end


"""
    get_SR_comatrix()

Generates the covariance matrix S, for Stochastic Reconfiguration

"""
function get_SR_comatrix()
    # Δk = local_determinantal_derivative(A, W, acceptance, model_geometry, parameters_to_optimize)
    # Δk_Δkprime = local_determinantal_derivative(A, W, acceptance, model_geometry, parameters_to_optimize) * local_determinantal_derivative(A, W, acceptance, model_geometry, parameters_to_optimize)

    S = Δk_Δkprime - Δk * transpose(Δk)
    
    return S
end


"""
    get_SR_forces()

Generates the force vector f, for Stochastic Reconfiguration.

"""
function get_SR_forces()
    Δk = local_determinantal_derivative(A, W, acceptance, model_geometry, parameters_to_optimize)
    E = measure_local_energy(model_geometry, tight_binding_model, jastrow, particle_positions)
    Δk_E = Δk * E

    f = Δk * E - Δk_E

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
    parameter_update!()

Update variational parameters.

"""
function parameter_update!()
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