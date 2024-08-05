# module Utilities

using LinearAlgebra

# export diagonalize
# export is_invertible
# export complex_zeros
# export is_openshell

"""
    diagonalize( H::Matrix{AbstractFloat} ) 

Returns all eigenenergies and all eigenstates of the mean-field Hamiltonian, 
the latter being stored in the columns of a matrix. Convention: H(diag) = U⁺HU.

"""
function diagonalize(H)
    # check if Hamiltonian is Hermitian
    @assert ishermitian(H) == true
    F = eigen(H_mf)    
    return F.values, F.vectors  
end


"""
    is_invertible( D::Matrix{AbstractFloat} ) 
    
Checks if given matrix is invertible by checking its rank.

"""
function is_invertible(D)
    try
        F = lu(D)
        # Check if there are any zeroes on the diagonal of the U matrix
        return all(abs.(diag(F.U)) .> eps(eltype(D)))
    catch e
        if isa(e, SingularException)
            return false
        else
            rethrow(e)
        end
    end
end

# function is_invertible_cond(A::AbstractMatrix; tol::Real = 1e12)
#     # tol is a threshold for the condition number
#     return cond(A) < tol
# end

# A = [1.0 2.0; 3.0 4.0]
# B = [1.0 2.0; 2.0 4.0]
# is_invertible_cond(A)
# is_invertible_cond(B)

# function is_invertible_lu(A::AbstractMatrix)
#     try
#         F = lu(A)
#         # Check if there are any zeroes on the diagonal of the U matrix
#         return all(abs.(diag(F.U)) .> eps(eltype(A)))
#     catch e
#         if isa(e, SingularException)
#             return false
#         else
#             rethrow(e)
#         end
#     end
# end

# # Example usage
# is_invertible_lu(A)
# is_invertible_lu(B)



"""
    complex_zeros{T<:AbstractFloat}(rows::Int, cols::Int)

Initialzies a zero matrix with complex type

"""
function complex_zeros( rows::Int, cols::Int)
    return zeros(Complex{AbstractFloat}, rows, cols)
end


"""
    is_openshell( Np::Int, ε::Vector{AbstractFloat} ) 

Checks whether configuration is open shell.

"""
function is_openshell(ε, Np)
    return ε[Np + 1] - ε[Np] < 0.0001
end


"""
    get_num_vpars(determinantal_parameters::DeterminantalParameters, 
    density_jastrow::Jastrow, spin_jastrow::Jastrow, eph_jastrow::Jastrow) 

Returns the total number of variational parameters.

"""
function get_num_vpars()
    num_jpars = 0

    num_detpars = determinantal_parameters.num_detpars

    # check if density Jastrow is defined
    if isdefined(Main, :density_jastrow)
        num_jpars += density_jastrow.num_jpars
    end

    # check if spin Jastrow is defined
    if isdefined(Main, :spin_jastrow)
        num_jpars += spin_jastrow.num_jpars
    end

    # check if e-ph Jastrow is defined
    if isdefined(Main, :eph_jastrow)
        num_jpars += eph_jastrow.num_jpars
    end

    num_vpars = num_detpars + num_jpars

    return num_vpars, num_detpars, num_jpars

end


"""
    cat_vpars(determinantal_parameters, jastrow)

Generates an array combining all values of determinantal and Jastrow parameters.

"""
function cat_vpars(determinantal_parameters, jastrow)
    jpars = jastrow.jpars
    detpars = determinantal_parameters.vals
    
    # concatenate variational parameters
    vpar_array = vcat(detpars, jpars)
    
    return vpar_array
end


"""
    get_parameter_indices()

Get indices of variational parameters from its respective matrix.

"""
function get_parameter_indices(par_matrix)
    nonzero_indices = findall(x -> x != 0, par_matrix)

    parameter_indices = sort(nonzero_indices, by=x->(x[1], x[2]))

    return parameter_indices
end 

# end # of module


