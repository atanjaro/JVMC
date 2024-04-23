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
    if size(D, 1) != size(D, 2)
        return false
    end

    return rank(D) == size(D, 1)
end


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

# end # of module


