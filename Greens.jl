# module Greens

using LinearAlgebra

# export get_equal_greens
# export update_equal_gree
# export recalc_equal_greens

"""
    get_equal_greens(M::Matrix{AbstractFloat}, D::Matrix{AbstractFloat}) 
    
Returns the equal-time Green's function (overlap ratios) matrix W by solving 
DᵀWᵀ = Mᵀ using full pivot LU decomposition.

"""
function get_equal_greens(M, D)
    # transpose M and D
    Dt = transpose(D)
    Mt = transpose(M)

    # perform LU decomposition of Dᵀ
    Ft = lu(Dt, Val(true))

    # define W matrix
    W = zeros(Np, 2*model_geometry.lattice.N)

    # solve the equation
    Wt = Dt \ Mt      

    # Update the entries of the W matrix
    W = transpose(Wt)
 
    return W                

    # return M * inv(D)       # Testing purposes only
end                         


"""
    update_equal_greens!( local_acceptance::LocalAcceptance ) 
    
Update equal-time Green's function. 

"""
function update_equal_greens!(local_acceptance)
    # get rₗ, the lth row of W
    rₗ = W[local_acceptance.fsite,:]

    # get cᵦ, the βth column vetor of W
    cᵦ = W[:,local_acceptance.particle]

    # subtract 1 from βth component of rₗ
    rₗ[local_acceptance.particle] =- 1

    # get W'
    W_new = W - (cᵦ * rₗ / W[local_acceptance.fsite,local_acecptance.particle])

    return W_new
end


"""
    recalc_equal_greens( Wᵤ::Matrix{AbstractFloat}, δW::AbstractFloat ) 
    
Checks floating point error accumulation in the equal-time Green's function
and if ΔW < δW, then the recalculated Green's function
Wᵣ replaces the updated Green's function Wᵤ.

"""
function recalc_equal_greens(Wᵤ, δW)
    Wᵣ = get_equal_greens(M, D)
    # ΔW = sqrt(sum( (Wᵤ-Wᵣ)^2 )/sum(Wᵣ))     
    if ΔW > δW
        if verbose == true
            println("WARNING! Green's function has been recalculated: ΔW = ", ΔW, " > δW = ", δW)
        end
        return Wᵣ, ΔW

    else # ΔW < δW
        if verbose == true
            println("Green's function is stable: ΔW = ", ΔW, " > δW = ", δW)
        end
        return Wᵤ, ΔW
    end  
end


# end # of module