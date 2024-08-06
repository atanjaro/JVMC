"""
    build_determinantal_state() 

Returns initial energies ε₀, matrix M, and Slater matrix D in the many-particle configuration 
basis.

"""
function build_determinantal_state()
    # Diagonalize Hamiltonian
    ε, U = diagonalize(H_mf)

    # Check for open shell configuration
    if is_openshell(ε, Np)
        verbose && println("WARNING! Open shell detected")
    else
        verbose && println("Generating shell...")
    end

    # Store energies
    ε₀ = ε[1:Np]

    # Store M matrix
    M = hcat(U[:,1:Np])

    # Build Slater determinant
    while true
        pconfig = generate_initial_fermion_configuration()
        config_indices = findall(x -> x == 1, pconfig)
        D = M[config_indices, :]

        # Check that starting configuration is not singular
        if is_invertible(D)
            # Write matrices to file if needed
            if write
                writedlm("H_mf.csv", H_mf)
                writedlm("D.csv", D)
                writedlm("M.csv", M)
                writedlm("U.csv", U)
            end

            return D, pconfig, ε, ε₀, M, U
        end
    end    
end


"""
    get_equal_greens( M::Matrix{Float64}, D::Matrix{Float64} )::Matrix{Float64}
    
Returns the equal-time Green's function by solving DᵀWᵀ = Mᵀ using full pivot LU decomposition.

"""
function get_equal_greens(M::Matrix{Float64}, D::Matrix{Float64})::Matrix{Float64}
    verbose && println("Getting equal-time Green's function...")

    # perform the linear solve directly
    Wt = D' \ M'     

    # transpose the result back to the original shape
    W = Wt'

    debug && println("W = $W")
 
    return W                
end          


"""
    update_equal_greens!( local_acceptance::LocalAcceptance, W::Matrix{Float64} ) 
    
Perform in-place update of the equal-time Green's function. 

"""
function update_equal_greens!(local_acceptance::LocalAcceptance, W::Matrix{Float64})
    # Get the indices
    fsite = local_acceptance.fsite
    particle = local_acceptance.particle
    
    # Get rₗ, the lth row of W
    rₗ = view(W, fsite, :)
    
    # Get cᵦ, the βth column vector of W
    cᵦ = view(W, :, particle)
    
    # Subtract 1 from the βth component of rₗ
    rₗ[particle] -= 1

    # Update W in place
    W .-= cᵦ * rₗ' / W[fsite, particle]

    return nothing
end


"""
    recalc_equal_greens( Wᵤ::Matrix{Float64}, δW::Float64, D::Matrix{Float64}, pconfig::Vector{Int64} ) 
    
Checks floating point error accumulation in the equal-time Green's function and if ΔW < δW, then the 
recalculated Green's function Wᵣ replaces the updated Green's function Wᵤ as well as the D matrix
for the current configuration.

"""
function recalc_equal_greens(Wᵤ::Matrix{Float64}, δW::Float64, D::Matrix{Float64}, pconfig::Vector{Int64})
    # L = model_geometry.lattice.N
    # Np = size(Wᵤ, 2)  # Assuming Np is the number of columns in Wᵤ

    # recalculate D for current configuration
    Dᵣ = M[findall(x -> x == 1, pconfig), :]
    
    # Recalculate Green's function from scratch
    Wᵣ = get_equal_greens(M, Dᵣ)
    
    # Difference in updated Green's function and recalculated Green's function
    diff = Wᵤ .- Wᵣ

    # Sum the absolute differences and the recalculated Green's function elements
    diff_sum = sum(abs.(diff))
    W_sum = sum(abs.(Wᵣ))

    ΔW = sqrt(diff_sum / W_sum)

    if ΔW > δW
        verbose && println("WARNING! Green's function has been recalculated: ΔW = ", ΔW, " > δW = ", δW)
        return Wᵣ, Dᵣ
    else
        verbose && println("Green's function is stable: ΔW = ", ΔW, " < δW = ", δW)
        return Wᵤ, D
    end  
end


