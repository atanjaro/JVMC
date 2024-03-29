# module Jastrow

using LatticeUtilities
using LinearAlgebra
using DelimitedFiles
using Distances

# export build_jastrow_factor
# export get_Tvec
# export update_Tvec

struct Jastrow
    # type of Jastrow parameter
    jastrow_type::AbstractString
    # T vector
    Tvec::Vector{AbstractFloat}
    # Jastrow parameter matrix
    jpar_matrix::Matrix{AbstractFloat}
    # number of Jastrow parameters
    num_jpars::Int
end

# a matrix of distances between each and every site is generated with
# each matrix element corresponds to a distance i.e. r₁₂ == distance between sites 1 and 2
# the matrix is symmetric i.e r₁₂ = r₂₁
#
#             [ r₁₁ r₁₂ r₁₃ r₁₄ ... ]
#             [ r₂₁ r₂₂ r₂₃ r₂₄ ... ]
# distances = [ r₃₁ r₃₂ r₃₃ r₃₄ ... ] 
#             [ r₄₁ r₄₂ r₄₃ r₄₄ ... ]
#             [ ... ... ... ... ... ]
#

"""
    get_distances() 

Returns a matrix of all possible distances between sites 'i' nd 'j' for a given lattice.

"""
function get_distances()
    dist = zeros(AbstractFloat, model_geometry.lattice.N, model_geometry.lattice.N)
    for i in 1:model_geometry.lattice.N
        for j in 1:model_geometry.lattice.N
            dist[i,j] = euclidean(site_to_loc(i,model_geometry.unit_cell,model_geometry.lattice)[1],
                                site_to_loc(j,model_geometry.unit_cell,model_geometry.lattice)[1])
        end
    end
 
    return dist
end


"""
    set_jpars!( dist_vec::Matrix{AbstractFloat}) 

Sets entries in the distance matrix to some initial Jastrow parameter value and 
sets parameters corresponding to the largest distance to 0.

"""
function set_jpars!(dist_matrix) 
    r_max = maximum(dist_matrix)
    for i in 1:model_geometry.lattice.N
        for j in 1:model_geometry.lattice.N
            if dist_matrix[i,j] == r_max
                dist_matrix[i,j] = 0
            elseif i == j
                dist_matrix[i,j] == 0  
            else
                dist_matrix[i,j] = 0.5
            end
        end
    end

    return dist_matrix
end


"""
    get_num_jpars( jpar_matrix::Matrix{AbstractFloat} ) 

Returns the number of Jastrow parameters.

"""
function get_num_jpars(jpar_matrix)
    return count(i->(i > 0),(jpar_matrix[tril!(trues(size(jpar_matrix)), -1)]))
end


"""
    get_Tvec( jpar_matrix::Matrix{AbstractFloat}, jastrow_type::AbstractString ) 

Returns vector of T with entries Tᵢ = ∑ⱼ vᵢⱼnᵢ(x) if using density Jastrow or 
Tᵢ = ∑ⱼ wᵢⱼSᵢ(x) if using spin Jastrow.

"""
function get_Tvec(jpar_matrix, jastrow_type)
    Tvec = Vector{AbstractFloat}(undef, model_geometry.lattice.N)
    for i in 1:model_geometry.lattice.N
        if jastrow_type == "density"
            if pht == true
                Tvec[i] = sum(jpar_matrix[i,:]) * (number_operator(i,pconfig)[1] - number_operator(i,pconfig)[2])  
            else
                Tvec[i] = sum(jpar_matrix[i,:]) * (number_operator(i,pconfig)[1] + number_operator(i,pconfig)[2])  
            end
        elseif jastrow_type == "spin"
            if pht == true
                Tvec[i] = sum(jpar_matrix[i,:]) * 0.5 * (number_operator(i,pconfig)[1] + number_operator(i,pconfig)[2])
            else
                Tvec[i] = sum(jpar_matrix[i,:]) * 0.5 * (number_operator(i,pconfig)[1] - number_operator(i,pconfig)[2])
            end
        elseif jastrow_type == "electron-phonon"
            # populate electron-phonon T vector
        end
    end

    return Tvec
end


"""
    update_Tvec!( tvec::Vector{AbstractFloat} )

Updates elements Tᵢ of the vector T after a Metropolis update.

"""
function update_Tvec!(l, k, Tvec)
    for i in 1:L
        Tvec[i] += (jpar_matrix[i,l] - jpar_matrix[i,k]) 
    end

    return Tvec
end


"""
    get_jastrow_ratio( l::Int, k::Int, Tₗ::Vector{AbstractFloat}, Tₖ::Vector{AbstractFloat}  )

Calculates ratio J(x₂)/J(x₁) = exp[-s(Tₗ - Tₖ) + vₗₗ - vₗₖ ] of Jastrow factors for particle configurations 
which differ by a single particle hopping from site 'l' (configuration 'x₁') to site 'k' (configuration 'x₂')
using the corresponding T vectors Tₗ and Tₖ, rsepctively.  

"""
function get_jastrow_ratio(l, k, Tₗ, Tₖ)
    jas_ratio = exp(-(Tₗ - Tₖ) + jpar_matrix[l,l] - jpar_matrix[l,k])

    return jas_ratio
end


"""
    build_jastrow_factor(jastrow_type::AbstractString)

Constructs relevant Jastrow factor and returns intitial T vector, matrix of Jastrow parameters, and
number of Jastrow parameters. 

"""
function build_jastrow_factor(jastrow_type)
    jpar_matrix = get_distances()
    set_jpars!(jpar_matrix)
    num_jpars = get_num_jpars(jpar_matrix)
    # report the number of Jastrow parameters
    if verbose == true
        println(num_jpars," Jastrow parameters initialized")
        println("Type: ", jastrow_type)
    end
    init_Tvec = get_Tvec(jpar_matrix,jastrow_type)

    return Jastrow(jastrow_type, init_Tvec, jpar_matrix, num_jpars)
end


"""
    recalc_Tvec(Tᵤ::Vector{AbstractFloat}, δT::AbstractFloat)

Checks floating point error accumulation in the T vector and if ΔT < δT, 
then the recalculated T vector Tᵣ replaces the updated T vector Tᵤ.

"""
function recalc_Tvec(Tᵤ, δT)
    Tᵣ = get_Tvec(jpar_matrix, jastrow_type)
    # ΔT = sqrt(sum( (Tᵤ-Tᵣ)^2 )/sum(Tᵣ))     
    if ΔT > δT
        if verbose == true
            println("WARNING! T vector has been recalculated: ΔT = ", ΔT, " > δT = ", δT)
        end
        return Tᵣ, ΔT

    else # ΔT < δT
        if verbose == true
            println("T vector is stable: ΔT = ", ΔT, " > δT = ", δT)
        end
        return Tᵤ, ΔT
    end  
end


# end # of module





