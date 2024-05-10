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
    # matrix of Jastrow parameters
    jpar_matrix::Matrix{AbstractFloat}
    # T vector
    Tvec::Vector{AbstractFloat}
    # Jastrow parameter dictionary
    jpar_map::Dict{Any,Any}
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
    set_jpars( dist_vec::Matrix{AbstractFloat}) 

Sets entries in the distance matrix to some initial Jastrow parameter value and 
sets parameters corresponding to the largest distance to 0.

"""
function set_jpars(dist_matrix) 
    jpar_matrix = copy(dist_matrix)

    r_max = maximum(dist_matrix)
    for i in 1:model_geometry.lattice.N
        for j in 1:model_geometry.lattice.N
            if jpar_matrix[i,j] == r_max
                jpar_matrix[i,j] = 0
            elseif i == j
                jpar_matrix[i,j] == 0  
            else
                jpar_matrix[i,j] = 0.5
            end
        end
    end

    return jpar_matrix
end


"""
    get_num_jpars( jpar_matrix::Matrix{AbstractFloat} ) 

Returns the number of Jastrow parameters.

"""
function get_num_jpars(jpar_matrix)
    return count(i->(i > 0),(jpar_matrix[tril!(trues(size(jpar_matrix)), -1)]))
end


"""
    map_jpars( jpar_matrix::Matrix{AbstractFloat}, dist_matrix::Matrix{AbstractFloat} ) 

Creates a dictionary of Jastrow parameters vᵢⱼ, and their distances rᵢⱼ

"""
function map_jpars(dist_matrix, jpar_matrix)
    upper_triangular_dist = UpperTriangular(dist_matrix)
    upper_triangular_jpars = UpperTriangular(jpar_matrix)
    nrows, ncols = size(dist_matrix)
    
    jpar_dict = Dict()
    
    for i in 1:nrows
        for j in i:ncols
            distance = upper_triangular_dist[i,j]
            jpar = upper_triangular_jpars[i, j]
            indices = (i, j)
            jpar_dict[indices] = (distance, jpar)
        end
    end

    jpar_dist_dict = Dict(key => value for (key, value) in jpar_dict if value[2] != 0.0)
    
    return jpar_dist_dict
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
function update_Tvec!(local_acceptance, jastrow, model_geometry)
    N = model_geometry.lattice.N

    l = local_acceptance.isite
    k = local_acceptance.fsite

    for i in 1:N
        jastrow.Tvec[i] += jastrow.jpar_matrix[i,l] - jastrow.jpar_matrix[i,k] 
    end

    return nothing
end


"""
    get_jastrow_ratio( l::Int, k::Int, Tₗ::Vector{AbstractFloat}, Tₖ::Vector{AbstractFloat}  )

Calculates ratio J(x₂)/J(x₁) = exp[-s(Tₗ - Tₖ) + vₗₗ - vₗₖ ] of Jastrow factors for particle configurations 
which differ by a single particle hopping from site 'l' (configuration 'x₁') to site 'k' (configuration 'x₂')
using the corresponding T vectors Tₗ and Tₖ, rsepctively.  

"""
function get_jastrow_ratio(l, k, jastrow)
    Tₗ = jastrow.Tvec[l]
    Tₖ = jastrow.Tvec[k]

    jas_ratio = exp(-(Tₗ - Tₖ) + jastrow.jpar_matrix[l,l] - jastrow.jpar_matrix[l,k]) #jastrow.jpar_map[(l,k)]

    return jas_ratio
end


"""
    build_jastrow_factor(jastrow_type::AbstractString)

Constructs relevant Jastrow factor and returns intitial T vector, matrix of Jastrow parameters, and
number of Jastrow parameters. 

"""
function build_jastrow_factor(jastrow_type)
    dist_matrix = get_distances()
    jpar_matrix = set_jpars(dist_matrix)
    num_jpars = get_num_jpars(jpar_matrix)
    jpar_map = map_jpars(dist_matrix, jpar_matrix)

    # report the number of Jastrow parameters
    if verbose == true
        println(num_jpars," Jastrow parameters initialized")
        println("Type: ", jastrow_type)
    end
    init_Tvec = get_Tvec(jpar_matrix,jastrow_type)

    return Jastrow(jastrow_type, jpar_matrix, init_Tvec, jpar_map, num_jpars)
end


"""
    recalc_Tvec(Tᵤ::Vector{AbstractFloat}, δT::AbstractFloat)

Checks floating point error accumulation in the T vector and if ΔT < δT, 
then the recalculated T vector Tᵣ replaces the updated T vector Tᵤ.

"""
#TODO: need to update method
function recalc_Tvec(jastrow, δT, model_geometry)
    N = model_geometry.lattice.N
    Tᵣ = get_Tvec(jpar_matrix, jastrow_type)
    diff = Tᵤ-Tᵣ
    diff_sum = 0.0
    T_sum = 0.0

    for i in 1:2*N
        diff_sum += sum(diff[i])
        T_sum += sum(Tᵣ[i])
    end

    ΔT = sqrt(diff_sum^2/T_sum)     

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





