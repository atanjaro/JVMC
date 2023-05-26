using LatticeUtilities
using LinearAlgebra
using DelimitedFiles
using Distances

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
    set_jpars( dist_vec::Matrix{AbstractFloat}, vpar_init::AbstractFloat ) 

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
                dist_matrix[i,j] = jpar_init
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
    get_Tvec( jpar_vec::Vector{AbstractFloat} ) 

Returns vector of T with entries Tᵢ = ∑ⱼ vᵢⱼnᵢ(x).

"""
function get_Tvec(jpar_matrix)
    Tvec = Vector{AbstractFloat}(undef, model_geometry.lattice.N)
    # if den_jastrow == true
        for i in 1:model_geometry.lattice.N
            Tvec[i] = sum(jpar_matrix[i,:]) * ( number_operator(i,pconfig)[1] + number_operator(i,pconfig)[2] )  
        end
    # elseif spn_jastrow == true
    #     # TODO: include spin Jastrow terms
    # end

    return Tvec
end


"""
    update_Tvec!( tvec::Vector{AbstractFloat} )

Updates elements Tᵢ of the vector T after a Metropolis update.

"""
function update_Tvec!(l::Int, k::Int, tvec)
    for i in 1:L
        Tvec[i] += (jpar_matrix[i,l] - jpar_matrix[i,k]) #TODO: include sgn for PHT
    end

    return Tvec
end


"""
    get_jastrow_ratio( l::Int, k::Int, tvec::Vector{AbstractFloat} )

Calculates ratio J(x₂)/J(x₁) of Jastrow factors for particle configurations which
differ by a single particle hopping from site 'l' (configuration 'x₁') to site 'k' (configuration 'x₂'). 

"""
function get_jastrow_ratio(l, k, tvec)
    # a particle hops from site l to site k
    # J(x₂)/J(x₁) = exp[-(Tₗ - Tₖ) - vₗₗ - vₗₖ ]
    return nothing
end


"""
    get_jastrow_factor()

Constructs relevant Jastrow factors.

"""
function get_jastrow_factor()
    jpar_matrix = get_distances()
    set_jpars!(jpar_matrix)
    num_jpars = get_num_jpars(jpar_matrix)
    # report the number of Jastrow parameters
    if verbose == true
        println(num_jpars," Jastrow parameters initialized")
    end
    Tvec = get_Tvec(jpar_matrix)

    return Tvec, jpar_matrix, num_jpars
end





