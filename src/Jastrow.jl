"""
    Jastrow( jastrow_type::AbstractString, jpar_matrix::Matrix{AbstractFloat}, 
             jpars::Vector{Float64}, Tvec::Vector{AbstractFloat}, jpar_map::Dict{Any,Any}, 
             num_jpars::Int )

A type defining quantities related to a Jastrow factor.

"""
mutable struct Jastrow
    # type of Jastrow parameter
    jastrow_type::String

    # matrix of Jastrow parameters
    jpar_matrix::Matrix{Float64}

    # vector of Jastrow parameters
    jpars::Vector{Float64}

    # T vector
    Tvec::Vector{Float64}

    # Jastrow parameter dictionary
    jpar_map::OrderedDict{Any, Any}

    # number of Jastrow parameters
    num_jpars::Int
end


"""
    get_distances() 

Returns a matrix of distances between each and every site 'i' and 'j', is generated with
each matrix element corresponds to a distance i.e. r₁₂ == distance between sites 1 and 2
the matrix is symmetric i.e r₁₂ = r₂₁

            [ r₁₁ r₁₂ r₁₃ r₁₄ ... ]
            [ r₂₁ r₂₂ r₂₃ r₂₄ ... ]
distances = [ r₃₁ r₃₂ r₃₃ r₃₄ ... ] 
            [ r₄₁ r₄₂ r₄₃ r₄₄ ... ]
            [ ... ... ... ... ... ]
"""
function get_distances()
    # extent of the lattice
    N = model_geometry.lattice.N

    # initialize matrix
    dist = zeros(Float64, N, N)
    
    # precompute locations
    locations = [site_to_loc(i, model_geometry.unit_cell, model_geometry.lattice)[1] for i in 1:N]

    # compute distances
    for i in 1:N
        for j in 1:N
            dist[i, j] = euclidean(locations[i], locations[j])
        end
    end

    return dist
end


# function set_jpars(model_geometry::ModelGeometry, rng::Xoshiro, readin_jpars::Bool)
#     if readin_jpars
#         # read columns of jpars of distances, indices, and parameter values
#     end

#     # extent of the lattice
#     L = model_geometry.lattice.N

#     # expected number of Jastrow parameters
#     num_jpars = div(L * (L + 1), 2)

#     # distance matrix
#     dist_matrix = get_distances()

#     # Create a dictionary of Jastrow parameters
#     jpar_map = Dict()
#     for row in 1:L
#         for col in row:L
#             # initialize all Jastrow parameters to random values
#             jpar_map[(dist_matrix[row, col], (row, col))] = rand(rng)
#         end
#     end

#     # sort dictionary
#     sorted_pairs = sort(collect(jpar_map), by=x->x[1][1])
#     jpar_map = OrderedDict(sorted_pairs)

#     # fix parameters for the largest distance to 0
#     first_elements = [key[1] for key in keys(jpar_map)]
#     unique_first_elements = Set(first_elements)
#     largest_first_element = maximum(unique_first_elements)
#     for key in keys(jpar_map)
#         if key[1] == largest_first_element
#             jpar_map[key] = 0.0
#         end
#     end

#     # verify that we have the correct number of parameters
#     length(jpar_map) == num_jpars

#     # vector of Jastrow parameters
#     jpars = collect(values(jpar_map))

#     return jpars, jpar_map, num_jpars
# end




"""
    get_jpar_matrix( model_geometry::ModelGeometry, rng::Xoshiro, readin_jpars::Bool )

Generates a matrix of distances and populates them with random initial Jastrow parameter 
values. Values can also be read-in from file.
"""
function get_jpar_matrix(model_geometry::ModelGeometry, rng::Xoshiro, readin_jpars::Bool)::Matrix{Float64}
    if readin_jpars
        # read columns of jpars of distances, indices, and parameter values
    end

    # extent of the lattice
    N = model_geometry.lattice.N

    # distance matrix
    jpar_matrix = get_distances()

    # find the largest distance and set it's parameter to 0
    max_dist = maximum(jpar_matrix)
    max_indices = findall(x -> x == max_dist, jpar_matrix)
    for index in max_indices
        jpar_matrix[Tuple(index)...] = 0.0    
    end

    for i in 1:N
        for j in i+1:N
            if jpar_matrix[i, j] != 0.0
                # initialize with random values
                value = rand(rng)
                jpar_matrix[i, j] = value
                jpar_matrix[j, i] = value
            end
        end
    end

    return jpar_matrix
end


# another version of get_jpar_matrix for recalculation during simulation
# function get_jpar_matrix(model_geometry, jpars)
    # do stuff
# end


"""
    get_jpar_map(model_geometry::ModelGeometry, jpar_matrix::Matrix{Float64})

Creates a dictionary of Jastrow parameters and their associated distances and indices (i, j),
such that the dictionary store keys and values of the form: (distance, (i, j)) => jpar_value.
Also returns vector of parameters 
"""
function get_jpar_map(model_geometry::ModelGeometry, jpar_matrix::Matrix{Float64})
    # extent of the lattice
    N = model_geometry.lattice.N

    #  # expected number of Jastrow parameters
    # num_jpars = div(N * (N + 1), 2) - N

    # distance matrix
    dist_matrix = get_distances()

    # Create dictionary of parameters
    jpar_map = Dict()

    for i in 1:N
        for j in 1:N
            if i <= j  
                jpar_map[dist_matrix[i, j], (i, j)] = jpar_matrix[i, j]
            end
        end
    end

    # Filter out entries with 0 values
    filtered_pairs = filter(x -> x[2] != 0, collect(jpar_map))

    # Sort the filtered pairs by the first element of the key tuple
    sorted_pairs = sort(filtered_pairs, by = x -> x[1][1])

    # Create the ordered dictionary from the sorted pairs
    jpar_map = OrderedDict(sorted_pairs)

    # check number of parameters
    num_jpars = length(jpar_map)

    # collect all parameters
    jpars = collect(values(jpar_map))

    return jpar_map, jpars, num_jpars
end






# """
#     set_jpars( dist_vec::Matrix{AbstractFloat}) 

# Sets entries in the distance matrix to some initial Jastrow parameter value and 
# sets parameters corresponding to the largest distance to 0.

# """
# function populate_jpars!(dist_matrix::Matrix{Float64})
#     jpar_matrix = copy(dist_matrix)
#     r_max = maximum(dist_matrix)

#     # Update jpar_matrix in place
#     for i in 1:size(dist_matrix, 1)
#         for j in 1:size(dist_matrix, 2)
#             if dist_matrix[i, j] == r_max
#                 jpar_matrix[i, j] = 0.0
#             elseif i == j
#                 jpar_matrix[i, j] = 0.0
#             else
#                 jpar_matrix[i, j] = 0.5
#             end
#         end
#     end

#     return jpar_matrix
# end

# function set_jpars(dist_matrix) 
#     jpar_matrix = copy(dist_matrix)

#     r_max = maximum(dist_matrix)
#     for i in 1:model_geometry.lattice.N
#         for j in 1:model_geometry.lattice.N
#             if jpar_matrix[i,j] == r_max
#                 jpar_matrix[i,j] = 0
#             elseif i == j
#                 jpar_matrix[i,j] == 0  
#             else
#                 jpar_matrix[i,j] = 0.5
#             end
#         end
#     end

#     return jpar_matrix
# end


# """
#     get_num_jpars( jpar_matrix::Matrix{AbstractFloat} ) 

# Returns the number of Jastrow parameters.

# """
# function get_num_jpars(jpar_matrix)
#     return count(i->(i > 0),(jpar_matrix[tril!(trues(size(jpar_matrix)), -1)]))
# end


# """
#     map_jpars( jpar_matrix::Matrix{AbstractFloat}, dist_matrix::Matrix{AbstractFloat} ) 

# Creates a dictionary of Jastrow parameters vᵢⱼ, and their distances rᵢⱼ

# """
# function map_jpars(dist_matrix, jpar_matrix)
#     upper_triangular_dist = UpperTriangular(dist_matrix)
#     upper_triangular_jpars = UpperTriangular(jpar_matrix)
#     nrows, ncols = size(dist_matrix)
    
#     jpar_dict = Dict()
    
#     for i in 1:nrows
#         for j in i:ncols
#             distance = upper_triangular_dist[i,j]
#             jpar = upper_triangular_jpars[i, j]
#             indices = (i, j)
#             jpar_dict[indices] = (distance, jpar)
#         end
#     end

#     jpar_dist_dict = Dict(key => value for (key, value) in jpar_dict if value[2] != 0.0)
    
#     return jpar_dist_dict
# end


"""
    get_Tvec( jpar_matrix::Matrix{AbstractFloat}, jastrow_type::AbstractString, pconfig::Vector{Int64} ) 

Returns vector of T with entries Tᵢ = ∑ⱼ vᵢⱼnᵢ(x) if using density Jastrow or 
Tᵢ = ∑ⱼ wᵢⱼSᵢ(x) if using spin Jastrow.

"""
function get_Tvec(jpar_matrix::Matrix{Float64}, jastrow_type::String, pconfig::Vector{Int64}, pht::Bool)
    # extent of the lattice
    N = model_geometry.lattice.N

    # initialize T vector
    Tvec = Vector{AbstractFloat}(undef, N)

    for i in 1:N
        # sum Jastrow parameters
        jpar_sum = sum(jpar_matrix[i,:])

        # check particle occupations
        num_up = number_operator(i,pconfig)[1]
        num_dn = number_operator(i,pconfig)[2]

        # check Jastrow type
        if jastrow_type == "e-den-den"  # electron density-density
            if pht
                Tvec[i] = jpar_sum * (num_up - num_dn)  
            else
                Tvec[i] = jpar_sum * (num_up + num_dn)  
            end
        elseif jastrow_type == "e-spn-spn"  # electron spin-spin
            if pht 
                Tvec[i] = 0.5 * jpar_sum * (num_up + num_dn)
            else
                Tvec[i] = 0.5 * jpar_sum * (num_up - num_dn)
            end
        elseif jastrow_type == "eph-den-den" # electron-phonon density-density
            # populate electron-phonon T vector
        elseif jastrow_type == "ph-dsp-dsp"  # phonon-displacement-displacement
            # populate
        elseif jastrow_type == "eph-den-dsp"  # electron-phonon density-displacement
            # populate
        end
    end

    return Tvec
end


"""
    update_Tvec!( local_acceptance::LocalAcceptance, jastrow::Jastrow, model_geometry::ModelGeometry )

Updates elements Tᵢ of the vector T after a Metropolis update.

"""
function update_Tvec!(local_acceptance, jastrow::Jastrow, model_geometry::ModelGeometry)
    # T vector
    Tvec = jastrow.Tvec

    # extent of the lattice
    N = model_geometry.lattice.N

    # initial and final sites from the particle hop
    l = local_acceptance.isite
    k = local_acceptance.fsite

    # update
    for i in 1:N
        Tvec[i] += jastrow.jpar_matrix[i,l] - jastrow.jpar_matrix[i,k] 
    end

    return nothing
end


"""
    get_jastrow_ratio( l::Int, k::Int, Tₗ::Vector{AbstractFloat}, Tₖ::Vector{AbstractFloat}  )

Calculates ratio J(x₂)/J(x₁) = exp[-s(Tₗ - Tₖ) + vₗₗ - vₗₖ ] of Jastrow factors for particle configurations 
which differ by a single particle hopping from site 'l' (configuration 'x₁') to site 'k' (configuration 'x₂')
using the corresponding T vectors Tₗ and Tₖ, rsepctively.  

"""
function get_jastrow_ratio(l::Int, k::Int, jastrow::Jastrow)
    # T vector
    Tvec = jastrow.Tvec

    # jpar matrix
    jpar_matrix = jastrow.jpar_matrix

    # jpar matrix
    vₗₗ =  jpar_matrix[l,l]
    vₗₖ = jpar_matrix[l,k]

    # select element for the initial and final sites
    Tₗ = Tvec[l]
    Tₖ = Tvec[k]

    # compute ratio
    jas_ratio = exp(-(Tₗ - Tₖ) + vₗₗ - vₗₖ) 

    return jas_ratio
end


"""
    build_jastrow_factor(jastrow_type::AbstractString)

Constructs relevant Jastrow factor and returns intitial T vector, matrix of Jastrow parameters, and
number of Jastrow parameters. 

"""
function build_jastrow_factor(jastrow_type::String, model_geometry::ModelGeometry, 
                                rng::Xoshiro, pconfig::Vector{Int64}, readin_jpars::Bool)
    # generate matrix of ALL Jastrow parameters
    jpar_matrix = get_jpar_matrix(model_geometry, rng, readin_jpars)

    # map Jastrow parameters
    jpar_map, jpars, num_jpars = get_jpar_map(model_geometry, jpar_matrix)

    # generate T vector
    init_Tvec = get_Tvec(jpar_matrix, jastrow_type, pconfig, pht)
   
    if verbose
        # report the number of Jastrow parameters initialized
        println(num_jpars," Jastrow parameters initialized")
        println("Type: ", jastrow_type)
    end

    return Jastrow(jastrow_type, jpar_matrix, jpars, init_Tvec, jpar_map, num_jpars)
end


"""
    update_jastrow!(jastrow::Jastrow, vpars::Vector{AbstractFloat})

Updates Jastrow parameters after Stochastic Reconfiguration.

"""
function update_jastrow!(jastrow::Jastrow, vpars::Vector{Float64})
    # number of Jastrow parameters
    num_jpars = jastrow.num_jpars;

    # number of Jastrow parameters
    jpars = jastrow.jpars;

    # map of Jastrow parameters
    jpar_map = jastrow.jpar_map;

    # number of Jastrow parameters
    jpar_matrix = jastrow.jpar_matrix;

    # get the new Jastrow parameters by getting the last num_jpars elements of vpars
    new_jpars = vpars[end-num_jpars+1:end]

    # overwrite current jpars 
    jpars .= new_jpars

    # overwrite jpars in map
    list_of_keys = collect(keys(jpar_map))
    for (i, key) in enumerate(list_of_keys)
        jpar_map[key] = new_jpars[i]
    end

    # push new elements to matrix
    for (key, value) in jpar_map
        # Extract the indices from the key
        _, (i, j) = key
        # Update the matrix at the specified indices
        jpar_matrix[i, j] = value
        jpar_matrix[j, i] = value
    end

    return nothing
end


"""
    recalc_Tvec(Tᵤ::Vector{AbstractFloat}, δT::AbstractFloat)

Checks floating point error accumulation in the T vector and if ΔT < δT, 
then the recalculated T vector Tᵣ replaces the updated T vector Tᵤ.

"""
function recalc_Tvec!(jastrow::Jastrow, δT::Float64)
    # Jastrow type
    jastrow_type = jastrow.jastrow_type

    # T vector that has been updated during MC cycles
    Tᵤ = jastrow.Tvec

    # get matrix
    jpar_matrix = jastrow.jpar_matrix
    # jpar_matrix = get_jpar_matrix(model_geometry, jpars)

    # recomputed T vector
    Tᵣ = get_Tvec(jpar_matrix, jastrow_type, pconfig, pht)

    # difference in updated T vector and recalculated T vector
    diff = Tᵤ .- Tᵣ

    # sum the absolute differences and the recalculated T vector elements
    diff_sum = sum(abs.(diff))
    T_sum = sum(abs.(Tᵣ))

    # rms difference
    ΔT = sqrt(diff_sum / T_sum)

    if ΔT > δT
        verbose && println("WARNING! T vector has been recalculated: ΔT = ", ΔT, " > δT = ", δT)

        # record new T vector
        jastrow.Tvec = Tᵣ

        return nothing
    else
        verbose && println("T vector is stable: ΔT = ", ΔT, " < δT = ", δT)
        return nothing
    end  
end


