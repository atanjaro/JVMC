# module Hamiltonian

using LatticeUtilities
using LinearAlgebra
using Test

# export ModelGeometry
# export TightBindingModel
# export initialize_determinantal_parameters
# export build_mean_field_hamiltonian
# export build_slater_determinant
# export get_Ak_matrices


"""
    ModelGeometry( unit_cell::UnitCell, lattice::Lattice )

A type defining model geometry

"""
struct ModelGeometry
    # unit cell
    unit_cell::UnitCell
    # extent of the lattice
    lattice::Lattice
    # lattice bonds
    bond::Vector{Vector{Bond{2}}}
end


"""
    TightBindingModel( t::Vector{AbstractFloat}, μ::AbstractFloat, 
                    model_geometry::ModelGeometry, nbr_table::Matrix{Int64} )

A type defining a non-interacting tight binding model

"""
struct TightBindingModel    # TODO: add onsite energy?
    # hopping amplitudes
    t::Vector{AbstractFloat}   # [t, t']    # TODO: change to [[t₁, t₂, t₃], t'] for SSH models
    # chemical potentials
    μ::AbstractFloat           # TODO: change to vector based on number of orbitals [μ₁,μ₂,μ₃,...]
end


"""
    DeterminantalParameters( pars::Vector{AbstractString}, 
                        vals::Vector{AbstractFloat}, num_detpars::Int )

A type defining a set of variational parameters.

"""
struct DeterminantalParameters
    # name of order parameter
    pars::Vector{AbstractString}
    # variational parameter values
    vals::Vector{AbstractFloat}
    # number of determinantal parameters
    num_detpars::Int
end


"""
    initialize_determinantal_parameters(pars:;Vector{AbstractString}, vals::Vector{AbstractFloat} ) 

Constructor for the variational parameters type.

"""
function initialize_determinantal_parameters(pars, vals)
    @assert length(pars) == length(vals) "Input vectors must have the same length"
    num_detpars = length(pars)

    return DeterminantalParameters(pars, vals, num_detpars)
end


"""
    map_determinantal_parameters( determinantal_parameters::DeterminantalParameters ) 

For a given set of variational parameters, returns a dictionary of 
that reports the value and optimization flag for a given parameter.

"""
function map_determinantal_parameters(determinantal_parameters)
    vparam_map = Dict()
    for i in 1:length(determinantal_parameters.vals)
       vparam_map[determinantal_parameters.pars[i]] = determinantal_parameters.vals[i]
    end
    return vparam_map
end


"""
    build_tight_binding_model( tight_binding_model::TightBindingModel ) 

Constructs a 2⋅n⋅N by 2⋅n⋅N Hamiltonian matrix, where n is the
number of orbitals per unit cell and N is the number of lattice sites,
given tight binding parameters t, t', and μ. TODO: change tight binding
parameters to tx,ty,t' for future SSH model functionality.

"""
function build_tight_binding_model(tight_binding_model, model_geometry, parameters_to_optimize, pht, bonds, Lx, Ly)
    dims = model_geometry.unit_cell.n * model_geometry.lattice.N
    nbr_table = build_neighbor_table(bonds[1], model_geometry.unit_cell, model_geometry.lattice)
    H_t = zeros(Float64, 2 * dims, 2 * dims)
    H_tp = zeros(Float64, 2 * dims, 2 * dims)
    μ_vec = Vector{Float64}(undef, 2 * dims)
    
    fill_μ_vec!(μ_vec, tight_binding_model.μ, dims, pht)
    
    fill_H_t!(H_t, nbr_table, tight_binding_model.t[1], model_geometry, Lx, Ly, pht)
    
    if tight_binding_model.t[2] != 0.0
        nbr_table_p = build_neighbor_table(bonds[2], model_geometry.unit_cell, model_geometry.lattice)
        fill_H_tp!(H_tp, nbr_table_p, tight_binding_model.t[2], model_geometry, Lx, Ly, pht)
    end

    if "μ_BCS" in parameters_to_optimize
        return H_t + H_tp
    else
        return H_t + H_tp + LinearAlgebra.Diagonal(μ_vec)
    end
end

function fill_μ_vec!(μ_vec, μ, dims, pht)
    for i in 1:dims
        μ_vec[i] = -μ
        μ_vec[dims + i] = pht ? μ : -μ
    end
end

function fill_H_t!(H_t, nbr_table, t1, model_geometry, Lx, Ly, pht)
    factor = pht ? -1 : 1
    for (i, j) in eachcol(nbr_table)
        H_t[i, j] += factor * t1
        if model_geometry.lattice.N > 2
            H_t[j, i] += factor * t1
        end
    end
    for (i, j) in eachcol(nbr_table .+ model_geometry.lattice.N)
        H_t[i, j] += -factor * t1
        if model_geometry.lattice.N > 2
            H_t[j, i] += -factor * t1
        end
    end
end

function fill_H_tp!(H_tp, nbr_table_p, t2, model_geometry, Lx, Ly, pht)
    factor = pht ? 1/2 : 1/2  # Halve the factor for 2x2 lattice to avoid double counting
    for (i, j) in eachcol(nbr_table_p)
        H_tp[i, j] += factor * t2
        H_tp[j, i] += factor * t2
    end
    for (i, j) in eachcol(nbr_table_p .+ model_geometry.lattice.N)
        H_tp[i, j] += -factor * t2
        H_tp[j, i] += -factor * t2
    end
end

# function build_tight_binding_model(tight_binding_model)
#     dims = model_geometry.unit_cell.n*model_geometry.lattice.N 
#     nbr_table = build_neighbor_table(bonds[1],
#                                     model_geometry.unit_cell,
#                                     model_geometry.lattice)
#     H_t = zeros(AbstractFloat, 2*dims, 2*dims)
#     H_tp = zeros(AbstractFloat, 2*dims, 2*dims)
#     μ_vec = Vector{AbstractFloat}(undef, 2*dims)
#     if pht == true
#         # particle-hole transformed chemical potential
#         if !("μ_BCS" in parameters_to_optimize)     
#             for i in 1:dims
#                 for j in dims+1:2*dims
#                     μ_vec[i] = -tight_binding_model.μ
#                     μ_vec[j] = tight_binding_model.μ
#                 end 
#             end
#         end
#         # particle-hole transformed nearest neighbor hopping
#         if Lx == 2 && Ly == 2 
#             for (i,j) in eachcol(nbr_table)
#                 H_t[i,j] += -tight_binding_model.t[1]
#             end
#             for (i,j) in eachcol(nbr_table .+ model_geometry.lattice.N)    
#                 H_t[i,j] += tight_binding_model.t[1]
#             end
#         # special case for 1D
#         elseif  Lx == 1 && Ly > Lx || Ly == 1 && Lx > Ly
#             for (i,j) in eachcol(nbr_table[:,1:model_geometry.lattice.N])
#                 H_t[i,j] += -tight_binding_model.t[1]
#                 if model_geometry.lattice.N > 2
#                     H_t[j,i] += -tight_binding_model.t[1]
#                 end
#             end
#             for (i,j) in eachcol(nbr_table[:,1:model_geometry.lattice.N] .+ model_geometry.lattice.N)    
#                 H_t[i,j] += tight_binding_model.t[1]
#                 if model_geometry.lattice.N > 2
#                     H_t[j,i] += tight_binding_model.t[1]
#                 end
#             end
#         # special case for Lx = 2 
#         elseif Lx == 2 && Ly > Lx
#             for (i,j) in eachcol(nbr_table[:,1:(size(nbr_table,2) - Ly)])
#                 H_t[i,j] += -tight_binding_model.t[1]
#                 H_t[j,i] += -tight_binding_model.t[1]
#             end
#             for (i,j) in eachcol(nbr_table[:,1:(size(nbr_table,2) - Ly)] .+ model_geometry.lattice.N)
#                 H_t[i,j] += tight_binding_model.t[1]
#                 H_t[j,i] += tight_binding_model.t[1]
#             end 
#         # special case for Ly = 2 
#         elseif Ly == 2 && Lx > Ly
#             for (i,j) in eachcol(nbr_table[:,1:(size(nbr_table,2) - Lx)])
#                 H_t[i,j] += -tight_binding_model.t[1]
#                 H_t[j,i] += -tight_binding_model.t[1]
#             end
#             for (i,j) in eachcol(nbr_table[:,1:(size(nbr_table,2) - Lx)] .+ model_geometry.lattice.N)
#                 H_t[i,j] += tight_binding_model.t[1]
#                 H_t[j,i] += tight_binding_model.t[1]
#             end 
#         else
#             for (i,j) in eachcol(nbr_table)
#                 H_t[i,j] += -tight_binding_model.t[1]
#                 if model_geometry.lattice.N > 2
#                     H_t[j,i] += -tight_binding_model.t[1]
#                 else
#                 end
#             end
#             for (i,j) in eachcol(nbr_table .+ model_geometry.lattice.N)    
#                 H_t[i,j] += tight_binding_model.t[1]
#                 if model_geometry.lattice.N > 2
#                     H_t[j,i] += tight_binding_model.t[1]
#                 else
#                 end
#             end
#         end
#         # particle-hole transformed next nearest neighbor hopping
#         if tight_binding_model.t[2] != 0.0
#             nbr_table_p = build_neighbor_table(bonds[2],
#                                             model_geometry.unit_cell,
#                                             model_geometry.lattice)
#             if Lx == 2 && Ly == 2
#                 for (i,j) in eachcol(nbr_table_p)
#                     H_tp[i,j] += tight_binding_model.t[2]/2
#                 end
#                 for (i,j) in eachcol(nbr_table_p .+ model_geometry.lattice.N)    
#                     H_tp[i,j] += -tight_binding_model.t[2]/2
#                 end
#             else
#                 for (i,j) in eachcol(nbr_table_p)
#                     H_tp[i,j] += tight_binding_model.t[2]
#                     H_tp[j,i] += tight_binding_model.t[2]
#                 end
#                 for (i,j) in eachcol(nbr_table_p .+ model_geometry.lattice.N)    
#                     H_tp[i,j] += -tight_binding_model.t[2]
#                     H_tp[j,i] += -tight_binding_model.t[2]
#                 end
#             end
#         else
#         end
#     else
#         # chemical potential
#         if !("μ_BCS" in parameters_to_optimize)
#             for i in 1:dims
#                 for j in dims+1:2*dims
#                     μ_vec[i] = -tight_binding_model.μ
#                     μ_vec[j] = -tight_binding_model.μ
#                 end 
#             end
#         end
#         # nearest neighbor hopping
#         if Lx == 2 && Ly == 2 
#             for (i,j) in eachcol(nbr_table)
#                 H_t[i,j] += -tight_binding_model.t[1]
#             end
#             for (i,j) in eachcol(nbr_table .+ model_geometry.lattice.N)    
#                 H_t[i,j] += -tight_binding_model.t[1]
#             end
#         # special case for 1D  
#         elseif  Lx == 1 && Ly > Lx || Ly == 1 && Lx > Ly
#             for (i,j) in eachcol(nbr_table[:,1:model_geometry.lattice.N])
#                 H_t[i,j] += -tight_binding_model.t[1]
#                 if model_geometry.lattice.N > 2
#                     H_t[j,i] += -tight_binding_model.t[1]
#                 end
#             end
#             for (i,j) in eachcol(nbr_table[:,1:model_geometry.lattice.N] .+ model_geometry.lattice.N)    
#                 H_t[i,j] += -tight_binding_model.t[1]
#                 if model_geometry.lattice.N > 2
#                     H_t[j,i] += -tight_binding_model.t[1]
#                 end
#             end
#         # special case for Lx = 2 
#         elseif Lx == 2 && Ly > Lx
#             for (i,j) in eachcol(nbr_table[:,1:(size(nbr_table,2) - Ly)])
#                 H_t[i,j] += -tight_binding_model.t[1]
#                 H_t[j,i] += -tight_binding_model.t[1]
#             end
#             for (i,j) in eachcol(nbr_table[:,1:(size(nbr_table,2) - Ly)] .+ model_geometry.lattice.N)
#                 H_t[i,j] += -tight_binding_model.t[1]
#                 H_t[j,i] += -tight_binding_model.t[1]
#             end 
#         # special case for Ly = 2 
#         elseif Ly == 2 && Lx > Ly
#             for (i,j) in eachcol(nbr_table[:,1:(size(nbr_table,2) - Lx)])
#                 H_t[i,j] += -tight_binding_model.t[1]
#                 H_t[j,i] += -tight_binding_model.t[1]
#             end
#             for (i,j) in eachcol(nbr_table[:,1:(size(nbr_table,2) - Lx)] .+ model_geometry.lattice.N)
#                 H_t[i,j] += -tight_binding_model.t[1]
#                 H_t[j,i] += -tight_binding_model.t[1]
#             end  
#         else
#             for (i,j) in eachcol(nbr_table)
#                 H_t[i,j] += -tight_binding_model.t[1]
#                 if model_geometry.lattice.N > 2
#                     H_t[j,i] += -tight_binding_model.t[1]
#                 else
#                 end
#             end
#             for (i,j) in eachcol(nbr_table .+ model_geometry.lattice.N)    
#                 H_t[i,j] += -tight_binding_model.t[1]
#                 if model_geometry.lattice.N > 2
#                     H_t[j,i] += -tight_binding_model.t[1]
#                 else
#                 end
#             end
#         end
#         # next nearest neighbor hopping
#         if tight_binding_model.t[2] != 0.0
#             nbr_table_p = build_neighbor_table(bonds[2],
#                                             model_geometry.unit_cell,
#                                             model_geometry.lattice)
#             if Lx == 2 && Ly ==2 
#                 for (i,j) in eachcol(nbr_table_p)
#                     H_tp[i,j] += tight_binding_model.t[2]/2
#                 end
#                 for (i,j) in eachcol(nbr_table_p .+ model_geometry.lattice.N)    
#                     H_tp[i,j] += tight_binding_model.t[2]/2
#                 end
#             else
#                 for (i,j) in eachcol(nbr_table_p)
#                     H_tp[i,j] += tight_binding_model.t[2]
#                     H_tp[j,i] += tight_binding_model.t[2]
#                 end
#                 for (i,j) in eachcol(nbr_table_p .+ model_geometry.lattice.N)    
#                     H_tp[i,j] += tight_binding_model.t[2]
#                     H_tp[j,i] += tight_binding_model.t[2]
#                 end
#             end
#         end
#     end

#     if !("μ_BCS" in parameters_to_optimize)
#         return H_t + H_tp + LinearAlgebra.Diagonal(μ_vec)
#     else
#         return H_t + H_tp
#     end
# end


"""
    build_variational_terms( determinantal_parameters::DeterminantalParameters ) 

Constructs a 2⋅n⋅N by 2⋅n⋅N matrices to be added to the non-interacting tight binding
Hamiltonian for each variational parameter. Returns a vector of the sum of
matrices and a vector of individual matrix terms.

"""
function build_variational_terms(determinantal_parameters)
    dims = model_geometry.unit_cell.n*model_geometry.lattice.N
    vparam_map = map_determinantal_parameters(determinantal_parameters) 
    Hs = zeros(AbstractFloat, 2*dims, 2*dims)
    Hd = zeros(AbstractFloat, 2*dims, 2*dims)    
    Ha = zeros(AbstractFloat, 2*dims, 2*dims)    
    Hc = zeros(AbstractFloat, 2*dims, 2*dims) 
    Hμ = zeros(AbstractFloat, 2*dims, 2*dims)       
    Hcs = zeros(AbstractFloat, 2*dims, 2*dims)
    Hss = zeros(AbstractFloat, 2*dims, 2*dims)    
    H_vpars = []
    V = []
    # check for s-wave order
    if haskey(vparam_map, "Δs") == true
        @assert pht == true
        bA = zeros(AbstractFloat, dims, dims)
        bD = zeros(AbstractFloat, dims, dims)
        Δ_vec_bB = Vector{AbstractFloat}(undef, dims)
        Δ_vec_bC = Vector{AbstractFloat}(undef, dims)
        for i in 1:dims
            Δ_vec_bB[i] = 1
            Δ_vec_bC[i] = 1
        end
        bB = LinearAlgebra.Diagonal(Δ_vec_bB)
        bC = LinearAlgebra.Diagonal(Δ_vec_bC)
        Vs = Matrix([bA bB; bC bD])
        Hs += vparam_map["Δs"][1]*Vs
        push!(H_vpars,Hs)
        push!(V, Vs)
    end
    # check for d-wave order
    if haskey(vparam_map, "Δd") == true
        @assert pht == true
    end
    # check for antiferromagnetic (Neél) order
    if haskey(vparam_map, "Δa") == true
        afm_vec = Vector{AbstractFloat}(undef, 2*dims)
        if pht == true
            # particle-hole transformed Neél order
            for i in 1:dims
                for j in dims+1:2*dims
                    for k in 1:ndims(model_geometry.lattice)
                        loc_up = site_to_loc(i,model_geometry.unit_cell,model_geometry.lattice)[1]
                        afm_vec[i] = (-1)^(sum(loc_up[k]))*(-1)
                        loc_dwn = site_to_loc(j-model_geometry.lattice.N,model_geometry.unit_cell,model_geometry.lattice)[1]
                        afm_vec[j] = (-1)^(sum(loc_dwn[k]))
                    end
                end 
            end
        else
            # Neél order
            for i in 1:dims
                for j in dims+1:2*dims
                    for k in 1:ndims(model_geometry.lattice)
                        loc_up = site_to_loc(i,model_geometry.unit_cell,model_geometry.lattice)[1]
                        afm_vec[i] = (-1)^(sum(loc_up[k]))*(-1)
                        loc_dwn = site_to_loc(j-model_geometry.lattice.N,model_geometry.unit_cell,model_geometry.lattice)[1]
                        afm_vec[j] = (-1)^(sum(loc_dwn[k]))*(-1) 
                    end
                end 
            end
        end
        Va = LinearAlgebra.Diagonal(afm_vec)
        Ha += vparam_map["Δa"][1]*Va
        push!(H_vpars,Ha)
        push!(V, Va)
    end
    # check for uniform charge order
    if haskey(vparam_map, "Δc") == true
        cdw_vec = Vector{AbstractFloat}(undef, 2*dims)
        if pht == true
            # particle-hole transformed charge density wave
            for i in 1:dims
                for j in dims+1:2*dims
                    for k in 1:ndims(model_geometry.lattice)
                        loc_up = site_to_loc(i,model_geometry.unit_cell,model_geometry.lattice)[1]
                        cdw_vec[i] = (-1)^(sum(loc_up[k]))
                        loc_dwn = site_to_loc(j-model_geometry.lattice,model_geometry.unit_cell,model_geometry.lattice)[1]
                        cdw_vec[j] = (-1)^(sum(loc_dwn[k]))
                    end
                end 
            end
        else
            # charge density wave
            for i in 1:dims
                for j in dims+1:2*dims
                    for k in 1:ndims(model_geometry.lattice)
                        loc_up = site_to_loc(i,model_geometry.unit_cell,model_geometry.lattice)[1]
                        cdw_vec[i] = (-1)^(sum(loc_up[k]))
                        loc_dwn = site_to_loc(j-model_geometry.lattice,model_geometry.unit_cell,model_geometry.lattice)[1]
                        cdw_vec[j] = (-1)^(sum(loc_dwn[k]))*(-1) 
                    end
                end 
            end
        end
        Vc = LinearAlgebra.Diagonal(cdw_vec)
        Hc += vparam_map["Δc"][1]*Vc
        push!(H_vpars,Hc)
        push!(V, Vc)
    end
    # check for BCS chemical potential
    if haskey(vparam_map, "μ_BCS") == true
        μ_vec = Vector{AbstractFloat}(undef, 2*dims)
        if pht == true
            # particle-hole transformed chemical potential
            for i in 1:dims
                for j in dims+1:2*dims
                    μ_vec[i] = -1
                    μ_vec[j] = 1
                end 
            end
        else
            # chemical potential
            for i in 1:dims
                for j in dims+1:2*dims
                    μ_vec[i] = -1
                    μ_vec[j] = -1
                end 
            end
        end
        Vμ = LinearAlgebra.Diagonal(μ_vec)
        Hμ += vparam_map["μ_BCS"][1]*Vμ
        push!(H_vpars,Hμ)
        push!(V, Vμ)
    end
    # check for charge stripes
    if haskey(vparam_map, "Δcs") == true
        # vector to be transformed to matrix
        cs_vec = Vector{AbstractFloat}(undef, 2*dims)
        for i in 1:dims
            for j in dims+1:2*dims
                cs_vec[i] = 1
                cs_vec[j] = 1
            end 
        end
        Vcs = LinearAlgebra.Diagonal(cs_vec)
        Hcs += vparam_map["Δcs"][1]*Vcs
        push!(H_vpars,Hcs)
        push!(V, Vcs)
    end
    # check for spin stripes
    if haskey(vparam_map, "Δss") == true
        # vector to be transformed to matrix
        ss_vec = Vector{AbstractFloat}(undef, 2*dims)
        for i in 1:dims
            for j in dims+1:2*dims
                ss_vec[i] = 1
                ss_vec[j] = -1
            end 
        end
        Vss = LinearAlgebra.Diagonal(ss_vec)
        Hss += vparam_map["Δcs"][1]*Vss
        push!(H_vpars,Hss)
        push!(V, Vss)
    end
    # pd_vec = Vector{AbstractFloat}(undef, 2*norbs*L)
    #     for i in 2:3:(2*norbs*L)-1
    #         for j in 3:3:2*norbs*L
    #             pd_vec[i] = 1
    #                 pd_vec[j] = 1
    #         end
    #     end
       
    #     return LinearAlgebra.Diagonal(pd_vec)
    return [sum(H_vpars),V]
end



"""
    build_mean_field_hamiltonian() 

Constructs a matrix by combining the non-interacting Hamiltonian with
matrix of variational terms.

"""
function build_mean_field_hamiltonian()
    if verbose
        println("Building mean-field Hamiltonian...")
    end
    return build_tight_binding_model(tight_binding_model) + build_variational_terms(determinantal_parameters)[1], build_variational_terms(determinantal_parameters)[2]
end


"""
    build_determinantal_state() 

Returns initial energies ε_init, matrix M, and Slater matrix D in 
the many-particle configuration basis with associated initial energies. 

"""
function build_determinantal_state()
    # diagonalize Hamiltonian
    ε, U = diagonalize(H_mf) 
    if is_openshell(ε,Np) == true
        if verbose 
            println("WARNING! Open shell detected")
        end
    else
        if verbose
            println("Generating shell...")
        end
    end
    # store energies and M matrix
    ε₀ = ε[1:Np]  
    M = hcat(U[:,1:Np])

    # build Slater determinant
    D = zeros(AbstractFloat, Np, Np)
    while true
        pconfig = generate_initial_electron_configuration()
        D = M[findall(x -> x == 1, pconfig), :]
        
        if is_invertible(D)
            # write matrices to file
            if write == true
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
    get_Ak_matrices( V::Vector{Matrix{AbstractFloat}}, U::Matrix{AbstractFloat}, ε::Vector{AbstractFloat}, model_geometry::ModelGeometry ) 
    
Returns variational parameter matrices Aₖ from the corresponding Vₖ. Computes Qₖ = (U⁺VₖU)_(ην) / (ε_η - ε_ν), for η > Nₚ and ν ≤ Nₚ and is 0 otherwise
(η and ν run from 1 to 2L)

"""
function get_Ak_matrices(V, U, ε, model_geometry)
    if verbose
        println("Building A matrices...")
    end

    dims = model_geometry.unit_cell.n * model_geometry.lattice.N
    A = Vector{Matrix{AbstractFloat}}()
    adjU = adjoint(U)

    # Generate matrix of perturbation theory energies
    pt_energies = 1.0 ./ ε[1:2*dims]' .- ε[1:2*dims]

    # iterate over each Vₖ for each variational parameter
    for i in V
        push!(A, U * (adjU * i * U * pt_energies) * adjU)
    end

    return A
end

# end # of module


















