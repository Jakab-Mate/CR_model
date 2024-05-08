module generative_functions

include("helper_functions.jl")
include("structs.jl")

using Distributions, StatsBase, Random

export create_species_pool

function create_metabolism(rng::Any; n_resources::Int64=10, n_levels::Int64=5, energy_yields::String="Uniform_1")

    levels::Vector{Int64} = vcat(1:n_levels, rand(rng, 2:n_levels-1, n_resources - n_levels))
    levels = sort(levels)
    print(levels, "\n")

    
    W = zeros(Float64, n_resources)
    for b in n_resources-1:-1:1
        if levels[b] == levels[b+1]
            W[b] = W[b+1]
        else
            if energy_yields == "Uniform_1"
                W[b] = 2 * W[b+1] + 1
            elseif energy_yields == "Random"
                W[b] = 2 * W[b+1] + rand() * 2
            end
        end
    end
    print("\n W: ", W)

    D = zeros(Int64, n_resources, n_resources)
    W_ba = zeros(Float64, n_resources, n_resources)
    for i in 1:n_resources
        for j in 1:n_resources
            if levels[i] > levels[j]
                D[i, j] = 2 ^ (levels[i] - levels[j])
                W_ba[i, j] = W[j] - W[i] * D[i, j]
                print("\n i: ", i)
                print("\n j: ", j)
                print("\n levels[i]: ", levels[i])
                print("\n levels[j]: ", levels[j])
                print("\n D_ij: ", D[i, j])
                print("\n W[i]: ", W[i])
                print("\n W[j]: ", W[j])
                print("\n W_ba[i, j]: ", W_ba[i, j])
            end
        end
    end
    
    print(typeof(D), typeof(W_ba))
    return D, W_ba
end

export create_metabolism

function create_species_pool(rng::Any, D::Matrix; n_families::Int64=5, 
    family_size::Int64=100, dirichlet_hyper::Real=100, between_family_var::Real=0.1, inside_family_var::Real=0.05, 
    h::Real=1, maintenance::Real=0.1, specialist::Real=1, generalist::Real=1)

    n_resources::Int64 = length(D[:, 1])
    
    specialist_prob::Int64 = round(Int64, specialist / (generalist + specialist) * 100)
    choices = vcat(fill("Spec", specialist_prob), fill("No", 100-specialist_prob))
    spec_gen = [rand(rng, choices) for _ in 1:n_families]

    pool = zeros(n_resources, n_resources, n_families * family_size)
    family_ids = Array{Int64}(undef, n_families * family_size)
    m = Array{Float64}(undef, n_families * family_size)
    n_reactions = Array{Int64}(undef, n_families * family_size)
    n_splits = Array{Float64}(undef, n_families * family_size)

    between_family_m_dist = Normal(maintenance, between_family_var)

    for family in 1:n_families
        family_ids[(family-1) * family_size + 1 : family * family_size] .= family
        family_idx_start = (family-1) * family_size  ### needed?
        value = rand(rng, between_family_m_dist)
        while value <= 0
            value = rand(rng, between_family_m_dist)
        end
        family_m = value
        inside_family_m_dist = Normal(family_m, inside_family_var)

        if spec_gen[family] == "Spec"
            family_number_of_reactions = rand(rng, [4, 5])
        else
            family_number_of_reactions = rand(rng, [2, 3])
        end

        print("\n Family $family has $family_number_of_reactions reactions \n")
        n_reactions[(family-1) * family_size + 1 : family * family_size] .= family_number_of_reactions

        prior_values = normalize(rand!(rng, zeros(family_number_of_reactions))) .* dirichlet_hyper
        Dir_dist = Distributions.Dirichlet(prior_values)
        family_reaction_indices = sample_reaction_indices(rng, D, family_number_of_reactions)

        for species in 1:family_size
            value = rand(rng, inside_family_m_dist)
            while value <= 0
                value = rand(rng, inside_family_m_dist)
            end
            m[family_idx_start + species] = value
            species_values = normalize(rand(rng, Dir_dist), h=h)
            for reaction in 1:family_number_of_reactions
                pool[family_reaction_indices[reaction]..., family_size*(family-1)+species] = species_values[reaction]
            end

            n_splits[family_idx_start + species] = sum((D .- 1) .*  pool[:, :, family_size*(family-1)+species])
        end

    end

    return pool_struct(pool, family_ids, m, n_reactions, n_splits)

end

end