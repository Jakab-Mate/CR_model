using LinearAlgebra

function equations(u, p, t)
    n_species, n_resources, present_species, C, D, W_ba, n_reactions, n_splits, m, phi, eta, tau, alpha = p
    #du[n_species+1:end] = (alpha .- u[n_species+1:end]) ./ tau
    #du[1:n_species] = zeros(Float64, n_species)
    species = zeros(Float64, n_species)
    consumption = zeros(Float64, n_resources)
    production = zeros(Float64, n_resources)

    # print("\n Shape u:", size(u), "Type:", eltype(u))
    # print("\n Shape C:", size(C), "Type:", eltype(C))
    # print("\n Shape D:", size(D), "Type:", eltype(D))
    # print("\n Shape tau:", size(tau), "Type:", eltype(tau))
    # print("\n Shape alpha:", size(alpha), "Type:", eltype(alpha))
    # print("\n Shape consumption:", size(consumption), "Type:", eltype(consumption))
    # print("\n resources", u[n_species+1:end])
    # print("\n Spahe D*res:", size(D .* u[n_species+1:end]), "Type:", eltype(D .* u[n_species+1:end]))
    # print("\n Shape res", size(u[n_species+1:end]), "Type:", eltype(u[n_species+1:end]))
    # print("\n D*res", print(D .* u[n_species+1:end]), "Type:", eltype(D .* u[n_species+1:end]))
    # print("\n t: ", t)

    for i in 1:n_species
        if i == 4
            if t > 100
                # print("\n W_ba=", W_ba, "\n")
                # print("\n D= ", D, "\n")
                # print("C*W_ba= ", C[:,:,i] .* W_ba, "\n")
                # print("C*W_ba*res= ", C[:,:,i] .* W_ba .* u[n_species+1:end], "\n")
                # print("sum= ", sum(C[:,:,i] .* W_ba .* u[n_species+1:end]), "\n")
            end
        end
        species[i] = (sum(C[:,:,i] .* W_ba .* u[n_species+1:end]') - (m[i] + phi * n_reactions[i] + eta * n_splits[i])) * u[i]
        #print(size(C[:,:,i]))
        #print(size(C[:,:,i] .* W_ba))
        #print(size(C[:,:,i] .* W_ba .* u[n_species+1:end]))
        #print(typeof(species[i]))
        #print(size(sum(C[:,:,i] .* W_ba .* u[n_species+1:end])))

        consumption += dropdims(sum(C[:, :, i] .* u[i] .* u[n_species+1:end]', dims=1)', dims=2)
        production += dropdims(sum(D .* C[:, :, i] .* u[i] .* u[n_species+1:end]', dims=2), dims=2)

        #print("\n Shape prod", size(production), "Type:", eltype(production))
        #print("\n Shape con", size(consumption), "Type:", eltype(consumption))

    end

    resources = (alpha .- u[n_species+1:end]) ./ tau - consumption + production
    
    #print("\n Shape res", size(resources), "Type:", eltype(resources))
    #print("\n Shape spec", size(species), "Type:", eltype(species))
    #print("\n Shape ret", size(vcat(species, resources)), "Type:", eltype(vcat(species, resources)))
    return vcat(species, resources)
end



function diffeq_wrapper(t, u, du)
    du .= equations2!(t, u, n_species, present_species, C, D, W_ba, n_reactions, n_splits, m, phi, eta, tau, alpha)
end