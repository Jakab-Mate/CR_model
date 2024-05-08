using Random, Dates, Plots, DifferentialEquations, LSODA, LinearAlgebra, MicrobiomeAnalysis

include("structs.jl")
include("generative_functions.jl")
include("sample_pool.jl")
include("equations.jl")

seed = now().instant.periods.value
rng = MersenneTwister(seed)
n_resources = 5
n_species = 10
D, W_ba = generative_functions.create_metabolism(rng, n_resources=n_resources)
pool = generative_functions.create_species_pool(rng, D,
    n_families=5, 
    family_size=10, 
    dirichlet_hyper=100, 
    between_family_var=0.1, 
    inside_family_var=0.05, 
    h=2, 
    maintenance=1, 
    specialist=1, 
    generalist=1)

s = sample_pool(rng, pool, n_resources, n_sampled=n_species)

present_species = 1:10 # calculate_present_species()
#counter = Ref(0)

phi = 0.1
eta = 0.05
tau = ones(Float64, n_resources)
alpha = vcat([100.0], zeros(Float64, n_resources-1))
C = s.C 
n_reactions = s.n_reactions
n_splits= s.n_splits
m = s.m

u0 = vcat(s.species_abundance, s.resource_abundance)
t_span = (0.0, 100)
params = (n_species, n_resources, present_species, C, D, W_ba, n_reactions, n_splits, m, phi, eta, tau, alpha)

prob = ODEProblem(equations, u0, t_span, params)

# solution = solve(prob,alg_hints=[:stiff], abstol = 1e-12, reltol = 1e-12; saveat=1) #, lsoda(), save_everystep = false
solution = solve(prob, lsoda())
t = solution.t
Xapp = hcat(map(u -> u[1:n_species], solution.u)...)
print("\nDims\n", size(Xapp))

assays = OrderedDict{String, AbstractArray}("sim" => Xapp);

rowdata = DataFrame(
    name = ["strain$i" for i in 1:n_species],
    genus = ["g$i" for i in 1:n_species],
    species = ["s$i" for i in 1:n_species]
);

coldata = DataFrame(
    name = ["t$i" for i in 1:length(t)],
    condition = rand(["lake", "ocean", "river"], length(t)),
    time = 1:length(t)
);

se = SummarizedExperiment(assays, rowdata, coldata);

ginisimpson_output = ginisimpson(se, "sim")
print("ginis", ginisimpson_output)
print("\nginis length: ", size(ginisimpson_output))

gr()
plot(solution)
savefig("myplot.png")

