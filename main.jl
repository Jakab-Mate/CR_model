using Random, Dates, Plots, DifferentialEquations, LSODA, LinearAlgebra, MicrobiomeAnalysis

include("structs.jl")
include("generative_functions.jl")
include("sample_pool.jl")
include("equations.jl")
include("callbacks.jl")

seed = now().instant.periods.value
rng = MersenneTwister(seed)
t_span = (0.0, 2000)
n_resources = 10
n_species = 5
n_invaders = 100
t_inv = 25
t_inv_0 = 100
D, W_ba = generative_functions.create_metabolism(rng, n_resources=n_resources)
pool = generative_functions.create_species_pool(rng, D,
    n_families=10, 
    family_size=20, 
    dirichlet_hyper=100, 
    between_family_var=0.1, 
    inside_family_var=0.05, 
    h=2, 
    maintenance=1, 
    specialist=1, 
    generalist=1)

s = sample_pool(rng, pool, n_resources, n_species, n_invaders)

present_species = 1:20 # calculate_present_species()

phi = 0.1
eta = 0.05
tau = ones(Float64, n_resources)
alpha = vcat([100.0], zeros(Float64, n_resources-1))
C = s.C 
n_reactions = s.n_reactions
n_splits= s.n_splits
m = s.m

u0 = vcat(s.species_abundance, s.resource_abundance)
print("\n Shape u0:", size(u0), "Type:", eltype(u0))
params = (n_species+n_invaders, n_resources, present_species, C, D, W_ba, n_reactions, n_splits, m, phi, eta, tau, alpha)

prob = ODEProblem(equations, u0, t_span, params)
cb = create_callbacks(t_inv, t_inv_0, n_invaders, n_species)

solution = solve(prob, KenCarp4(autodiff=false), abstol = 1e-10, reltol = 1e-10; saveat=1, callback=cb) #, lsoda(), save_everystep = false
# solution = solve(prob, alg_hints=[:stiff], abstol = 1e-8, reltol = 1e-8; saveat=1, callback=cb) #, lsoda(), save_everystep = false
# solution = solve(prob, Tsit5(), callback=cb)
t = solution.t
Xapp = hcat(map(u -> u[1:(n_species+n_invaders)], solution.u)...)
print("\nDims\n", size(Xapp))

assays = OrderedDict{String, AbstractArray}("sim" => Xapp);

rowdata = DataFrame(
    name = ["strain$i" for i in 1:(n_species+n_invaders)],
    genus = ["g$i" for i in 1:(n_species+n_invaders)],
    species = ["s$i" for i in 1:(n_species+n_invaders)]
);

coldata = DataFrame(
    name = ["t$i" for i in 1:length(t)],
    condition = rand(["lake", "ocean", "river"], length(t)),
    time = 1:length(t)
);

se = SummarizedExperiment(assays, rowdata, coldata);

ginisimpson_output = ginisimpson(se, "sim")
print("\nginis length: ", size(ginisimpson_output))

gr()
plot(solution)
savefig("/home/jakab/vs_code/CR_model/figures/myplot.png")

