
mutable struct population
    x::Float64
    y::Float64
    k::Float64
    efs::Array{Float64}
end

mutable struct dispersal_kernel
    n_pops::Int64
    D::Array{Float64,2}
end

mutable struct metapop
    populations::Array{population}
    diskern::dispersal_kernel
    n_indivs::Int64
end

mutable struct fits
    mp::metapop
    n_pops::Int64
    n_alleles::Int64
    n_gen::Int64
    log_freq::Int64
    migration_rate::Float64
    mutation_rate::Float64
    ct_map::Array{Float64}
    rng::MersenneTwister
    fits(mp::metapop; n_alleles::Int64=5, mutation_rate::Float64 = 10^-5, migration_rate::Float64=0.01, n_gen=300, log_freq=20, rseed=1) = new(mp, length(mp.populations), n_alleles, n_gen, log_freq, migration_rate, mutation_rate, zeros(Float64, length(mp.populations), n_alleles), MersenneTwister(rseed))
end


mutable struct chromosome
    n_loci::Int64
    length::Float64
    selection_str::Array{Float64}
    ef_map::Array{Int64}
    poly_ct::Array{Int64}
end

mutable struct genome
    n_loci::Int64
    n_haplo::Int64
    mutation_rate::Float64
    allele_dict::Array{Array{Float64,1},1}
    chromosomes::Array{chromosome}
end

mutable struct ibm
    mp::metapop
    g::genome
    dk::dispersal_kernel
    genotypes::Array{Float64, 3}
    population_map::Array{Int64, 1}
    fitness_map::Array{Float64,1}
    rng::MersenneTwister
end

macro ecological_selection() return :(1) end

macro background_selection() return :(2) end

macro neutral_selection() return :(3) end

macro extant() return :(1) end
macro extinct() return :(0) end

macro rect_diskern() return :(1) end
macro gauss_diskern() return :(2) end
macro exponential_diskern() return :(3) end
macro uniform_diskern() return :(4) end
