
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
    ct_map::Array{Float64}
    rng::MersenneTwister
    fits(mp::metapop; n_alleles::Int64=5, migration_rate::Float64=0.01, n_gen=300, log_freq=20, rseed=1) = new(mp, length(mp.populations), n_alleles, n_gen, log_freq, migration_rate, zeros(Float64, length(mp.populations), n_alleles), MersenneTwister(rseed))
end


macro ecological_selection() return :(1) end

macro background_selection() return :(2) end

macro neutral_selection() return :(3) end
