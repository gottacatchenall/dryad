
mutable struct population
    x::Float64
    y::Float64
    k::Float64
    efs::Array{Float64}
    population(x,y,k) = new(x,y,k,[1.0])
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
    migration_rate::Float64
    ct_map::Array{Float64}
    fits(mp::metapop, n_alleles::Int64, migration_rate::Float64) = new(mp, length(mp.populations), n_alleles, migration_rate, zeros(Float64, length(mp.populations), n_alleles))
end
