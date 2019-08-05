include("../types.jl")

function init_random_metapop(;num_indivs::Int64=1000, num_populations::Int64=20)::metapop
    pops::Array{population} = get_random_populations(num_indivs, num_populations)
    diskern::dispersal_kernel = init_uniform_diskern(num_populations)
    mp::metapop = metapop(pops,diskern,num_indivs)
    return(mp)
end

function get_random_populations(num_indivs::Int64, num_populations::Int64)::Array{population}
    pops::Array{population} = []
    k::Float64 = num_indivs/num_populations
    for p = 1:num_populations
        tmp::population = population(rand(Uniform()), rand(Uniform()), k)
        push!(pops, tmp)
    end
    return pops
end

function init_uniform_diskern(num_populations::Int64)::dispersal_kernel
    val::Float64 = 1.0 / num_populations
    D::Array{Float64, 2} = zeros(num_populations, num_populations)
    for i = 1:num_populations
        for j = 1:num_populations
            if (i != j)
                D[i,j] = val
            end
        end
    end
    diskern::dispersal_kernel = dispersal_kernel(num_populations, D)
    return diskern
end
