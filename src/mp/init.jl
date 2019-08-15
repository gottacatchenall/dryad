include("../types.jl")

function init_random_metapop(;num_indivs::Int64=1000, num_populations::Int64=20, selection_type=@background_selection, n_efs::Int64=1)::metapop
    pops::Array{population} = get_random_populations(num_indivs, num_populations, selection_type, n_efs)
    diskern::dispersal_kernel = init_uniform_diskern(num_populations)
    mp::metapop = metapop(pops,diskern,num_indivs)
    return(mp)
end

# this is where to add diff types of selection as a keyword argument
function get_random_populations(num_indivs::Int64, num_populations::Int64, selection_type, n_efs::Int64)::Array{population}
    pops::Array{population} = []
    k::Float64 = num_indivs/num_populations

    for p = 1:num_populations
        efs = init_efs(selection_type, n_efs)
        tmp::population = population(rand(Uniform()), rand(Uniform()), k, efs)
        push!(pops, tmp)
    end
    return pops
end

function init_efs(selection_type, n_efs)
    efs::Array{Float64} = []
    if selection_type == @ecological_selection
        efs = rand(Uniform(), n_efs)
    elseif selection_type == @background_selection
        efs = fill(1.0, n_efs)
    elseif selection_type == @neutral_selection
        efs = []
    end
    return efs
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

function set_mp_total_k(mp::metapop, k::Int64)
    n_pops::Int64 = length(mp.populations)
    base::Float64 = (k/n_pops)
    for p = 1:n_pops
        mp.populations[p].k = base
    end
    mp.n_indivs = k
end
