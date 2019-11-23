include("../types.jl")

function init_random_metapop(;num_indivs::Int64=1000, num_populations::Int64=20, selection_type=@background_selection, diskern_type=@uniform_diskern, ibd_decay = 0.1, n_efs::Int64=1)::metapop
    pops::Array{population} = get_random_populations(num_indivs, num_populations, selection_type, n_efs)

    diskern::dispersal_kernel = init_uniform_diskern(num_populations)

    if diskern_type == @exponential_diskern
        diskern = init_exponential_diskern(pops, ibd_decay)
    elseif diskern_type == @gauss_diskern
        diskern = init_gauss_diskern(pops, ibd_decay)
    elseif diskern_type == @rect_diskern
        diskern = init_rect_diskern(pops, ibd_decay)
    end

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

function init_exponential_diskern(pops::Array{population}, diskern_strength)::dispersal_kernel
    num_populations = length(pops)

    D::Array{Float64, 2} = zeros(num_populations, num_populations)
    for i = 1:num_populations
        x1 = pops[i].x
        y1 = pops[i].y

        row_sum = 0.0
        for j = 1:num_populations
            if i != j
                x2 = pops[j].x
                y2 = pops[j].y

                dist = sqrt((x2-x1)^2 + (y2-y1)^2)

                kern = exp(-1*dist *diskern_strength)
                row_sum += kern
                D[i,j] = kern
            end
        end

        for j = 1:num_populations
            D[i,j] = D[i,j] / row_sum
        end
    end
    diskern::dispersal_kernel = dispersal_kernel(num_populations, D)
    return diskern
end

function init_gauss_diskern(pops::Array{population}, diskern_strength)::dispersal_kernel
    num_populations = length(pops)

    D::Array{Float64, 2} = zeros(num_populations, num_populations)
    for i = 1:num_populations
        x1 = pops[i].x
        y1 = pops[i].y

        row_sum = 0.0
        for j = 1:num_populations
            if i != j
                x2 = pops[j].x
                y2 = pops[j].y

                dist = sqrt((x2-x1)^2 + (y2-y1)^2)

                kern = exp(-1*dist*dist*diskern_strength*diskern_strength)
                row_sum += kern
                D[i,j] = kern
            end
        end

        for j = 1:num_populations
            D[i,j] = D[i,j] / row_sum
        end
    end
    diskern::dispersal_kernel = dispersal_kernel(num_populations, D)
    return diskern
end

function init_rect_diskern(pops::Array{population}, ibd_str)::dispersal_kernel
    num_populations = length(pops)

    D::Array{Float64, 2} = zeros(num_populations, num_populations)
    for i = 1:num_populations
        x1 = pops[i].x
        y1 = pops[i].y

        row_sum = 0.0
        for j = 1:num_populations
            if i != j
                x2 = pops[j].x
                y2 = pops[j].y

                dist = sqrt((x2-x1)^2 + (y2-y1)^2)

                if dist < (1.0/ibd_str)
                    row_sum += 1
                    D[i,j] = 1
                end
            end
        end

        if row_sum > 0
            for j = 1:num_populations
                D[i,j] = D[i,j] / row_sum
            end
        end
    end
    diskern::dispersal_kernel = dispersal_kernel(num_populations, D)
    return diskern
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
