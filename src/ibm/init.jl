function init_ef_map(n_chromo::Int64, loci_per_chromo::Int64, fitness_loci::Array{Int64,2})
    ef_map::Array{Int64, 2} = zeros(n_chromo, loci_per_chromo)

    n_ef::Int64 = size(fitness_loci)[1]
    n_loci_per_ef::Int64 = size(fitness_loci)[2]

    n_loci::Int64 = n_ef*n_loci_per_ef

    locus_id::Int64 = 0
    chr_num::Int64 = 0
    chr_position::Int64 = 0

    for ef = 1:n_ef
        for lpef = 1:n_loci_per_ef
            locus_id = fitness_loci[ef,lpef]
            chr_num = convert(Int64, ceil(locus_id / loci_per_chromo))
            chr_position = (locus_id % loci_per_chromo) + 1
            ef_map[chr_num, chr_position] = ef
        end
    end

    return ef_map
end

function init_sel_strs(n_chromo::Int64, loci_per_chromo::Int64)
    sel_strs::Array{Float64, 2} = zeros(n_chromo, loci_per_chromo)
end

function init_poly_array(n_chromo::Int64, loci_per_chromo::Int64, init_poly_ct::Float64, random_poly::Bool)
    poly_array::Array{Int64, 2} = zeros(n_chromo, loci_per_chromo)
    if (random_poly)
        for c = 1:n_chromo
            for l = 1:loci_per_chromo
                poly_array[c,l] = rand(Poisson(init_poly_ct))
            end
        end


        return poly_array
    else
        poly_array = fill(init_poly_ct, n_chromo, loci_per_chromo)
        return poly_array
    end
end

function init_fitness_loci(n_loci::Int64, n_loci_per_ef::Int64, n_ef::Int64)::Array{Int64, 2}
    n_fitness_loci::Int64 = n_loci_per_ef*n_ef
    fitness_loci::Array{Int64} = randperm(n_loci)[1:n_fitness_loci]
    fit_loci_matrix::Array{Int64, 2} = zeros(n_ef, n_loci_per_ef)

    fit_index::Int64 = 1
    l_ct::Int64 = 1
    for ef = 1:n_ef
        for lpef = 1:n_loci_per_ef
            fit_loci_matrix[ef,lpef] = fitness_loci[fit_index]
            fit_index += 1
        end
        fit_loci_matrix[ef,:] = sort(fit_loci_matrix[ef,:])
    end
    return fit_loci_matrix
end

function init_random_genome(; n_loci::Int64=200, n_chromo::Int64=5, n_haplo::Int64=2, n_ef::Int64=1, n_loci_per_ef::Int64=20, init_poly_ct::Float64=5.0, mutation_rate::Float64=0.00001, genome_length::Float64=100.0, random_poly::Bool=false)::genome
    loci_per_chromo::Int64 = convert(Int64, floor(n_loci / n_chromo))

    fitness_loci::Array{Int64,2} = init_fitness_loci(n_loci, n_loci_per_ef, n_ef)
    ef_map::Array{Int64, 2} = init_ef_map(n_chromo, loci_per_chromo, fitness_loci)
    sel_strs::Array{Float64, 2} = init_sel_strs(n_chromo, loci_per_chromo)
    poly_ct::Array{Int64, 2} = init_poly_array(n_chromo, loci_per_chromo, init_poly_ct, random_poly )


    chromo_length = genome_length / n_chromo
    chromos::Array{chromosome, 1} = []

    # init new chromo
    for c = 1:n_chromo
        tmp_chromo::chromosome = chromosome(loci_per_chromo, chromo_length, sel_strs[c,:], ef_map[c,:], poly_ct[c,:])
        push!(chromos, tmp_chromo)
    end

    allele_dict::Array{Array{Float64,1},1} = fill([], n_loci)

    return genome(n_loci, n_haplo, mutation_rate, allele_dict, chromos)
end

function get_max_indivs_per_pop(mp::metapop)
    maxK::Int64 = 0
    for p = 1:length(mp.populations)
        thisK::Int64 = round(mp.populations[p].k)
        if (thisK > maxK)
            maxK = thisK
        end
    end
    return(2*maxK)
end


function init_ibm(mp::metapop, g::genome, b; min=0.0, max=1.0, init_condition=false, rs=1)::ibm
    max_n_indivs::Int64 = 5*b*mp.n_indivs
    n_loci::Int64 =  g.n_loci
    n_haplo::Int64 = g.n_haplo
    n_pops::Int64 = length(mp.populations)

    pops::Array{population} = mp.populations

    #dk::dispersal_kernel = init_uniform_diskernel(mp)
    df::dispersal_kernel = mp.diskern

    # pop map stores indecies to genome and fitness
    pop_map::Array{Int64, 1} = zeros(max_n_indivs)

    # split n_indivs evenly across pops
    ind_ct::Int64 = 1
    for p = 1:n_pops
        ni::Int64 = round(pops[p].k/10)
        for i = 1:ni
            pop_map[ind_ct] = p
            ind_ct += 1
        end
    end

    genotypes::Array{Float64, 3} = zeros(max_n_indivs, n_loci, n_haplo)
    fitness_map::Array{Float64, 1} = zeros(max_n_indivs)

    if (init_condition)
        init_fixed_genotypes(g, genotypes, ind_ct)
        return(ibm(mp,g,dk,genotypes, pop_map, fitness_map, MersenneTwister(rs)))

    else
        init_random_uniform_genotypes(g, genotypes, ind_ct, min=min, max=max)
        return(ibm(mp,g,dk,genotypes, pop_map, fitness_map,MersenneTwister(rs)))
    end
end

# write function to init genotypes all to 1.0
function init_fixed_genotypes(g::genome, genotypes::Array{Float64, 3},  n_indivs::Int64)
    n_loci_per_chromo::Int64 = g.chromosomes[1].n_loci
    n_loci::Int64 = size(genotypes)[2]
    # put data in genotypes

    # check fitness map, if fitness, set val = 1.0
    for l = 1:n_loci
        c::Int64 = convert(Int64, ceil(l/n_loci_per_chromo))
        index_in_chromo::Int64 = (l % n_loci_per_chromo)+1

        # store this in a thing, every time there is a mutation push it on to this
        #==if g.chromosomes[c].ef_map[index_in_chromo] != 0
            for i = 1:n_indivs
                genotypes[i, l, 1] = 1.0
                genotypes[i, l, 2] = 1.0
            end
            g.allele_dict[l] = [1.0]
        else
            this_locus_init_poly::Int64 = g.chromosomes[c].poly_ct[index_in_chromo]
            this_locus_alleles = rand(Uniform(), this_locus_init_poly)
            g.allele_dict[l] = this_locus_alleles

            for i = 1:n_indivs
                rindex::Array{Int64} = rand(DiscreteUniform(1, this_locus_init_poly), 2)
                genotypes[i, l, 1] = this_locus_alleles[rindex[1]]
                genotypes[i, l, 2] = this_locus_alleles[rindex[2]]
            end
        end
        ==#
        for i = 1:n_indivs
            genotypes[i, l, 1] = 1.0
            genotypes[i, l, 2] = 1.0
        end
        g.allele_dict[l] = [1.0]
    end
end

function init_random_uniform_genotypes(g::genome, genotypes::Array{Float64, 3},  n_indivs::Int64; min=0.0, max=1.0)
    n_loci_per_chromo::Int64 = g.chromosomes[1].n_loci
    n_loci::Int64 = size(genotypes)[2]
    # put data in genotypes
    for l = 1:n_loci
        c::Int64 = convert(Int64, ceil(l/n_loci_per_chromo))
        index_in_chromo::Int64 = (l % n_loci_per_chromo)+1
        this_locus_init_poly::Int64 = g.chromosomes[c].poly_ct[index_in_chromo]

        # store this in a thing, every time there is a mutation push it on to this
        # this

        this_locus_alleles = rand(Uniform(min, max), this_locus_init_poly)
        g.allele_dict[l] = this_locus_alleles

        for i = 1:n_indivs
            rindex::Array{Int64} = rand(DiscreteUniform(1, this_locus_init_poly), 2)
            genotypes[i, l, 1] = this_locus_alleles[rindex[1]]
            genotypes[i, l, 2] = this_locus_alleles[rindex[2]]
        end
    end
end

function init_uniform_diskernel(mp::metapop)::dispersal_kernel
    n_pops::Int64 = length(mp.populations)
    base::Float64 = 1.0 / (n_pops - 1)
    D::Array{Float64, 2} = fill(base, n_pops, n_pops)
    return dispersal_kernel(n_pops, D)
end
