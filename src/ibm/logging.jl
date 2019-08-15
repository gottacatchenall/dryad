
function logging(instance::ibm, gen::Int64, genomeDF::DataFrame, popDF::DataFrame, log_pt::Int64)
    allele_freq_map::Array{Array{Float64,2}} = get_allele_freq_map(instance)
    id::Int64 = 1

    log_genome_data(allele_freq_map, genomeDF, gen)
    log_pop_data(instance, allele_freq_map, popDF, gen)
end

function log_pop_data(instance::ibm, allele_freq_map::Array{Array{Float64,2}}, popDF::DataFrame, gen::Int64)
    n_pops = length(instance.fitness_map)
    pop_cts::Array{Float64, 1} = get_cts(instance)

    for p = 1:n_pops
        push!(popDF.gen, gen)
        push!(popDF.pop, p)
        push!(popDF.w_mean, instance.fitness_map[p])

        prop_k::Float64 = pop_cts[p] / instance.mp.populations[p].k
        push!(popDF.prop_of_k, prop_k)

        mean_poly::Float64 = get_mean_poly_ct_in_pop(p, allele_freq_map)
        push!(popDF.mean_poly_ct, mean_poly)
    end
end

function get_mean_poly_ct_in_pop(p::Int64, allele_freq_map::Array{Array{Float64,2}})
    n_loci::Int64 = length(allele_freq_map)

    sum::Float64 = 0
    for l = 1:n_loci
        n_alleles::Int64 = size(allele_freq_map[l])[2]
        for i = 1:n_alleles
            if (allele_freq_map[l][p,i] > 0)
                sum += 1.0
            end
        end
    end
    mean_poly::Float64 = sum / n_loci
    return mean_poly
end

function log_genome_data(allele_freq_map::Array{Array{Float64,2}}, genomeDF::DataFrame, gen::Int64)
    n_loci::Int64 = length(allele_freq_map)
    for l = 1:n_loci
        state::Array{Float64,2} = (allele_freq_map[l])
        gst::Float64 = calc_gst(state)
        jostd::Float64 = calc_jost_d(state)
        mean_poly_ct::Float64 = calc_mean_poly_ct(state)
        # log poly ct
        push!(genomeDF.mean_poly_ct_per_pop, mean_poly_ct)
        push!(genomeDF.gen, gen)
        push!(genomeDF.gst, gst)
        push!(genomeDF.jostd, jostd)
        push!(genomeDF.locus, l)
    end

end

function get_allele_freq_map(instance::ibm)
    n_pops::Int64 = length(instance.mp.populations)
    n_loci::Int64 = instance.g.n_loci
    allele_freq_map::Array{Array{Float64,2}} = fill(zeros(1,1), n_loci)

    max_n_indiv::Int64 = length(instance.population_map)

    for l = 1:n_loci
        n_poly_this_locus::Int64 = length(instance.g.allele_dict[l])
        allele_freq_map_this_locus::Array{Float64, 2} = zeros(n_pops, n_poly_this_locus)

        for i = 1:max_n_indiv
            pop::Int64 = instance.population_map[i]
            if (pop != 0)
                for h = 1:2
                    index::Int64 = get_allele_index(instance.g.allele_dict[l], instance.genotypes[i, l, h])
                    allele_freq_map_this_locus[pop, index] += 1
                end
            end
        end
        allele_freq_map[l] = allele_freq_map_this_locus
    end


    pop_cts::Array{Float64, 1} = get_cts(instance)

    for l = 1:n_loci
        n_poly_this_locus::Int64 = length(instance.g.allele_dict[l])
        for p = 1:n_pops
            eff_pop_size::Float64 = pop_cts[p]*2
            for al = 1:n_poly_this_locus
                if (eff_pop_size > 0)
                    allele_freq_map[l][p, al] = allele_freq_map[l][p, al] / (eff_pop_size)
                end
            end
        end
    end

    return allele_freq_map
end


function get_allele_index(allele_list::Array{Float64, 1}, allele::Float64)
    index::Array{Int64} = findall(x -> x == allele, allele_list)
    return index[1]
end



function update_ibm_metadata(metadata::DataFrame,id_ct::Int64, m::Float64, s::Float64, k::Int64, n_ef::Int64, n_chromo::Int64, genome_length::Float64)
    push!(metadata.id, id_ct)
    push!(metadata.m, m)
    push!(metadata.s, s)
    push!(metadata.k, k)
    push!(metadata.n_ef, n_ef)
    push!(metadata.n_chromo, n_chromo)
    push!(metadata.genome_length, genome_length)
end

function update_fits_metadata(metadata::DataFrame, id_ct::Int64, m::Float64, k::Int64, init_poly_ct::Float64)
    push!(metadata.id, id_ct)
    push!(metadata.m, m)
    push!(metadata.k, k)
    push!(metadata.init_poly_ct, init_poly_ct)
end
