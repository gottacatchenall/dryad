
# ========================================================
#  compute
#  summary statisitics
#
# ========================================================
function calc_ht(state::Array{Float64,2})
    n_pops::Int64 = size(state)[1]
    n_alleles::Int64 = size(state)[2]

    outer_sum::Float64 = 0.0
    for al = 1:n_alleles
        inner_s::Float64 = 0.0
        for p = 1:n_pops
            eff_pop_size::Float64 = sum(state[p,:])
            if (eff_pop_size > 0)
                inner_s += (state[p,al] / eff_pop_size)
            end
        end

        inner_s = inner_s / n_pops
        outer_sum += (inner_s^2)
    end
    return(1.0 - outer_sum)
end

function calc_hs(state::Array{Float64,2})
    n_pops::Int64 = size(state)[1]
    n_alleles::Int64 = size(state)[2]

    s::Float64 = 0.0
    for p = 1:n_pops
        eff_pop_size::Float64 = sum(state[p,:])
        inner_s::Float64 = 0.0

        if (eff_pop_size > 0)
            inner_s = 0.0
            for al = 1:n_alleles
                inner_s += (state[p,al] / eff_pop_size)^2
            end
        else
            inner_s = 0.0
        end
        s += (1.0 - inner_s)

    end

    s = s / n_pops
    return s
end


function calc_jost_d(state::Array{Float64,2})
    n_pops::Int64 = size(state)[1]
    n_alleles::Int64 = size(state)[2]

    H_T::Float64 = calc_ht(state)
    H_S::Float64 = calc_hs(state)


    non_empty_pops::Int64 = 0
    for p = 1:n_pops
        eff_pop_size::Float64 = sum(state[p,:])
        if (eff_pop_size > 0)
            non_empty_pops += 1
        end
    end

    jostD::Float64 = 1.0
    if (non_empty_pops > 1)
        jostD = ((H_T - H_S)*non_empty_pops) / ((1.0 - H_S)*(non_empty_pops-1))
    end

    return(jostD)
end

function calc_gst(state::Array{Float64,2})
    n_pops::Int64 = size(state)[1]
    n_alleles::Int64 = size(state)[2]

    H_T::Float64 = calc_ht(state)
    H_S::Float64 = calc_hs(state)


    # check total polymorphism count
    # if all fixed, call it Gst 0


    if H_S == NaN
        return 0.0
    elseif H_T == 0 || H_S == 0 || H_T == NaN
        return 1.0
    end

    gst::Float64 = (H_T - H_S) / (H_T)

    return(gst)
end

function calc_mean_poly_ct(state::Array{Float64,2})::Float64
    n_pops::Int64 = size(state)[1]
    n_alleles::Int64 = size(state)[2]

    s::Float64 = 0.0

    non_empty_pops::Int64 = 0

    for p = 1:n_pops
        if (sum(state[p,:]) > 0)
            non_empty_pops += 1
        end
        for i = 1:n_alleles
            if (state[p,i] > 0)
                s += 1
            end
        end
    end
    mean_poly_ct::Float64 = s/non_empty_pops
    return mean_poly_ct
end
