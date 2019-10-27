function create_treatments(param_dict::Dict)
	names::Array{String} = []
	vals = []
	treatment_df = DataFrame(names...)
	treatment_df.treatment_ct = []

	for (param_name, param_values) in param_dict
		push!(names, param_name)
		push!(vals, param_values)
		treatment_df[!,  Symbol(param_name)] = []
	end

	ct = 1
	for p in Iterators.product(vals...)
		for (i, val) in enumerate(p)
			param = names[i]
			push!(treatment_df[!, Symbol(param)], val)
		end
		push!(treatment_df.treatment_ct, ct)
		ct += 1
	end

	return(treatment_df)
end
