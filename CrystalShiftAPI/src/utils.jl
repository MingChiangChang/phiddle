using CrystalShift: get_crystal

function CrystalPhase_helper(dict::Dict{String, Any})
    crystal = get_crystal([dict["a"], dict["b"], dict["c"], dict["α"], dict["β"], dict["γ"]])
    peak_hkl = dict["peaks"]
    peaks = Vector{Peak}(undef, length(peak_hkl))
    for i in eachindex(peak_hkl)
        peaks[i] = Peak(peak_hkl[i]["h"], peak_hkl[i]["k"], peak_hkl[i]["l"], 0., peak_hkl[i]["I"])
    end

    CrystalPhase(crystal, crystal, peaks, dict["param_num"],
                 dict["id"], dict["name"], dict["act"], dict["σ"],
                 FixedPseudoVoigt(dict["peak_profile"]["α"]),
                 dict["norm_constant"])
end


function get_cs_with_names(cs::AbstractVector, names::AbstractVector)
    cs_idx = Int64[]
    for i in eachindex(names)
        idx = findfirst(x->x.name==names[i], cs)
        if !isnothing(idx)
            push!(cs_idx, idx)
        end
    end
    return cs_idx
end

function get_opt_stn_from_dict(dict::Dict)
    optimization_mode = get_opt_mode_enum(dict["optimization_mode"]) 
    optimize_method = get_opt_method_enum(dict["method"])
    mean_θ = convert.(Float64, dict["mean_θ"])
    std_θ = convert.(Float64, dict["std_θ"])
    OptimizationSettings{Float64}(dict["std_noise"], mean_θ, std_θ,
                                  dict["maxiter"], dict["regularization"], optimize_method, dict["objective"],
                                  optimization_mode, dict["em_loop_num"], dict["λ"], dict["verbose"],
                                  dict["tol"])
end

function get_opt_mode_enum(opt_mode::String)
    if opt_mode == "Simple"
        Simple
    elseif opt_mode == "EM"
        EM
    elseif opt_mode == "Uncer"
        WithUncer
    end
end


function get_opt_method_enum(opt_mode::String)
    if opt_mode == "LM"
        LM
    elseif opt_mode == "bfgs"
        bfgs
    end
end


function get_ts_stn_from_dict(dict::Dict, phases)
    opt_stn = get_opt_stn_from_dict(dict)
    if dict["default_phase_idx"] == 0
        return TreeSearchSettings{Float64}(dict["depth"],
                                           dict["k"],
                                           dict["amorphous"],
                                           dict["background"],
                                           dict["background_length"],
                                           opt_stn)
    end

    return TreeSearchSettings{Float64}(dict["depth"],
                                       dict["k"],
                                       dict["amorphous"],
                                       dict["background"],
                                       dict["background_length"],
                                       phases[dict["default_phase_idx"]],
                                       opt_stn)

end


function struct_instance_from_dict(struct_type::Type, dict::Dict{String, Any})
    # Convert keys to Symbols if they're Strings to match field names
    symbolized_dict = Dict(Symbol(k) => v for (k, v) in dict)
    # Use keyword argument splatting to pass the dictionary values to the constructor
    return struct_type(; symbolized_dict...)
end

function get_phase_with_name(cs::AbstractVector, name::String)
    return cs[findfirst(x->x.name==name, cs)]
end
