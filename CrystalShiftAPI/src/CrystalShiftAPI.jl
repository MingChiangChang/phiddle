using CrystalShift
using CrystalShift: OptimizationSettings, FixedPseudoVoigt, get_fraction
using CrystalTree
using CrystalTree: TreeSearchSettings, Lazytree, get_probabilities
using HTTP
using Oxygen
using JSON
using BackgroundSubtraction
using Base.Threads
using CovarianceFunctions: EQ

include("utils.jl")

# Get everything calls through pyjulia in phiddle to be here
# If necessary, can strip the unneeded data in tree nodes and keep only the essentials
#
# TODO: 

mutable struct ServerState
    csv_file::String
    phases::AbstractVector
    fit_result::AbstractVector
    bg::AbstractVector
end

global server_state = ServerState("", Vector{CrystalPhase}(), [], [])

@get "/last_background" function ()
    return server_state.bg
end

@get "/peak_info" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    phase = get_phase_with_name(server_state.phases, dict["phase_name"])
    # if all(getproperty.(phase.peaks, :q) .==0)

    phase.peaks
end

@get "/peak_info/{idx}" function (req::HTTP.Request, idx::Int)
    server_state.phases[idx+1].peaks # Use python indexing
end

@get "/phase_names" function (req::HTTP.Request)
    # dict = JSON.parse(String(req.body))
    # phase = get_phase_with_name(server_state.phases, dict["phase_name"])
    [phase.name for phase in server_state.phases]
end

@put "/set_csv" function set_csv(req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    server_state.csv_file = dict["csv_file"]

    open(server_state.csv_file, "r") do f # Takes 1-2 ms
        global server_state.phases = CrystalPhase(f, 0.2, FixedPseudoVoigt(0.5))
    end

    return server_state.csv_file
end


@post "/label" function label(req::HTTP.Request)
    dict = JSON.parse(String(req.body))

    data = convert.(Float64, dict["data"])
    q = convert.(Float64, dict["q"])
    mean_θ = convert.(Float64, dict["mean_θ"])
    std_θ = convert.(Float64, dict["std_θ"])

    min_data = minimum(data)
    n = maximum(data) - min_data
    data .-= min_data
    data ./= n
    if dict["background_option"] == "MCBL"
        bg = mcbl(data, q, convert(Float64, dict["background_length"]))
        data .-= bg
        data[data .< 0] .= 1E-5
    else
        bg = zero(data)
    end

    ts_stn = get_ts_stn_from_dict(dict)
    lt = Lazytree(server_state.phases, q)
    res = search!(lt, q, data, ts_stn)

    res = reduce(vcat, res[2:end])
    probs = get_probabilities(res, q, data, dict["std_noise"], mean_θ, std_θ)
    idx = sortperm(probs, rev=true)
    server_state.fit_result = [res[i].phase_model for i in idx]  
    println(idx)
    println(server_state.fit_result[1])

    return res[idx], probs[idx], bg
end


@post "/label_with_phase" function label_with_phase(req::HTTP.Request)
    dict = JSON.parse(String(req.body))

    data = convert.(Float64, dict["data"])
    q = convert.(Float64, dict["q"])
    mean_θ = convert.(Float64, dict["mean_θ"])
    std_θ = convert.(Float64, dict["std_θ"])

    phase_idx = get_cs_with_names(server_state.phases, dict["names"])

    min_data = minimum(data)
    n = maximum(data) - min_data
    data .-= min_data
    data ./= n
    if dict["background_option"] == "MCBL"
        bg = mcbl(data, q, convert(Float64, dict["background_length"]))
        data .-= bg
        data[data .< 0] .= 1E-5
    else
        bg = zeros(data)
    end

    ts_stn = get_ts_stn_from_dict(dict)
    lt = Lazytree(server_state.phases[phase_idx], q)
    res = search!(lt, q, data, ts_stn)
    
    res = reduce(vcat, res[2:end])
    probs = get_probabilities(res, q, data, dict["std_noise"], mean_θ, std_θ)
    idx = sortperm(probs, rev=true)
    server_state.fit_result = [res[i].phase_model for i in idx]  

    return res[idx], probs[idx], bg
end


@post "/fit_with_phase" function fit_with_phase(req::HTTP.Request)
    dict = JSON.parse(String(req.body))

    data = convert.(Float64, dict["data"])
    q = convert.(Float64, dict["q"])

    phase_idx = get_cs_with_names(server_state.phases, dict["names"])

    min_data = minimum(data)
    n = maximum(data) - min_data
    data .-= min_data
    data ./= n
    if dict["background_option"] == "MCBL"
        bg = mcbl(data, q, convert(Float64, dict["background_length"]))
        data .-= bg
        data[data .< 0] .= 1E-5
        pm = PhaseModel(server_state.phases[phase_idx], nothing, nothing)
    elseif dict["background_option"] == "Default"
        _bg = BackgroundModel(q, EQ(), dict["background_length"], 10.)
        pm = PhaseModel(server_state.phases[phase_idx], nothing, _bg)
        bg = zero(q)
    end

    opt_stn = get_opt_stn_from_dict(dict)
    result = optimize!(pm, q, data, opt_stn)
    server_state.fit_result = [result]
    return result, [1.], bg
end

@get "/evaluate" function evaluate(req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    q = convert.(Float64, dict["q"])
    if !(dict["phase_name"] isa AbstractVector)
        dict["phase_name"] = [dict["phase_name"]]
    end
    cps = server_state.phases[get_cs_with_names(server_state.phases, dict["phase_name"])]
    evaluate!(zero(q), cps, q)
end



@get "/evaluate_fitted_background/{idx}" function evaluate(req::HTTP.Request, idx::Int)
    (length(server_state.fit_result) < (idx+1)) && error(BoundsError, ": Index $(idx) is out of range") 
    dict = JSON.parse(String(req.body))
    q = convert.(Float64, dict["q"])
    evaluate!(zero(q), server_state.fit_result[idx+1].background, q)
end

@get "/evaluate/{idx}" function evaluate(req::HTTP.Request, idx::Int)
    dict = JSON.parse(String(req.body))
    q = convert.(Float64, dict["q"])
    evaluate!(zero(q), server_state.phases[idx+1], q)
end


@get "/evaluate_by_name/{phase_name}" function evaluate(req::HTTP.Request, phase_name::String)
    dict = JSON.parse(String(req.body))
    q = convert.(Float64, dict["q"])
    cps = server_state.phases[get_cs_with_names(server_state.phases, [phase_name])]
    evaluate!(zero(q), cps, q)
end


@get "/evaluate_fitted/{idx}" function evaluate_fitted(req::HTTP.Request, idx::Int)
    (length(server_state.fit_result) < (idx+1)) && error(BoundsError, ": Index $(idx) is out of range") 
    dict = JSON.parse(String(req.body))
    q = convert.(Float64, dict["q"])

    pm = server_state.fit_result[idx+1]
    return_dict = Dict{String, Any}()
    for i in eachindex(pm.CPs)
        return_dict[pm.CPs[i].name] = evaluate!(zero(q), pm.CPs[i], q)
    end
    return_dict["background"] = evaluate!(zero(q), pm.background, q)

    return_dict
end


@get "/evaluate_fitted" function evaluate_fitted(req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    (length(server_state.fit_result) < (dict["idx"]+1)) && error(BoundsError, ": Index $(dict["idx"]) is out of range") 
    q = convert.(Float64, dict["q"])
    evaluate!(zero(q), server_state.fit_result[dict["idx"]+1], q)
end

@get "/get_fitted_phase_name/{idx}" function (req::HTTP.Request, idx::Int)
    (length(server_state.fit_result) < (idx+1)) && error(BoundsError, ": Index $(idx) is out of range") 
    return [cp.name for cp in server_state.fit_result[idx+1]]
end



serve()
