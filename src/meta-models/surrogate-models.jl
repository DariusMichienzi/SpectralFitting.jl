"""
    SurrogateSpectralModel <: AbstractSpectralModel
    SurrogateSpectralModel(modelkind, surrogate, params, params_symbols)

Used to wrap a surrogate function into an [`AbstractSpectralModel`](@ref).

# Example

Creating a surrogate function using [`make_surrogate_harness`](@ref):
```julia
# build and optimize a surrogate model
surrogate = make_surrogate_harness(model, lower_bounds, upper_bounds)

# create surrogate spectral model
sm = SurrogateSpectralModel(
    Multiplicative(),
    surrogate,
    (FitParam(1.0),),
    (:ηH,)
)
```

The `lower_bounds` and `upper_bounds` must be tuples in the form `(E, params...)`, where `E`
denotes the bounds on the energy range to train over.
"""
struct SurrogateSpectralModel{T,K,N,S,Symbols} <: AbstractSpectralModel{T,K}
    surrogate::S
    params::NTuple{N,T}
end

Base.copy(s::SurrogateSpectralModel) = typeof(s)(deepcopy(s.surrogate), deepcopy(s.params))

function SurrogateSpectralModel(
    ::K,
    surrogate::S,
    params::Tuple,
    symbols,
) where {K<:AbstractSpectralModelKind,S}
    P = typeof(params)
    N = length(P.parameters)
    T = first(P.parameters)
    SurrogateSpectralModel{T,K,N,S,symbols}(surrogate, params)
end

unpack_parameters_as_named_tuple(
    model::SurrogateSpectralModel{T,K,N,S,Symbols},
) where {T,K,N,S,Symbols} = NamedTuple{Symbols}(model.params)
remake_with_number_type(
    model::SurrogateSpectralModel{<:FitParam{T},K,N,S,Symbols},
) where {T,K,N,S,Symbols} = SurrogateSpectralModel{T,K,N,S,Symbols}(
    model.surrogate,
    ((get_value(f) for f in model.params)...,),
)

@fastmath function invoke!(output, domain, model::SurrogateSpectralModel{T}) where {T}
    # TODO: rebin if different domains
    output .= model.surrogate(model.params)
end

"""
    optimize_accuracy!(
        surr::AbstractSurrogate,
        obj::Function,
        lb,
        ub;
        sample_type::SamplingAlgorithm = SobolSample(),
        maxiters = 200,
        N_truth = 5000,
        verbose = false,
    )

Improve accuracy (faithfullness) of the surrogate model in recreating the objective function.

Samples a new space of `N_truth` points between `lb` and `ub`, and calculates the objective
function `obj` at each. Finds the point with largest MSE between surrogate and objective, and
adds the point to the surrogate pool. Repeats `maxiters` times, adding `maxiters` points to
surrogate model.

Optionally print to stdout the MSE and iteration count with `verbose = true`.

Note that upper- and lower-bounds should be in the form `(E, params...)`, where `E` is a single
energy and `params` are the model parameters.
"""
function optimize_accuracy!(
    surr::AbstractSurrogate,
    obj::Function,
    lb,
    ub;
    sample_type::SamplingAlgorithm = SobolSample(),
    maxiters = 50,
    verbose = false,
)
    samples = _initial_space(obj, lb, ub, sample_type, maxiters * 10)
    X = samples.x
    y = samples.y
    for epoch = 1:maxiters
        ŷ = surr.(X)
        σ = map(eachindex(y)) do i
            surrogate_error(y[i], ŷ[i])
        end
        ℳσ, i = findmax(σ)
        if verbose
            println("$(rpad(epoch, 5)): ", ℳσ)
        end
        Surrogates.update!(surr, X[i], y[i])
    end
    surr
end

function surrogate_error(y, ŷ)
    sum(j -> (y[j] - ŷ[j])^2, eachindex(y))
end

"""
    wrap_model_as_objective(model::AbstractSpectralModel; ΔE = 1e-1)
    wrap_model_as_objective(M::Type{<:AbstractSpectralModel}; ΔE = 1e-1)

Wrap a spectral model into an objective function for building/optimizing a surrogate model.
Returns an anonymous function taking the tuple `(E, params...)` as the argument, and
returning a single flux value.
"""
function wrap_model_as_objective(model::AbstractSpectralModel, domain)
    outputs = allocate_model_output(model, domain)
    function _surr_objective(params)
        invokemodel!(outputs, domain, model, [params...]) |> copy
    end
end

function _initial_space(obj, lb::NTuple{N}, ub::NTuple{N}, sample_type, n) where {N}
    xs = sample(n, lb, ub, sample_type)
    _xs = if N === 1
        map(i -> Tuple(i), xs)
    else
        xs
    end
    (; x = _xs, y = obj.(_xs))
end

struct SurrogateHarness{M,O,S,L}
    model::M
    obj::O
    surrogate::S
    lower_bounds::L
    upper_bounds::L
end

function optimize_accuracy!(harness::SurrogateHarness, ; kwargs...)
    optimize_accuracy!(
        harness.surrogate,
        harness.obj,
        harness.lower_bounds,
        harness.upper_bounds;
        kwargs...,
    )
end

"""
    make_surrogate_harness(
        model::M,
        lowerbounds::T,
        upperbounds::T;
        optimization_samples = 200,
        seed_samples = 50,
        S::Type = RadialBasis,
        sample_type = SobolSample(),
        verbose = false,
    )

Creates and optimizes a surrogate model of type `S` for `model`, using [`wrap_model_as_objective`](@ref)
and  [`optimize_accuracy!`](@ref) for `optimization_samples` iterations. Model is initially
seeded with `seed_samples` points prior to optimization.

!!! warning
    Additive models integrate energies to calculate flux, which surrogate models are currently not
    capable of. Results for Additive models likely to be inaccurate. This will be patched in a future
    version.
"""
function make_surrogate_harness(
    make_surrogate::Function,
    domain,
    model::M,
    lowerbounds::T,
    upperbounds::T,
    ;
    seed_samples = 50,
    sample_type = SobolSample(),
) where {M<:AbstractSpectralModel,T<:NTuple}
    # do tests here to make sure lower and upper bounds are okay
    obj = wrap_model_as_objective(model, domain)

    K = modelkind(M)
    if K === Additive()
        @warn "Additive models integrate energies to calculate flux, which spectral surrogate models currently do not do. Surrogate is consequently likely to have poor accuracy for Additive models."
    end

    # build surrogate
    inits = _initial_space(obj, lowerbounds, upperbounds, sample_type, seed_samples)

    SurrogateHarness(model, obj, make_surrogate(inits.x, inits.y), lowerbounds, upperbounds)
end

function make_model(harness::SurrogateHarness)
    params = unpack_parameters_as_named_tuple(harness.model)
    SurrogateSpectralModel(
        modelkind(harness.model),
        harness.surrogate,
        Tuple(params),
        keys(params),
    )
end

export SurrogateSpectralModel,
    optimize_accuracy!, wrap_model_as_objective, make_surrogate_harness, make_model
