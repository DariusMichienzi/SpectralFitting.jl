export AbstractSpectralModel,
    AbstractSpectralModelKind,
    Multiplicative,
    Additive,
    Convolutional,
    modelkind,
    AbstractSpectralModelImplementation,
    XSPECImplementation,
    JuliaImplementation,
    implementation,
    AbstractSpectralModelClosureType,
    WithClosures,
    WithoutClosures,
    closurekind,
    has_closure_params,
    invokemodel,
    invokemodel!,
    objective_cache_count,
    modelparameters,
    updateparameters!,
    make_parameter_cache

"""
    abstract type AbstractSpectralModel{T,K}

Supertype of all spectral models. Sub-types must implement the following interface
- [`modelkind`](@ref)
- [`SpectralFitting.invoke!`](@ref)

The parametric type parameter `T` is the number type of the model and `K` defines the [`AbstractSpectralModelKind`](@ref).
"""
abstract type AbstractSpectralModel{T,K} end
numbertype(::AbstractSpectralModel{T}) where {T<:Number} = T
numbertype(::AbstractSpectralModel{FitParam{T}}) where {T<:Number} = T

supports_contiguosly_binned(::Type{<:AbstractSpectralModel}) = true

"""
    abstract type AbstractSpectralModelKind

Abstract type of all model kinds. The algebra of models is as follows
```julia
A + A = A
M * M = M
M * A = A
C(A)  = A
```
where `A` is [`Additive`](@ref), `M` is [`Multiplicative`](@ref), and `C` is [`Convolutional`](@ref).
All other operations are prohibited, e.g. `C(M)` or `M * C`. To obtain `M * C` there must be an
additive component, e.g. `M * C(A)`.
"""
abstract type AbstractSpectralModelKind end
"""
    Additive <: AbstractSpectralModelKind
    Additive()

Additive models are effectively the sources of photons, and are the principle building blocks
of composite models. Every additive model has a normalisation parameter which re-scales the
output by a constant factor `K`.

!!! note
    Defining custom additive models requires special care. See [Defining new models](@ref).
"""
struct Additive <: AbstractSpectralModelKind end
"""
    Multiplicative <: AbstractSpectralModelKind
    Multiplicative()

Multiplicative models act on [`Additive`](@ref) models, by element-wise
multiplying the output in each energy bin of the additive model by a different factor.
"""
struct Multiplicative <: AbstractSpectralModelKind end
"""
    Convolutional <: AbstractSpectralModelKind
    Convolutional()

Convolutional models act on the output generated by [`Additive`](@ref) models, similar to
[`Multiplicative`](@ref) models, however may convolve kernels through the output also.
"""
struct Convolutional <: AbstractSpectralModelKind end

"""
    modelkind(M::Type{<:AbstractSpectralModel})
    modelkind(::AbstractSpectralModel)

Return the kind of model given by `M`: either `Additive`, `Multiplicative`, or `Convolutional`.
"""
modelkind(::Type{<:AbstractSpectralModel{T,K}}) where {T,K} = K()
modelkind(::M) where {M<:AbstractSpectralModel} = modelkind(M)

abstract type AbstractSpectralModelImplementation end
struct XSPECImplementation <: AbstractSpectralModelImplementation end
struct JuliaImplementation <: AbstractSpectralModelImplementation end

implementation(::Type{<:AbstractSpectralModel}) = JuliaImplementation()
implementation(model::AbstractSpectralModel) = implementation(typeof(model))

abstract type AbstractSpectralModelClosureType end
struct WithClosures <: AbstractSpectralModelClosureType end
struct WithoutClosures <: AbstractSpectralModelClosureType end

closurekind(::Type{<:AbstractSpectralModel}) = WithoutClosures()

has_closure_params(::WithClosures) = true
has_closure_params(::WithoutClosures) = false
has_closure_params(M::Type{<:AbstractSpectralModel}) = has_closure_params(closurekind(M))
has_closure_params(::M) where {M<:AbstractSpectralModel} = has_closure_params(M)

Δoutput_length(::Type{<:AbstractSpectralModel}) = -1
Δoutput_length(::M) where {M<:AbstractSpectralModel} = Δoutput_length(M)

# interface for ConstructionBase.jl
function ConstructionBase.setproperties(
    model::M,
    patch::NamedTuple{names},
) where {M<:AbstractSpectralModel,names}
    symbols = all_parameter_symbols(model)
    args = (s in names ? getproperty(patch, s) : getproperty(model, s) for s in symbols)
    M(args...)
end
ConstructionBase.constructorof(::Type{M}) where {M<:AbstractSpectralModel} = M

# implementation interface
# never to be called directly
# favour `invokemodel!` instead
"""
    SpectralFitting.invoke!(output, energy, M::Type{<:AbstractSpectralModel}, params...)

Used to define the behaviour of models. Should calculate output of the model and write in-place
into `output`.

!!! warning
    This function should not be called directly. Use [`invokemodel`](@ref) instead.

Parameters are passed in in-order as defined in the model structure. For example
```julia
struct MyModel{F1,F2,F3,...} <: AbstractSpectralModel
    p1::F1
    p2::F2
    p3::F3
    # ...
end
```
would have the arguments passed to `invoke!` as
```julia
function SpectralFitting.invoke!(output, energy, ::Type{<:MyModel}, p1, p2, p3, ...)
    # ...
end
```

The only exception to this are [`Additive`](@ref) models, where the normalisation parameter
`K` is not passed to `invoke!`.
"""
invoke!(output, energy, M::AbstractSpectralModel) = error("Not defined for $(M).")

"""
    invokemodel(energy, model)
    invokemodel(energy, model, free_params)

Invoke the [`AbstractSpectralModel`](@ref) given by `model`, optionally overriding the free
parameters with values given in `free_params`. `free_params` may be a vector or tuple with element
type [`FitParam`](@ref) or `Number`.

This function, unlike [`SpectralFitting.invoke!`](@ref) used to define models, is sensitive to performing
any normalisation or post-processing tasks that a specific model kind may require.

!!! note
    Users should always call models using [`invokemodel`](@ref) or [`invokemodel!`](@ref) to ensure
    normalisations and closures are accounted for.

`invokemodel` allocates the needed output arrays based on the element type of `free_params` to allow
automatic differentation libraries to calculate parameter gradients.

In-place non-allocating variants are the [`invokemodel!`](@ref) functions.

# Example

```julia
model = XS_PowerLaw()
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, model)

p0 = [0.1, 2.0] # change K and a
invokemodel(energy, model, p0)
```
"""
function invokemodel(e, m::AbstractSpectralModel)
    output = construct_objective_cache(m, e) |> vec
    invokemodel!(output, e, m)
    output
end

"""
    invokemodel!(output, energy, model)
    invokemodel!(output, energy, model, free_params)
    invokemodel!(output, energy, model, free_params, frozen_params)

In-place variant of [`invokemodel`](@ref), calculating the output of an [`AbstractSpectralModel`](@ref)
given by `model`, optionally overriding the free and/or frozen parameter values. These arguments
may be a vector or tuple with element type [`FitParam`](@ref) or `Number`.

The number of fluxes to allocate for a model may change if using any [`CompositeModel`](@ref)
as the `model`. It is generally recommended to use [`objective_cache_count`](@ref) to ensure the correct number
of output arrays are allocated with [`construct_objective_cache`](@ref) when using composite models.

Single spectral model components should use [`make_flux`](@ref) instead.

# Example

```julia
model = XS_PowerLaw()
energy = collect(range(0.1, 20.0, 100))
output = make_flux(model, energy)
invokemodel!(output, energy, model)

p0 = [0.1, 2.0] # change K and a
invokemodel!(output, energy, model, p0)
```
"""

@inline function invokemodel!(f, e, m::AbstractSpectralModel{<:FitParam})
    # need to extract the parameter values
    model = remake_with_number_type(m)
    invokemodel!(view(f, :, 1), e, model)
end
@inline function invokemodel!(f, e, m::AbstractSpectralModel, cache::ParameterCache)
    invokemodel!(f, e, m, cache.parameters)
end
@inline function invokemodel!(f, e, m::AbstractSpectralModel, parameters::AbstractArray)
    invokemodel!(view(f, :, 1), e, remake_with_parameters(m, parameters))
end

invokemodel!(
    f::AbstractVector,
    e::AbstractVector,
    m::AbstractSpectralModel{<:Number,K},
) where {K} = invokemodel!(f, e, K(), m)
@inline function invokemodel!(
    output::AbstractVector,
    domain::AbstractVector,
    ::Additive,
    model::AbstractSpectralModel,
)
    invoke!(output, domain, model)
    # perform additive normalisation
    @. output *= model.K
    output
end
@inline function invokemodel!(
    output::AbstractVector,
    domain::AbstractVector,
    ::AbstractSpectralModelKind,
    model::AbstractSpectralModel,
)
    invoke!(output, domain, model)
    output
end

# printing

function _printinfo(io::IO, m::M) where {M<:AbstractSpectralModel}
    param_tuple = all_parameters_to_named_tuple(m)
    params = [String(s) => p for (s, p) in zip(keys(param_tuple), param_tuple)]
    print(io, "$(FunctionGeneration.model_base_name(M))\n")

    pad = maximum(i -> length(first(i)), params) + 1

    for (s, val) in params
        print(io, "   $(rpad(s, pad)) => ")
        println(io, val)
    end
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(model::AbstractSpectralModel))
    buff = IOBuffer()
    _printinfo(buff, model)
    s = String(take!(buff))
    print(io, encapsulate(s))
end

modelparameters(model::AbstractSpectralModel{T}) where {T} =
    T[model_parameters_tuple(model)...]

# todo: this function could be cleaned up with some generated hackery 
function remake_with_number_type(model::AbstractSpectralModel{P}, T::Type) where {P}
    M = typeof(model).name.wrapper
    params = model_parameters_tuple(model)
    new_params = if P <: FitParam
        convert.(T, get_value.(params))
    else
        convert.(T, param)
    end
    M{T}(new_params...)
end
remake_with_number_type(model::AbstractSpectralModel{FitParam{T}}) where {T} =
    remake_with_number_type(model, T)

"""
    updatemodel(model::AbstractSpectralModel; kwargs...)
    updatemodel(model::AbstractSpectralModel, patch::NamedTuple)

Modify parameters in a given model by keyvalue, or with a named tuple.
"""
updatemodel(model::AbstractSpectralModel, patch::NamedTuple) =
    ConstructionBase.setproperties(model, patch)
updatemodel(model::AbstractSpectralModel; kwargs...) =
    ConstructionBase.setproperties(model; kwargs...)

@inline function updatefree(model::AbstractSpectralModel, free_params)
    patch = free_parameters_to_named_tuple(free_params, model)
    updatemodel(model, patch)
end

@inline function updateparameters(model::AbstractSpectralModel, params)
    patch = all_parameters_to_named_tuple(params, model)
    updatemodel(model, patch)
end

@inline function updateparameters(model::AbstractSpectralModel; params...)
    updatemodel(model; params...)
end

# for modifying models with FitParams
function updateparameters!(model::AbstractSpectralModel{<:FitParam}, params::AbstractVector)
    for (i, s) in enumerate(all_parameter_symbols(model))
        v = params[i]
        if typeof(v) <: FitParam
            set!(getproperty(model, s), v)
        else
            set_value!(getproperty(model, s), get_value(v))
        end
    end
end

function updateparameters!(model::AbstractSpectralModel{<:FitParam}; params...)
    for (s, v) in params
        if typeof(v) <: FitParam
            set!(getproperty(model, s), v)
        else
            set_value!(getproperty(model, s), get_value(v))
        end
    end
    model
end

_allocate_free_parameters(model::AbstractSpectralModel) =
    filter(isfree, modelparameters(model))

function make_parameter_cache(model::AbstractSpectralModel)
    parameters = modelparameters(model)
    ParameterCache(parameters)
end

function make_diff_parameter_cache(
    model::AbstractSpectralModel;
    param_diff_cache_size = nothing,
)
    parameters = modelparameters(model)
    free_mask = _make_free_mask(parameters)

    vals = map(get_value, parameters)
    N = isnothing(param_diff_cache_size) ? length(vals) : param_diff_cache_size
    diffcache = DiffCache(vals, ForwardDiff.pickchunksize(N))

    # embed current parameter values inside of the dual cache
    # else all frozens will be zero
    get_tmp(diffcache, ForwardDiff.Dual(one(eltype(vals)))) .= vals

    ParameterCache(free_mask, diffcache)
end
