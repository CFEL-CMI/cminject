using LabelledArrays


## Abstract Particle definition

abstract type AbstractParticle end


# copied from LabelledArrays.jl, since it doesn't export this function
@inline _symnames(::Type{SLArray{S,T,N,L,Syms}}) where {S,T,N,L,Syms} = Syms
# helper for getting the type of the fields of a struct, from the struct type
@inline _get_fieldtype(t::Type, f::Symbol) = t.types[findfirst(n->n==f, fieldnames(t))]

@generated function get_at(x::P,::Val{s}) where {s, P<:AbstractParticle}
    _utype = _get_fieldtype(x, :_u)
    _ptype = _get_fieldtype(x, :_p)
    if s ∈ _symnames(_utype)
        _eltype = eltype(_utype)
        quote
            Base.@_propagate_inbounds_meta
            Base.getproperty(x._u, s)::$_eltype
        end
    elseif s ∈ fieldnames(_ptype)
        _ftype = _get_fieldtype(_ptype, s)
        quote
            Base.@_propagate_inbounds_meta
            Base.getfield(x._p, s)::$_ftype
        end
    else
        quote
            Base.@_propagate_inbounds_meta
            Base.getfield(x, s)
        end
    end
end

Base.@propagate_inbounds function Base.getproperty(particle::P, s::Symbol) where P<:AbstractParticle
    get_at(particle,Val(s))
end

@generated function accelerate!(du, acc)
    # Why is this written the way it is?
    #
    # - We want the compiler to generate a specific sequence of add-assignment (+=) operations,
    #    for a specific combination of the types of `du` and `acc`
    # - We want to infer the properties to add-assign purely based on the fieldnames of the `acc` value,
    #   so that a Field does not need to care about specific particle types, but only needs to expect
    #   some set of phase-space position properties that it will affect
    # - The exact sequence must be known *at compile time* to give us type-stability
    # - There should be no for-loop showing up during runtime!
    #
    # We will therefore, at compile time, generate a list of assignment operations purely from the
    # fieldnames of the acceleration object's type, and then splat (...) this list into the returned
    # full expression that should make up the actual runtime function.
    #
    # The result in generated code is that this version here, for a specific combination
    # of types of du and acc, will generate a function like:
    #
    # function accelerate!(du, acc)
    #   # assuming that fieldnames(typeof(acc)) == (:vx, :vz)
    #   du.vx += acc.vx
    #   du.vz += acc.vz
    # end
    #
    # rather than something like:
    #
    # function accelerate!(du, acc)
    #   for name in (:vx, :vz)
    #     setproperty!(du, name, getproperty(du, name) + getproperty(acc, name))
    #   end
    # end
    #
    # and only the former is type-stable.
    names = fieldnames(acc)
    setter_actions = [:(du.$name += acc.$name) for name in names]
    quote
        Base.@_propagate_inbounds_meta
        $(setter_actions...)
    end
end


## Generic particle function definitions

function mass(particle::P) where P <: AbstractParticle
    particle.ρ * (4/3)particle.r^3 * π
end


## SphericalParticle2D definition

# TODO: auto-generate the following complex declarations via a readable macro DSL,
# e.g. a declaration like
#
# @particle struct SphericalParticle2D{T}
#    @phase_space T (:x, :z, :vx, :vz)
#    @velocities (:x => :vx, :z => :vz)
#    r::T
#    ρ::T
# end


# These serve as no-op/error function declarations, for concrete methods to be attached to
associated_particle_type(x) = error("Unknown particle props type $(typeof(x))")
carry_velocities!(du, u) = nothing


macro phase_space(type, names)
    type, names
end

macro velocities(mapping)
    mapping
end

function find_macrocall(args, macro_name)
    idx = findfirst(
        (expr) -> typeof(expr) == Expr && expr.head == :macrocall && expr.args[1] == Symbol(macro_name),
        args
    )
    isnothing(idx) ? nothing : args[idx]
end

macro declare_particle(expr::Expr)
    # TODO FIXME the declared type variables are not considered explicitly!
    @assert (expr.head == :struct) "Particle definition needs to be a struct, check the docs"
    @assert (length(expr.args) == 3) "Malformed particle struct definition, check the docs"
    _, typedef, block = expr.args
    ptype = typedef.args[1]
    typeargs = typedef.args[2:end]

    phase_declaration = find_macrocall(block.args, "@phase_space")
    @assert !isnothing(phase_declaration) "@phase_space declaration is required!"
    velocity_declaration = find_macrocall(block.args, "@velocities")
    @assert !isnothing(velocity_declaration) "@velocities declaration is required!"

    velocity_map = eval(velocity_declaration)
    velocity_assignments = [:(du.$lhs = u.$rhs) for (lhs, rhs) in velocity_map]

    phase_type = gensym("PhaseSpace")
    props_type = gensym("ParticleProps")
    phase_typevar, phase_names_tuple = eval(phase_declaration)
    phase_names_tuple = eval(phase_names_tuple)
    phase_fieldnames = collect(phase_names_tuple)
    @assert all([
        lhs ∈ phase_fieldnames && rhs ∈ phase_fieldnames
        for (lhs, rhs) in velocity_map
    ]) "@velocity declaration must refer to phase-space position properties only!"

    # Retrieve the struct's property declarations, which will later actually be put in a "sub-struct"
    prop_declarations = [expr for expr in block.args if typeof(expr) == Expr && expr.head == Symbol("::")]
    prop_fieldnames = [expr.args[1] for expr in prop_declarations]
    # Required for the generated function force-realization block in the type constructor declaration
    all_prop_accesses = [:(p.$prop) for prop in vcat([:_u, :_p], phase_fieldnames, prop_fieldnames)]

    esc(quote
        # Declare the phase-space type
        $phase_type{$phase_typevar} = @SLVector $phase_typevar $phase_names_tuple
        # Declare the props type (as a mutable struct)
        mutable struct $props_type{$(typeargs...)}
            $(prop_declarations...)
        end

        struct $ptype{$(typeargs...)} <: CMInject.AbstractParticle
            _u::$phase_type{$phase_typevar}  # TODO carry over type variables from outer-level declaration!!
            _p::$props_type{$(typeargs...)}

            $ptype{$(typeargs...)}(u, p) where {$(typeargs...)} = let p = new{$(typeargs...)}(u, p)
                # Force-realize the generated functions for access to all properties.
                # If we don't do this, using solve() will compile to slow versions (for whatever reason), and those
                # will stick. Note that these field accesses are compiled out: The compiler understands that they
                # don't do anything during runtime.
                $(all_prop_accesses...)
                p
            end
        end
        $ptype{$(typeargs...)}($(phase_fieldnames...), $(prop_fieldnames...)) where {$(typeargs...)} = $ptype{$(typeargs...)}(
            $phase_type{$phase_typevar}($(phase_fieldnames...)),
            $props_type{$(typeargs...)}($(prop_fieldnames...))
        )
        $ptype{$(typeargs...)}(; $(phase_fieldnames...), $(prop_fieldnames...)) where {$(typeargs...)} = $ptype{$(typeargs...)}(
            $(phase_fieldnames...),
            $(prop_fieldnames...)
        )

        # The block below attaches concrete methods to the utility functions
        # `associated_particle_type` and `carry_velocities!`.

        # Maps the props type (which f! and g! receive) to the particle type
        # (which the acceleration functions should receive)
        function CMInject.associated_particle_type(::Type{$props_type{$(typeargs...)}}) where {$(typeargs...)}
            $ptype{$(typeargs...)}
        end
        # Based on the @velocities declaration, defines a mapping of position->velocity properties of the phase-space
        # position. Can be used to realize any second-order dynamics between pairs of phase-space terms.
        @inline function CMInject.carry_velocities!(du, u::$phase_type{$phase_typevar}) where $phase_typevar
            @inbounds begin
                $(velocity_assignments...)
            end
            nothing
        end

        nothing
    end)
end

@declare_particle struct SphericalParticle2D{T}
    @phase_space T (:x, :z, :vx, :vz)
    @velocities (:x => :vx, :z => :vz)
    r::T
    ρ::T
end


## Example instances

const example_particle = SphericalParticle2D{Float64}(r=100e-9, ρ=1050.0, x=1e-4, z=-0.1284, vx=7e-5, vz=3.63)
const example_u = example_particle._u
const example_p = example_particle._p