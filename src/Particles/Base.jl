"""
Abstract supertype for all particle types in CMInject.

Subtypes `MyParticle <: AbstractParticle` should have (at least) two fields:

- `_u`: The phase-space position, a *named* structure (LabelledArray / NamedTuple / struct /...)
- `_p`: A *named* structure containing all particle properties

Subtypes should also implement the following methods:

- mass(particle::MyParticle)
- associated_u_type(::Type{MyParticle})
- associated_props_type(::Type{MyParticle})
- associated_particle_type(::Type{PropsType}),
  where `PropsType` is the type of the `_p` field and
  `associated_particle_type(associated_props_type(MyParticle)) == MyParticle` must be ensured
- transfer_velocities!(du, particle::MyParticle)

A comparatively convenient way to define new particle types is using the [`@declare_particle`](@ref) macro: it
will define all of the above fields and methods automatically based on a declarative syntax.
"""
abstract type AbstractParticle end

## Logic for passing through particle property accesses to the nested fields _u and _p

# (Helper functions for `get_at`)
# copied from LabelledArrays.jl, since it doesn't export this function
@inline _symnames(::Type{SLArray{S,T,N,L,Syms}}) where {S,T,N,L,Syms} = Syms
# helper for getting the type of the fields of a struct, from the struct type
@inline _get_fieldtype(t::Type, f::Symbol) = t.types[findfirst(n->n==f, fieldnames(t))]

"""
    get_at(x::P,::Val{s}) where {s, P<:AbstractParticle}

Given a particle of some specific particle type and a `Val(s)` indicating a property of the particle,
e.g., Val(:x) and a SphericalParticle2D instance, returns this property of the particle.

This is a generated function, meaning to generate the most efficient code possible for each specific field. It's used
to define a specific `Base.getproperty` for all particle types. See also the blog post by Chris Rackauckas on
zero-overhead abstractions, which this is heavily based on:

http://www.stochasticlifestyle.com/zero-cost-abstractions-in-julia-indexing-vectors-by-name-with-labelledarrays/
"""
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

"""
    Base.getproperty(particle::AbstractParticle, s::Symbol)

Returns the property `s` of the particle `particle`, where `s` may be part of the phase-space properties or the
particle properties. This implementation allows you to write code like

    particle.x * particle.TSubtypes `MyField <: Field` should implement:

rather than

    particle._u.x * particle._p.T
"""
Base.@propagate_inbounds function Base.getproperty(particle::AbstractParticle, s::Symbol)
    get_at(particle,Val(s))
end


## Functions for implementing the main loop

"""
    associated_particle_type(props_type)

Must return the associated particle type given the particle props type that's nested within this particle type.
When using the @declare_particle macro, a implementation for each new particle type will automatically be created.
"""
associated_particle_type(props_type) = error("Don't know the particle type of $(props_type)")

"""
    associated_props_type(particle_type)

Must return the associated particle props type given the particle type that's nesting this props type.
When using the @declare_particle macro, a implementation for each new particle type will automatically be created.
"""
associated_props_type(particle_type) = error("Don't know the particle props type of $(particle_type)")

"""
    associated_u_type(particle_type)

Must return the associated particle phase-space type given the particle type that's nesting this props type.
When using the @declare_particle macro, a implementation for each new particle type will automatically be created.
"""
associated_u_type(particle_type) = error("Don't know the phase-space type of $(particle_type)")

"""
    transfer_velocities(du, particle::AbstractParticle)

Enforces second-order dynamics for each particle type as declared by `@velocities` in the `@declare_particle` macro, by
writing each current velocity (e.g., `particle.vx`) to the differential of each associated position (e.g., `du.x`).
"""
transfer_velocities!(du, particle::AbstractParticle) = nothing


"""
    accelerate!(du, acc)

A generated function that applies the acceleration `acc` returned from a field to the differential `du`, in an
efficient way that specializes on the types of `acc` and `du`.

Unless you change the internals of how CMInject simulations are run, you most likely won't need to call this method.
"""
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

"""
    mass(particle::AbstractParticle)

Calculates the mass of the particle. The default implementation returns

    particle.ρ * (4/3)particle.r^3 * π

and may be overridden for a concrete particle type via multiple dispatch.
"""
function mass(particle::P) where P <: AbstractParticle
    particle.ρ * (4/3)particle.r^3 * π
end


## @declare_particle macro definition

# A pseudo-macro to be used within @declare_particle
# We may want to do more metaprogramming in the future, so keep it a macro
macro phase_space(type, names)
    type, names
end

# A pseudo-macro to be used within @declare_particle
# We may want to do more metaprogramming in the future, so keep it a macro
macro velocities(mapping)
    mapping
end

# Finds a macro-call with a particular name within an Expr
function _find_macrocall(args, macro_name)
    idx = findfirst(
        (expr) -> typeof(expr) == Expr && expr.head == :macrocall && expr.args[1] == Symbol(macro_name),
        args
    )
    isnothing(idx) ? nothing : args[idx]
end

# Extracts only the typevar name from a type declaration expression,
# e.g., when called on :(A <: AbstractArray), will return just the symbol :A.
# When called with a Symbol as `expr`, will return exactly that symbol (no-op).
function _extract_typename(expr::Expr)
    if expr.head == :<:
        expr.args[1]
    else
        expr
    end
end
_extract_typename(expr::Symbol) = expr


# TODO require some additional constraints about field names, e.g. forbid :t (reserved for time, relevant for storing)
"""
A macro to conveniently define a new particle struct and let all required methods for a particle type be written
automatically for you.

# Example

    @declare_particle struct MyParticle
        @phase_space Float64 (:x, :y, :z, :vx, :vy, :vz)
        @velocities (:x => :vx, :y => :vy, :z => :vz)
        r::Float64
        ρ::Float64
        T::Float64
    end

Here, `@phase_space Float64 (:x, :y, ...)` declares that the phase-space position of this particle type consists of
the named components `x`, `y` and so on, and that these components have the type `Float64`.
`@velocities (:x => :vx, ...)` declares that `vx` is the velocity for `x`, and so on.
The declarations for `r`, `ρ` and `T` are just as you would write in a regular struct, and will all become part
of the particle properties.
"""
macro declare_particle(expr::Expr)
    # TODO perhaps force all type vars to be unique through gensym, so we won't get clashes with fieldnames
    # (expected example for error: T as fieldname for temperature, but also as the main type variable {T})
    @assert (expr.head == :struct) "Particle definition needs to be a struct, check the docs"
    @assert (length(expr.args) == 3) "Malformed particle struct definition, check the docs"
    _, typedef, block = expr.args
    particle_typename = isa(typedef, Symbol) ? typedef : typedef.args[1]
    typeargs = isa(typedef, Symbol) ? () : typedef.args[2:end]
    typeargs_names = _extract_typename.(typeargs)

    phase_declaration = _find_macrocall(block.args, "@phase_space")
    @assert !isnothing(phase_declaration) "@phase_space declaration is required!"
    phase_typevar, phase_fieldnames_tuple = eval(phase_declaration)
    phase_fieldnames_tuple = eval(phase_fieldnames_tuple)
    phase_fieldnames = collect(phase_fieldnames_tuple)

    # Parse and process the @velocities declaration
    velocity_declaration = _find_macrocall(block.args, "@velocities")
    @assert !isnothing(velocity_declaration) "@velocities declaration is required!"
    velocity_map = eval(velocity_declaration)
    velocity_map = velocity_map isa Tuple ? velocity_map : (velocity_map,)
    @assert all([
        lhs ∈ phase_fieldnames && rhs ∈ phase_fieldnames
        for (lhs, rhs) in velocity_map
    ]) "@velocity declaration must refer to phase-space properties only!"

    # Parse and process the struct's property declarations, which will later actually be put in a "sub-struct"
    prop_declarations = [expr for expr in block.args if typeof(expr) == Expr && expr.head == Symbol("::")]
    prop_fieldnames = [expr.args[1] for expr in prop_declarations]
    clashing_names = intersect(Set(phase_fieldnames), Set(prop_fieldnames))
    @assert isempty(clashing_names) "Property names must be unique; I found these duplicated between phase-space and other props: $(collect(clashing_names))"
    all_fieldnames = vcat(phase_fieldnames, prop_fieldnames)

    # For the generated function force-realization block in the type constructor declaration:
    all_prop_accesses = [:(particle.$prop) for prop in vcat([:_u, :_p], all_fieldnames)]
    # For defining `transfer_velocities!`:
    velocity_assignments = [:(du.$lhs = particle.$rhs) for (lhs, rhs) in velocity_map]

    # Generate new unique type names for the 'inner' phase-space and props types
    phase_typename = gensym("ParticlePhaseSpace")
    props_typename = gensym("ParticleProps")
    phase_type = :($phase_typename{$phase_typevar})
    props_type = typeargs === () ? :($props_typename) : :($props_typename{$(typeargs_names...)})
    particle_type = typeargs === () ? :($particle_typename) : :($particle_typename{$(typeargs_names...)})

    # TODO: move the esc inwards as much as possible, ensuring macro hygiene
    esc(quote
        # Declare the phase-space type
        $phase_typename{$phase_typevar} = CMInject.@SLVector $phase_typevar $phase_fieldnames_tuple

        # Declare the props type (as a mutable struct)
        mutable struct $props_typename{$(typeargs...)}
            $(prop_declarations...)
        end

        # Declare the particle type
        struct $particle_typename{$(typeargs...)} <: CMInject.AbstractParticle  # TODO allow a custom supertype!
            _u::$phase_type  # TODO carry over type variables from outer-level declaration!!
            _p::$props_type

            $particle_type(u, p) where {$(typeargs...)} = let particle = new{$(typeargs...)}(u, p)
                # Force-realize the generated functions for access to all properties.
                # If we don't do this, using solve() will compile to slow versions (for whatever reason), and those
                # will stick. Note that these field accesses are compiled out: The compiler understands that they
                # don't do anything during runtime.
                $(all_prop_accesses...)
                particle
            end

            $particle_type($(phase_fieldnames...), $(prop_fieldnames...)) where {$(typeargs...)} = $particle_type(
                $phase_type($(phase_fieldnames...)),
                $props_type($(prop_fieldnames...))
            )
            $particle_type(; $(phase_fieldnames...), $(prop_fieldnames...)) where {$(typeargs...)} = $particle_type(
                $(phase_fieldnames...),
                $(prop_fieldnames...)
            )
        end

        function Base.show(io::IO, ::MIME"text/plain", particle::$particle_type) where {$(typeargs...)}
            field_list = join([
                "$(sym)="*repr(MIME("text/plain"), getproperty(particle, sym))
                for sym in $all_fieldnames
            ], ", ")
            print(io, "$(typeof(particle))($(field_list))")
        end

        # The block below attaches concrete methods to the utility functions
        # `associated_particle_type` and `carry_velocities!`.

        # Maps the props type (which f! and g! receive) to the particle type
        # (which the acceleration functions should receive)
        function CMInject.associated_particle_type(::Type{$props_type}) where {$(typeargs...)}
            $particle_type
        end

        # Maps the particle type to the wrapped particle props type
        function CMInject.associated_props_type(::Type{$particle_type}) where {$(typeargs...)}
            $props_type
        end

        # Maps the particle type to the wrapped particle props type
        function CMInject.associated_u_type(::Type{$particle_type}) where {$(typeargs...)}
            $phase_typename{$phase_typevar}
        end

        # Based on the @velocities declaration, defines a mapping of position->velocity properties of the phase-space
        # position. Can be used to realize any second-order dynamics through pairs of phase-space terms.
        @inline function CMInject.transfer_velocities!(du, particle::$particle_type) where {$(typeargs...)}
            @inbounds begin
                $(velocity_assignments...)
            end
            nothing
        end

        nothing
    end)
end