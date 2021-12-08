# [Developer guide](@id develoeprguide)


## Setup

You should set up your machine the same way as described in the [User guide](userguide.md).

!!! note
    You can skip the Sysimage generation part until you've developed new features, but after making a new release, refer to [Sysimage generation for performant simulations by end-users](@ref) and then run this step yourself to verify that it works as intended.


## Sysimage generation for performant simulations by end-users

For the purpose of allowing users quick execution of simulations without the overhead of JIT compilation, we provide users with a way to generate sysimages (see [the docs from PackageCompiler.jl](https://julialang.github.io/PackageCompiler.jl/v2.0/sysimages.html) about this topic). The scripts to do this are located in `build/create_sysimage.sh` (for graphical usage, like on users' laptops) and `build/create_sysimage_noplot.sh` (for non-graphical usage, like running on compute clusters).

Generating a sysimage that *actually* prevents JIT compilation for all things a user might want to do, depends on us developers providing an example execution file. `PackageCompiler.jl` executes this file and traces which methods were called, to know which concrete methods should be compiled and stored in the sysimage (see [the PackageCompiler.jl docs](https://julialang.github.io/PackageCompiler.jl/v2.0/) for a more detailed explanation).

Therefore, if you add new features (field types, particle types, user options, plotting commands, ... almost anything really), and make a new release of the software, make sure that you add example usage of these features in this example execution file, and enforce that users re-generate this sysimage.

!!! danger
    At CFEL, we store a commonly used sysimage on `/gpfs/cfel/group/cmi/labs/CMInject.jl/build/CMInject.so`. This sysimage *must* be updated when a new release is to be provided to users on Maxwell. Otherwise, all new functionalities you've release will be JIT compiled for every simulation, which will likely be terribly slow.

!!! todo
    Automate this process for users, e.g., through git hooks.


## Code documentation

Please see the APIdocs for documentation on the functions exported from `CMInject`. For guidance with writing documentation, check [the Documentation section in the Julia manual](https://docs.julialang.org/en/v1/manual/documentation/).

For generating documentation, we follow the established directory structure and generation procedure from [Documenter.jl](https://juliadocs.github.io/Documenter.jl/v0.27/man/guide/):

- Documentation Markdown files ("pages") are in `docs/src`
- The documentation is generated by running `julia docs/make.jl`
- The output will be written to `docs/build` - open `docs/build/index.html` as an entry point