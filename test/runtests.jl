# Pin the CondaPkg environment (Python + symfem, see CondaPkg.toml) to a stable
# path next to this file. Without this it lands in the throwaway `Pkg.test`
# sandbox project and every test run re-provisions it from scratch. Must be set
# before PythonCall loads.
ENV["JULIA_CONDAPKG_ENV"] = joinpath(@__DIR__, ".conda_env")

# Default to plain pip for the pip dependencies: the conda-forge `uv` binary
# that CondaPkg otherwise uses has been observed to crash with SIGILL (e.g. on
# aarch64 Linux VMs). Overridable through the environment.
get!(ENV, "JULIA_CONDAPKG_PIP_BACKEND", "pip")

# CondaPkg discovers CondaPkg.toml in the *active project*, which under
# `Pkg.test` is the throwaway sandbox, not this directory -- copy ours there so
# the symfem pip dependency is actually declared.
let src = joinpath(@__DIR__, "CondaPkg.toml"), dst = joinpath(dirname(Base.active_project()), "CondaPkg.toml")
    src == dst || cp(src, dst; force = true)
end

using FerriteInterpolations
using ParallelTestRunner: find_tests
using Test

const TESTDIR = @__DIR__

# `find_tests` auto-discovers every `.jl` file in `test/` (one test file per
# element). Each file is self-contained (carries its own `using` and
# `include("test_utils.jl")`) and runs in its own module below.
testsuite = find_tests(TESTDIR)

# Shared helpers, `include`d by the tests that need them:
delete!(testsuite, "test_utils")

# Optional positional args select a subset of test files (by prefix), e.g.
# `Pkg.test(test_args = ["test_bernstein"])`.
names = sort!(collect(keys(testsuite)))
if !isempty(ARGS)
    filter!(name -> any(arg -> startswith(name, arg), ARGS), names)
end

# TODO: Run the test files in parallel worker processes again
# (`ParallelTestRunner.runtests`). Currently blocked: PythonCall inside a Malt
# worker leaks every Python temporary (the Julia-side `Py` wrappers are never
# collected, even after `GC.gc(true)`; verified with a minimal repro without
# Test/Ferrite, while the identical loop in a plain process is flat), so the
# symfem cross-checks balloon by hundreds of MB per test file and workers get
# OOM-killed on small machines. Needs upstream investigation (Malt.jl /
# PythonCall.jl). Until then, run each file serially in this process.
@testset "FerriteInterpolations" begin
    for name in names
        @testset "$name" begin
            mod = Module(Symbol(name))
            # `Module(...)` does not auto-define `include` like the `module`
            # keyword does; the test files need it for test_utils.jl.
            Core.eval(mod, :(include(x) = Base.include($mod, x)))
            Base.include(mod, joinpath(TESTDIR, name * ".jl"))
        end
    end
end
