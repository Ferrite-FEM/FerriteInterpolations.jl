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
using ParallelTestRunner

const TESTDIR = @__DIR__

# `find_tests` auto-discovers every `.jl` file in `test/` (one test file per
# element). Each file runs in its own isolated worker process, so files must be
# self-contained: they carry their own `using` and `include("test_utils.jl")`.
testsuite = find_tests(TESTDIR)

# Shared helpers, `include`d by the tests that need them:
delete!(testsuite, "test_utils")

# Auto CPU thread count detection in ParallelTestRunner is bad
push!(ARGS, "--jobs=$(Sys.CPU_THREADS)")

runtests(
    FerriteInterpolations, ARGS;
    testsuite,
    init_code = :(using FerriteInterpolations, Ferrite),
)
