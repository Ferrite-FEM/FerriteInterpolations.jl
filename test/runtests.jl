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
