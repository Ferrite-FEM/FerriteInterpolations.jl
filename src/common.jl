# Shared helpers. Keep this file minimal: basis functions, dof index tuples and
# reference coordinates always live in the individual element files.

@noinline function throw_out_of_range(ip, i::Int)
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end
