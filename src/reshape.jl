function Base.reshape(x::Number, dims...)
    prod(dims) == 1 || throw(DimensionMismatch("new dimensions $(dims) must be consistent with array size 1"))
    return fill(x, dims...)
end
