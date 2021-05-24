function Base.reshape(x::Number, dims...)
    prod(dims) == 1 || throw(DimensionMismatch("new dimensions $(dims) must be consistent with array size 1"))
    return fill(x, dims...)
end

function slice(A::AbstractArray, id::Int)
    s=size(A,id)
    d=ndims(A)
    [A[fill(:,id-1)...,i,fill(:,d-id)...] for i=1:s]
end
