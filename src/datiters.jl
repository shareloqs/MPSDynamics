mutable struct VarT{T}
    m::T
end
mutable struct VarX{T}
    m::T
end

function Base.getindex(iter::VarT{Vector{Array{T,2}}}, x::Int) where T
    numprec = length(iter.m)
    return [iter.m[p][x,:] for p=1:numprec]
end
function Base.getindex(iter::VarT{Vector{Array{T,3}}}, x::Int, y::Int) where T
    numprec = length(iter.m)
    return [iter.m[p][x,y,:] for p=1:numprec]
end
function Base.getindex(iter::VarT{Array{T,2}}, x::Int) where T
    return iter.m[x,:]
end
function Base.getindex(iter::VarT{Array{T,3}}, x::Int, y::Int) where T
    return iter.m[x,y,:]
end
function Base.getindex(iter::VarX{Vector{Array{T,2}}}, t::Int) where T
    numprec = length(iter.m)
    return [iter.m[p][:,t] for p=1:numprec]
end
function Base.getindex(iter::VarX{Vector{Array{T,3}}}, t::Int) where T
    numprec = length(iter.m)
    return [iter.m[p][:,:,t] for p=1:numprec]
end
function Base.getindex(iter::VarX{Array{T,2}}, t::Int) where T
    return iter.m[:,t]
end
function Base.getindex(iter::VarX{Array{T,3}}, t::Int) where T
    return iter.m[:,:,t]
end

Base.length(iter::VarT{Vector{Array{T,2}}}) where T = size((iter.m[1]))[1]
Base.length(iter::VarT{Vector{Array{T,3}}}) where T = prod(size((iter.m[1]))[1:2])
Base.size(iter::VarT{Vector{Array{T,3}}}) where T = size((iter.m[1]))[1:2]
Base.length(iter::VarX{Vector{Array{T,2}}}) where T = size((iter.m[1]))[2]
Base.length(iter::VarX{Vector{Array{T,3}}}) where T = size((iter.m[1]))[3]
Base.length(iter::VarT{Array{T,2}}) where T = size(iter.m)[1]
Base.length(iter::VarT{Array{T,3}}) where T = prod(size(iter.m)[1:2])
Base.size(iter::VarT{Array{T,3}}) where T = size(iter.m)[1:2]
Base.length(iter::VarX{Array{T,2}}) where T = size(iter.m)[2]
Base.length(iter::VarX{Array{T,3}}) where T = size(iter.m)[3]
Base.firstindex(iter::VarT) = 1
Base.firstindex(iter::VarT{Array{T,3}}, dim::Int) where T = 1
Base.firstindex(iter::VarX) = 1
Base.lastindex(iter::VarT) = length(iter)
Base.lastindex(iter::VarT{Array{T,3}}, dim::Int) where T = size(iter)[dim]
Base.lastindex(iter::VarX) = length(iter)

Base.iterate(iter::Union{VarT,VarX}) = (iter[1], 2)
function Base.iterate(iter::Union{VarT,VarX}, i)
    if i > length(iter)
        return nothing
    else
        return (iter[i], i+1)
    end
end
