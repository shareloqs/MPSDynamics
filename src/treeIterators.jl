struct Radial{T}
    iter::T
    start::Int
end
Radial(iter) = Radial(iter, findheadnode(iter))
struct Walk{T}
    iter::T
end
struct Traverse{T}
    iter::T
end
struct Path <: AbstractArray{Tuple, 1}#integrate with directional
    iter::TreeNetwork
    firstsite::Int
    lastsite::Int
end
Path(iter, lastsite::Int) = Path(iter, findheadnode(iter), lastsite)

struct Directional{T<:Union{Traverse{Tree}, Walk{Tree}}}
    iter::T
end

function Base.iterate(iter::Directional)
    (id1, state1) = iterate(iter.iter)
    next = iterate(iter.iter, state1)
    if next == nothing
        return ((nothing, id1), next)
    else
        id2 = next[1]
        node = iter.iter.iter[id1]
        dir = isconnected(node, id2) ? findbond(node, id2) : 0
        return ((dir, id1), next)
    end
end
function Base.iterate(iter::Directional, current)
    if current == nothing
        return nothing
    else
        id1 = current[1]
        next = iterate(iter.iter, current[2])
        if next == nothing
            return ((nothing, id1), next)
        else
            id2 = next[1]
            node = iter.iter.iter[id1]
            dir = isconnected(node, id2) ? findbond(node, id2) : 0
            return ((dir, id1), next)
        end
    end
end

function Base.iterate(iter::Path)
    pth = path(iter.iter, iter.firstsite, iter.lastsite)
    dir = findbond(iter.iter.tree[iter.firstsite], pth[2])
    return ((iter.iter[iter.firstsite], iter.firstsite, dir), (pth, 2))
end
function Base.iterate(iter::Path, state)
    (pth, i) = state
    if i > length(pth)
        return nothing
    else
        dir = i+1 > length(pth) ? nothing : findbond(iter.iter.tree[pth[i]], pth[i+1])
        return ((iter.iter[pth[i]], pth[i], dir), (pth, i+1))
    end
end

Base.length(iter::Path) = length(path(iter.iter, iter.firstsite, iter.lastsite))
Base.size(iter::Path) = (length(iter),)
Base.IndexStyle(::Type{<:Path}) = IndexLinear()
function Base.getindex(iter::Path, i::Int)
    if i != length(iter)
        ids = path(iter.iter, iter.firstsite, iter.lastsite)[i:i+1]
        dir = findbond(iter.iter.tree[ids[1]], ids[2])
        return (iter.iter[ids[1]], ids[1], dir)
    else
        id = path(iter.iter, iter.firstsite, iter.lastsite)[end]
        return (iter.iter[id], id, nothing)
    end
end

function Base.getindex(iter::Path, i::UnitRange{Int64})
    length(collect(i)) == 0 && return Any[]
    1 <= i.stop <= length(iter) || throw(BoundsError(iter, i.stop))
    1 <= i.start <= length(iter) || throw(BoundsError(iter, i.start))
    
    pth = path(iter.iter, iter.firstsite, iter.lastsite)
    len = i.stop - i.start + 1
    res = Vector{Any}(undef, len)
    for j in 1:len-1
        n = j + i.start - 1
        dir = findbond(iter.iter.tree[pth[n]], pth[n+1])
        res[j] = (iter.iter[pth[n]], pth[n], dir)
    end
    if i.stop == length(pth)
        res[end] = (iter.iter[pth[end]], pth[end], nothing)
    else
        n = len + i.start - 1
        dir = findbond(iter.iter.tree[pth[n]], pth[n+1])
        res[end] = (iter.iter[pth[n]], pth[n], dir)
    end
    return res
end

function Base.getindex(iter::Union{Tree,TreeNetwork}, I::Tuple{Int,Int})
    1 <= I[2] <= length(iter) || throw(BoundsError(iter, I[2]))
    1 <= I[1] <= length(iter) || throw(BoundsError(iter, I[1]))
    tmp = setheadnode(iter.tree, I[2])
    return pathtohead(tmp, I[1])
end


Base.firstindex(iter::Union{Tree,TreeNetwork}) = 1
Base.lastindex(iter::Union{Tree,TreeNetwork}) = length(iter)
Base.eltype(iter::Union{Walk{Tree},Traverse{Tree}}) = Int
Base.eltype(iter::Radial{Tree}) = Vector{Int}
Base.eltype(iter::Union{Walk{TreeNetwork{T}},Traverse{TreeNetwork{T}}}) where {T} = Tuple{Int,T}
Base.eltype(iter::Radial{TreeNetwork{T}}) where {T} = Vector{Tuple{Int,T}}

Base.length(iter::Traverse) = length(iter.iter)
function Base.length(iter::Union{Walk,Radial})
    len=0
    for i in iter
        len += 1
    end
    return len
end
function Base.getindex(iter::Union{Radial,Walk,Traverse}, I::Int)
    1 <= I <= length(iter) || throw(BoundsError(iter, I))
    for (i, val) in enumerate(iter)
        i==I && return val
    end
end
function Base.getindex(iter::Union{Radial,Walk,Traverse}, I::UnitRange{Int})
    1 <= I.start <= length(iter) || throw(BoundsError(iter, I))
    1 <= I.stop <= length(iter) || throw(BoundsError(iter, I))
    len = I.stop-I.start+1
    r = Vector{eltype(iter)}(undef, len)
    i=1
    (val, state) = iterate(iter)
    while i < I.start
        (val, state) = iterate(iter, state)
        i+=1
    end
    r[i-I.start+1] = val
    while i < I.stop
        i+=1
        (val, state) = iterate(iter, state)
        r[i-I.start+1] = val
    end
    return r
end
Base.lastindex(iter::Union{Radial,Walk,Traverse,Path}) = length(iter)
Base.firstindex(iter::Union{Radial,Walk,Traverse,Path}) = 1

function Base.iterate(iter::Radial{Tree}, state::Vector{Int})
    children = state
    length(children) == 0 && return nothing
    tree = iter.iter
    grandchildren = cat(map(x->x.children, tree[children])..., dims=1)
    return (deepcopy(children), grandchildren)
end
function Base.iterate(iter::Radial{Tree})
    tree = iter.iter
    hn = iter.start
    children = tree[hn].children
    return ([hn], children)
end

function Base.iterate(iter::Radial{TreeNetwork{T}}, state::Vector{Int}) where {T}
    net = iter.iter
    start = iter.start
    next = iterate(Radial(net.tree, start), state)
    if next==nothing
        return nothing
    else
        (i, state) = next
        return ([(i[k],net[i[k]]) for k=1:length(i)], state)
    end
end
function Base.iterate(iter::Radial{TreeNetwork{T}}) where {T}
    net = iter.iter
    start = iter.start
    (i, state) = iterate(Radial(net.tree, start))
    return ([(i[k],net[i[k]]) for k=1:length(i)], state)
end

function Base.iterate(iter::Walk{Tree}, state::Tuple{Tree,Int})
    tree = state[1]
    id = state[2]
    children = tree[id].children
    par = id
    if length(children) == 0
        par = tree[par].parent
        par == 0 && return nothing
        popfirst!(tree[par].children)
        return (par, (tree, par))
    end
    child = tree[par].children[1]
    return (child, (tree, child))
end
function Base.iterate(iter::Walk{Tree})
    tree_ = iter.iter
    hn = findheadnode(tree_)
    tree = deepcopy(tree_)
    return (hn, (tree, hn))
end

function Base.iterate(iter::Walk{TreeNetwork{T}}, state::Tuple{Tree,Int}) where {T}
    net = iter.iter
    next = iterate(Walk(net.tree), state)
    if next==nothing
        return nothing
    else
        (i, state) = next
        return ((i, net[i]), state)
    end
end
function Base.iterate(iter::Walk{TreeNetwork{T}}) where {T}
    net = iter.iter
    (i, state) = iterate(Walk(net.tree))
    return ((i, net[i]), state)
end

function Base.iterate(iter::Traverse{Tree}, state::Tuple{Tree,Int})
    tree = state[1]
    id = state[2]
    children = tree[id].children
    par = id
    while length(children) == 0
        par = tree[par].parent
        par == 0 && return nothing
        popfirst!(tree[par].children)
        children = tree[par].children
    end
    child = tree[par].children[1]
    return (child, (tree, child))
end
function Base.iterate(iter::Traverse{Tree})
    tree_ = iter.iter
    hn = findheadnode(tree_)
    tree = deepcopy(tree_)
    return (hn, (tree, hn))
end

function Base.iterate(iter::Traverse{TreeNetwork{T}}, state::Tuple{Tree,Int}) where {T}
    net = iter.iter
    next = iterate(Traverse(net.tree), state)
    if next==nothing
        return nothing
    else
        (i, state) = next
        return ((i, net[i]), state)
    end
end
function Base.iterate(iter::Traverse{TreeNetwork{T}}) where {T}
    net = iter.iter
    (i, state) = iterate(Traverse(net.tree))
    return ((i, net[i]), state)
end

