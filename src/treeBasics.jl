mutable struct TreeNode
    parent::Int
    children::Vector{Int}
end

mutable struct Tree
    nodes::Vector{TreeNode}
end
Tree() = Tree([TreeNode(0, Vector{Int}())])
function Tree(len::Int)
    tree = Tree()
    for i=1:len-1
        addchild!(tree, i)
    end
    return tree
end

mutable struct TreeNetwork{T}
    tree::Tree
    sites::Vector{T}
end
TreeNetwork(sites::Vector{T}) where {T} = TreeNetwork{T}(Tree(length(sites)), sites)

Base.length(tree::Tree) = length(tree.nodes)
Base.length(net::TreeNetwork) = length(net.tree.nodes)
Base.iterate(iter::Tree, state) = iterate(iter.nodes, state)
Base.iterate(iter::Tree) = iterate(iter.nodes)
Base.iterate(iter::TreeNetwork, state) = iterate(iter.sites, state)
Base.iterate(iter::TreeNetwork) = iterate(iter.sites)
Base.getindex(iter::Tree, i::Union{Int,Vector{Int},UnitRange}) = getindex(iter.nodes, i)
Base.getindex(iter::TreeNetwork, i::Union{Int,Vector{Int},UnitRange}) = getindex(iter.sites, i)

Base.setindex!(iter::Tree, v::TreeNode, i) = setindex!(iter.nodes, v, i)
Base.setindex!(iter::TreeNetwork{T}, v::T, i) where {T} = setindex!(iter.sites, v, i)

Base.eltype(iter::Tree) = TreeNode
Base.eltype(iter::TreeNetwork{T}) where {T} = T


"""
    addchild!(tree::Tree, id::Int)

Add child to node `id` of `tree`.

"""
function addchild!(tree::Tree, id::Int)
    1 <= id <= length(tree) || throw(BoundsError(tree, id))
    push!(tree.nodes, TreeNode(id, Vector{}()))
    child = length(tree)
    push!(tree.nodes[id].children, child)
    return tree
end

"""
    addchildren!(tree::Tree, id::Int, n::Int)

Add `n` children to node `id` of `tree`.

"""
function addchildren!(tree::Tree, id::Int, n::Int)
    for i=1:n
        addchild!(tree, id)
    end
end
function addchild!(network::TreeNetwork, id::Int, site)
    addchild!(network.tree, id)
    push!(network.sites, site)
    return network
end
function removechild!(tree::Tree, id::Int)
    loopcheck(tree)
    1 <= id <= length(tree) || throw(BoundsError(tree, id))
    parid = tree[id].parent
    parid == 0 && throw(ErrorException("attempt to remove the head-node"))
    for child in tree[id].children
        removechild!(tree, child)
    end
    filter!(x->x!=id, tree[parid].children)
    deleteat!(tree.nodes, id)
    for (i, node) in enumerate(tree)
        if node.parent > id
            tree[i].parent -= 1
        end
        for (j, child) in enumerate(node.children)
            if (child > id)
                tree[i].children[j] -= 1
            end
        end
    end
    return tree
end

#adds tree2 to tree1 at node id of tree1
function addtree(tree1::Tree, id::Int, tree2::Tree)
    t1=deepcopy(tree1)
    t2=deepcopy(tree2)
    len1 = length(t1)
    push!(t1[id].children, len1+1)
    t2[1].parent = id
    t2[1].children .+= len1
    for i=2:length(t2)
        t2[i].parent += len1
        t2[i].children .+= len1
    end
    Tree([t1.nodes..., t2.nodes...])
end
#same as above but modifies tree1
function addtree!(tree1::Tree, id::Int, tree2::Tree)
    t2=deepcopy(tree2)
    len1 = length(tree1)
    push!(tree1[id].children, len1+1)
    t2[1].parent = id
    t2[1].children .+= len1
    for i=2:length(t2)
        t2[i].parent += len1
        t2[i].children .+= len1
    end
    push!(tree1.nodes, t2.nodes...)
    return tree1
end

function addtree(net1::TreeNetwork, id::Int, net2::TreeNetwork)
    TreeNetwork(addtree(net1.tree, id, net2.tree), [net1.sites..., net2.sites...])
end
function addtree!(net1::TreeNetwork, id::Int, net2::TreeNetwork)
    addtree!(net1.tree, id, net2.tree)
    push!(net1.sites, net2.sites...)
    return net1
end

function bonds(tree::Tree)
    N = length(tree)
    h = findheadnode(tree)
    r = []
    for i=1:h-1
        push!(r, (tree[i].parent, i))
    end
    for i=h+1:N
        push!(r, (tree[i].parent, i))
    end
    r
end
bonds(net::TreeNetwork) = bonds(net.tree)

function bondview(tree::Tree)
    N = length(tree)
    mat = zeros(Int, N, N)
    for bond in bonds(tree)
        mat[bond...] = 1
        mat[reverse(bond)...] = 1
    end
    mat
end
bondview(net::TreeNetwork) = bondview(net.tree)

#finds the leaves of a tree, ie the sites which are only connected to one other site
function leaves(tree::Tree)
    r=[]
    for (i, node) in enumerate(tree)
        if length(node.children)==0 || (length(node.children)==1 && node.parent==0)
            push!(r, i)
        end
    end
    r
end
leaves(net::TreeNetwork) = leaves(net.tree)

#return list of node ids starting with id and ending with the head-node such that each element is the parent of the one to its left
function pathtohead(tree::Tree, id::Int)
    pth = [id]
    hn = findheadnode(tree)
    i=id
    while i != hn
        i = tree[i].parent
        push!(pth, i)
    end
    return pth
end
pathtohead(net::TreeNetwork, id::Int) = pathtohead(net.tree, id)

function pathfromhead(tree::Tree, id::Int)
    pth = [id]
    hn = findheadnode(tree)
    i=id
    while i != hn
        i = tree[i].parent
        push!(pth, i)
    end
    return reverse(pth)
end
pathfromhead(net::TreeNetwork, id::Int) = pathfromhead(net.tree, id)

function path(tree, firstsite::Int, lastsite::Int)
    return tree[(firstsite, lastsite)]
end
function path(tree, lastsite::Int)
    firstsite = findheadnode(tree)
    return tree[(firstsite, lastsite)]
end
    
function loopcheck(tree::Tree)
    len = length(tree)
    ids = [findheadnode(tree)]
    for i=1:len
        for child in tree[i].children
            in(child, ids) && throw(ErrorException("loop found in tree!"))
            push!(ids, child)
        end
    end
end
loopcheck(net::TreeNetwork) = loopcheck(net.tree)

function findheadnode(tree::Tree)
    for i in 1:length(tree)
        if tree[i].parent==0
            return i
        end
    end
end
findheadnode(net::TreeNetwork) = findheadnode(net.tree)

"""
    findchild(node::TreeNode, id::Int)

Return integer corresponding to the which number child site `id` is of `node`.

"""
function findchild(node::TreeNode, id::Int)
    return findfirst(x->x==id, node.children)
end
function findchild(children::Vector{Int}, id::Int)
    return findfirst(x->x==id, children)
end
isconnected(node::TreeNode, id::Int) = in(id, node.children) || id==node.parent    
function findbond(node::TreeNode, id::Int)
    if in(id, node.children)
        return findchild(node, id) + 1
    elseif node.parent == id
        return 1
    else
        throw("node $id is not connected to $node")
    end
end

function setheadnode!(tree::Tree, id::Int)
    par = tree[id].parent
    if par != 0
        #set the parent to the head-node
        setheadnode!(tree, par)
        #println("setting parent of $par to $id")
        tree[par].parent=id
        #println("removing $id as child of $(par)")
        filter!(x->x!=id, tree[par].children)
        #println("adding $par as child of $id")
        pushfirst!(tree[id].children, par)
        #println("setting parent of $id to 0")
        tree[id].parent=0
    end
end
function setheadnode!(net::TreeNetwork, id::Int)
    par = net.tree[id].parent
    if par != 0
        #set the parent to the head-node
        setheadnode!(net, par)
        #println("setting parent of $par to $id")
        net.tree[par].parent=id
        #println("removing $id as child of $par")
        childpos=findchild(net.tree[par].children, id)
        filter!(x->x!=id, net.tree[par].children)
        #println("adding $par as child of $id")
        pushfirst!(net.tree[id].children, par)
        #println("setting parent of $id to 0")
        net.tree[id].parent=0

        Apar = net[par]
        dpar = size(Apar)
        #println("remove dummy index of old head-node")
        Apar = reshape(Apar, dpar[2:end])
        nd = ndims(Apar)
        IA = collect(1:nd)
        IC = collect(1:nd)
        deleteat!(IC, childpos)
        pushfirst!(IC, childpos)
        net[par] = tensorcopy(Apar, IA, IC)
        #println("add dummy index to new head-node")
        Ahead = net[id]
        net[id] = reshape(Ahead, 1, size(Ahead)...)        
    end
end
function setheadnode(tree_::Tree, id::Int)
    tree = deepcopy(tree_)
    setheadnode!(tree, id)
    return tree
end
function setheadnode(net_::TreeNetwork, id::Int)
    net = deepcopy(net_)
    setheadnode!(net, id)
    return net
end

import Base: print, println, show
function print(tree::Tree, id::Int)
    loopcheck(tree)
    print(id)
    nochild = length(tree[id].children)
    if nochild == 1
        print("->")
        print(tree, tree[id].children[1])
    elseif nochild > 1
        print("->")
        printstyled("(",color=:yellow)
        for child in tree[id].children[1:end-1]
            print(tree, child)
            printstyled(";", color=:green)
        end
        print(tree, tree[id].children[end])
        printstyled(")", color=:yellow)
    end
end

print(tree::Tree) = print(tree, findheadnode(tree))
println(tree::Tree) = (print(tree); println())
print(net::TreeNetwork) = (println(net.tree); print(net.sites))
println(net::TreeNetwork) = (print(net); println())

show(io::IO, tree::Tree) = println(tree)
show(io::IO, net::TreeNetwork) = println(net)

import Plots.plot
function plot(tree::Tree; ids=true)
    N = length(tree)
    length(Radial(setheadnode(tree, leaves(tree)[1]))) > 100 && error("tree is too large to plot, consider printing instead")
    hn = findheadnode(tree)
    nodecolor=[RGBA(0,0.6,0.8,0) for i=1:N]
    nodecolor[hn] = RGBA(1,0,0,0)
    names = ids ? collect(1:N) : []
    graphplot(bondview(tree), curves=false, nodecolor=nodecolor, root=:top, names=names)
end
function plot(net::TreeNetwork; bonds=true, ids=true, thickness=false)
    tree=net.tree
    bdims = bonddims(net)
    N = length(tree)
    length(Radial(setheadnode(tree, leaves(tree)[1]))) > 100 && error("tree is too large to plot, consider printing instead")
    hn = findheadnode(tree)
    nodecolor=[RGBA(0,0.6,0.8,0) for i=1:N]
    nodecolor[hn] = RGBA(1,0,0,0)
    edgelabel = bonds ? bdims : nothing
    names = ids ? collect(1:N) : []
    edgewidth = thickness ? bdims./(0.2*findmax(bdims)[1]) : (s,d,w) -> 1
    graphplot(bondview(tree), curves=false, nodecolor=nodecolor, root=:top, edgewidth=edgewidth, edgelabel=edgelabel, names=names)
end

"""
    randtree(numnodes::Int, maxdegree::Int)

Construct a random tree with `nummodes` modes and max degree `maxdegree`.

"""
function randtree(numnodes::Int, maxdegree::Int)
    tree = Tree()
    count = 1
    while count < numnodes
        for leaf in leaves(tree)
            if count==1 || leaf!=1
                r = rand(collect(1:maxdegree))
                if count+r >= numnodes
                    r = numnodes-count
                    addchildren!(tree, leaf, r)
                    count = numnodes
                    break
                end
                count += r
                addchildren!(tree, leaf, r)
            end
        end
    end
    return tree
end
