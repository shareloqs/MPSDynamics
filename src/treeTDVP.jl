drop=Iterators.drop

mutable struct TreeLightCone <: LightCone
    ref::Vector
    edge::Vector{Int}
    thresh::Float64
    TreeLightCone(ref::TreeNetwork, edge::Vector{Int}, thresh::Float64) = new(orthcentersmps(ref), edge, thresh)
end
TreeLightCone(A0::TreeNetwork, rad::Int=2, thresh::Float64=DEFLCTHRESH) = TreeLightCone(A0, Radial(A0.tree)[rad], thresh)
LightCone(A0::TreeNetwork, rad::Int=2, thresh::Float64=DEFLCTHRESH) = TreeLightCone(A0, rad, thresh)

#returns list of the orthoganality centres of A, assumes A is right-normalised
function orthcentersmps(A::TreeNetwork)
    B = deepcopy(A.sites)
    for (id, X) in drop(Traverse(A), 1)
        nc = length(A.tree[id].children)
        par = A.tree[id].parent
        dir = findbond(A.tree[par], id)
        AL, C = QR(B[par], dir)
        B[id] = contractC!(B[id], C, 1)
    end
    return B
end

function physdims(M::TreeNetwork)
    N = length(M)
    res = Vector{Int}(undef, N)
    for (i, site) in enumerate(M)
        res[i] = size(site)[end]
    end
    return Dims(res)
end

function mpsrightnorm!(net::TreeNetwork, id::Int)
    loopcheck(net)
    children = net.tree[id].children
    nc = length(children)
    for (i, child) in enumerate(children)
        length(net.tree[child].children) >= 1 && mpsrightnorm!(net, child)
        dchild = size(net[child])
        dpar = size(net[id])
        C, AR = lq(reshape(net[child], dchild[1], :))
        net[child] = reshape(Matrix(AR), dchild)
        IC=collect(1:nc+2)
        IA=collect(1:nc+2)
        IC[i+1]=-1

        net[id] = tensorcontract(net[id], IA, C, [i+1,-1], IC)
    end
end

"""
    mpsrightnorm!(A::TreeNetwork)

When applied to a tree-MPS, right normalise towards head-node.

"""
mpsrightnorm!(net::TreeNetwork) = mpsrightnorm!(net, findheadnode(net))

"""
    mpsmixednorm!(A::TreeNetwork, id::Int)

Normalise tree-MPS `A` such that orthogonality centre is on site `id`.

"""
function mpsmixednorm!(net::TreeNetwork, id::Int)
    1 <= id <= length(net) || throw(BoundsError(net.tree, id))
    setheadnode!(net, id)
    mpsrightnorm!(net, id)
end

"""
    mpsmoveoc!(A::TreeNetwork, id::Int)

Move the orthogonality centre of right normalised tree-MPS `A` to site `id`.

This function will be more efficient than using `mpsmixednorm!` if the tree-MPS is already right-normalised.

"""
function mpsmoveoc!(A::TreeNetwork, id::Int)
    for site in drop(pathfromhead(A.tree, id), 1)
        mpsshiftoc!(A, site)
    end
end

"""
    mpsshiftoc!(A::TreeNetwork, newhd::Int)

Shift the orthogonality centre by one site, setting new head-node `newhd`.

"""
function mpsshiftoc!(A::TreeNetwork, newhd::Int)
    oldhd = findheadnode(A)
    in(newhd, A.tree[oldhd].children) || throw("site $id is not child of head-node")

    setheadnode!(A, newhd)
    AL, C = QR(A[oldhd], 1)
    A[oldhd] = AL
    A[newhd] = contractC!(A[newhd], C, findchild(A.tree[newhd], oldhd)+1)
end

function calcbonddims!(tree::Tree, physdims::Dims, Dmax::Int, M::Array{Int, 2}, id::Int)
    for child in tree[id].children
        length(tree[child].children) >= 1 && calcbonddims!(tree, physdims, Dmax, M, child)
        D = physdims[child]
        for grandchild in tree[child].children
            D *= M[child, grandchild]
        end
        M[id, child] = min(D, Dmax)
        M[child, id] = min(D, Dmax)
    end
end
function calcbonddims(tree::Tree, physdims::Dims, Dmax::Int)
    loopcheck(tree)
    N = length(tree)
    M = zeros(Int, N, N)
    calcbonddims!(tree, physdims, Dmax, M, findheadnode(tree))
    return M
end

"""
    randmps(tree::Tree, physdims, Dmax::Int, T::Type{<:Number} = Float64)
Construct a random, right-normalised, tree-MPS, with structure given by tree and max bond-dimension given by `Dmax`.

The local Hilbert space dimensions are specified by physdims which can either be of type `Dims{length(tree)}`, specifying
the dimension of each site, or of type `Int`, in which case the same local dimension is used for every site.

"""#check bond-dims correct
function randmps(tree_::Tree, physdims::Dims, Dmax::Int, T::Type{<:Number} = Float64)
    tree = deepcopy(tree_)
    hn = findheadnode(tree)
    leafnodes = leaves(tree)
    N = length(tree)
    setheadnode!(tree, leafnodes[1])
    bonddims1 = calcbonddims(tree, physdims, Dmax)
    setheadnode!(tree, leafnodes[2])
    bonddims2 = calcbonddims(tree, physdims, Dmax)
    bonddims = min.(bonddims1, bonddims2)
    setheadnode!(tree, hn)
    
    A = Vector{AbstractArray}(undef, N)
    for (id, node) in enumerate(tree)
        if id != hn
            Dpar = bonddims[id, node.parent]
            d = physdims[id]
            Dchildren = bonddims[id, node.children]
            A[id] = reshape(randisometry(T, Dpar, prod(Dchildren)*d), Dpar, Dchildren..., d)
        else
            d = physdims[id]
            Dchildren = bonddims[id, node.children]
            A[id] = reshape(randisometry(T, 1, prod(Dchildren)*d), 1, Dchildren..., d)
        end
    end
    TreeNetwork(tree, A)
end
randmps(tree_::Tree, d::Int, Dmax::Int, T::Type{<:Number} = Float64) = randmps(tree_, ntuple(i -> d, length(tree_)), Dmax, T)

function normmps(net::TreeNetwork, id::Int)
    nc = length(net.tree[id].children)
    IA = collect(1:nc+2)
    IB = collect(nc+3:2*nc+4)
    IA[end] = -1 #contract physical indices
    IB[end] = -1 #contract physical indices 
    ρ = tensorcontract(net[id], IA, conj(net[id]), IB)
    for (i, child) in enumerate(net.tree[id].children)
        ρchild = normmps(net, child)
        nd = (nc+2-i)*2
        halfnd = div(nd,2)
        IA = collect(1:nd)
        IA[2] = -1
        IA[halfnd+2] = -2
        ρ = tensorcontract(ρ, IA, ρchild, [-1, -2])
    end
    return ρ
end

"""
    normmps(net::TreeNetwork; mpsorthog=:None)
    
When applied to a tree-MPS `mpsorthog=:Left` is not defined.

"""
function normmps(net::TreeNetwork; mpsorthog=:None)
    loopcheck(net)
    if mpsorthog==:Right
        OC = findheadnode(net)
        AC = net[OC]
        nd = ndims(AC)
        IA = collect(1:nd)
        return real(scalar(tensorcontract(AC, IA, conj(AC), IA)))
    elseif typeof(mpsorthog)<:Int
        OC = mpsorthog
        AC = net[OC]
        nd = ndims(AC)
        IA = collect(1:nd)
        return real(scalar(tensorcontract(AC, IA, conj(AC), IA)))
    elseif mpsorthog==:None
        hn = findheadnode(net)
        ρ = normmps(net, hn)
        return real(ρ[1])
    end
end

function initenvs!(A::TreeNetwork, M::TreeNetwork, F::Vector, id::Int)
    for child in A.tree[id].children
        F = initenvs!(A, M, F, child)
    end
    F[id] = updaterightenv(A[id], M[id], F[A.tree[id].children]...)
    return F
end
function initenvs(A::TreeNetwork, M::TreeNetwork, F::Nothing)
    hn = findheadnode(A)
    N = length(A)
    F = Vector{Any}(undef, N)
    for child in A.tree[hn].children
        F = initenvs!(A, M, F, child)
    end
    return F
end
function initenvs(A::TreeNetwork, M::TreeNetwork, F::Vector)
    return F
end

function tdvp1sweep_test!(dt, A::TreeNetwork, M::TreeNetwork, F=nothing; verbose=false, kwargs...)
    F = initenvs(A, M, F)

    for (dir, i) in Directional(Walk(A.tree))
        children = A.tree[i].children
        par = A.tree[i].parent
        if par != 0
            Fpar = F[par]
        else
            Fpar = fill!(similar(M[1], (1,1,1)), 1)
        end
            
        A[i], info = exponentiate(x->applyH1(x, M[i], Fpar, F[children]...), -im*dt/2, A[i]; ishermitian=true)

        if verbose
            E = real(dot(A[i], applyH1(A[i], M[i], Fpar, F[children]...)))
            println("Sweep L->R: AC on site $i, energy = $E")
        end

        dir == nothing && break
        
        A[i], C = QR(A[i], dir)
        F[i] = updateenv(A[i], M[i], dir, Fpar, F[children]...)

        legs = [par, children...]
        next = legs[dir]
        C, info = exponentiate(x->applyH0(x, F[i], F[next]), im*dt/2, C; ishermitian=true)
        if verbose
            E = real(dot(C, applyH0(C, F[i], F[next])))
            println("Sweep L->R: C between sites $i and $next, energy = $E")
        end

        idir = findbond(A.tree[next], i)
        A[next] = contractC!(A[next], C, idir)
    end

    return A, F
end

tdvp1sweep!(dt, A::TreeNetwork, M::TreeNetwork, F=nothing; verbose=false, kwargs...) =
    tdvp1sweep!(dt, A, M, initenvs(A, M, F), findheadnode(A); verbose=verbose, kwargs...)
function tdvp1sweep!(dt, A::TreeNetwork, M::TreeNetwork, F::Vector, id::Int; verbose=false, kwargs...)

    children = A.tree[id].children
    parent = A.tree[id].parent
    nc = length(children)
    AC = A[id]

    F0 = parent==0 ? fill!(similar(M[1], (1,1,1)), 1) : F[parent]

    #(OC begins on node)
    #evolve node forward half a time step
    AC, info = exponentiate(x->applyH1(x, M[id], F0, F[children]...), -im*dt/2, AC; ishermitian=true)
    if verbose
        E = real(dot(AC, applyH1(AC, M[id], F0, F[children]...)))
        println("Sweep L->R: AC on site $id, energy = $E")
    end

    for (i, child) in enumerate(children)
        grandchildren = A.tree[child].children
        otherchildren = filter(x->x!=child, children)

        #extract C from node
        AL, C = QR(AC, i+1)
        F[id] = updateleftenv(AL, M[id], i, F0, F[otherchildren]...)

        #evolve C backwards half a time step
        C, info = exponentiate(x->applyH0(x, F[id], F[child]), im*dt/2, C; ishermitian=true)
        if verbose
            E = real(dot(C, applyH0(C, F[id], F[child])))
            println("Sweep L->R: C between sites $id and $child, energy = $E")
        end

        #contract C with child
        A[child] = contractC!(A[child], C, 1)
        #(OC is now on child)
        
        #evolve child forward one full time step
        A, F = tdvp1sweep!(dt, A, M, F, child; verbose=verbose, kwargs...)

        #extract C from child
        AR, C = QR(A[child], 1)
        A[child] = AR
        F[child] = updaterightenv(AR, M[child], F[grandchildren]...)
        
        #evolve C backwards half a time step
        C, info = exponentiate(x->applyH0(x, F[child], F[id]), im*dt/2, C; ishermitian=true)
        if verbose
            E = real(dot(C, applyH0(C, F[child], F[id])))
            println("Sweep R->L: C between sites $child and $id, energy = $E")
        end

        #contract C with node
        AC = contractC!(AL, C, i+1)
        #(OC is now on node)
    end

    #evolve node forward half a time step
    AC, info = exponentiate(x->applyH1(x, M[id], F0, F[children]...), -im*dt/2, AC; ishermitian=true)
    if verbose
        E = real(dot(AC, applyH1(AC, M[id], F0, F[children]...)))
        println("Sweep R->L: AC on site $id, energy = $E")
    end

    A[id] = AC
    #node id has now been evolved one full time step
    return A, F
end

function tdvp1sweep_lc!(dt, A::TreeNetwork, M::TreeNetwork, lc::TreeLightCone, F=nothing; verbose=false, kwargs...)

    hn = findheadnode(A)
    if hn != findheadnode(M)
        setheadnode!(M, hn)
    end

    F = initenvs(A, M, F)
    F0 = fill!(similar(M[1], (1,1,1)), 1)
    children = A.tree[hn].children
    nc = length(children)
    AC = A[hn]

    #(OC begins on headnode)
    #evolve headnode forward half a time step
    AC, info = exponentiate(x->applyH1(x, M[hn], F0, F[children]...), -im*dt/2, AC; ishermitian=true)

    #expand light-cone if necessary
    if in(hn, lc.edge)
        v = tensorcontract(conj(AC), collect(1:nc+2), lc.ref[hn], collect(1:nc+2))[1]
        if 1-norm(v) > lc.thresh
            filter!(x->x!=hn, lc.edge)
            push!(lc.edge, children...)
        else
            #evolve node forward half a time step
            AC, info = exponentiate(x->applyH1(x, M[hn], F[parent], F[children]...), -im*dt/2, AC; ishermitian=true)
            if verbose
                E = real(dot(AC, applyH1(AC, M[hn], F[parent], F[children]...)))
                println("Sweep R->L: AC on site $hn, energy = $E")
            end

            A[hn] = AC
            #node id has now been evolved one full time step
            return A, F
        end
    end

    if verbose
        E = real(dot(AC, applyH1(AC, M[hn], F0, F[children]...)))
        println("Sweep L->R: AC on site $hn, energy = $E")
    end

    for (i, child) in enumerate(children)
        grandchildren = A.tree[child].children
        otherchildren = filter(x->x!=child, children)
        ngc = length(grandchildren)

        #extract C from headnode
        AL, C = QR(AC, i+1)
        F[hn] = updateleftenv(AL, M[hn], i, F0, F[otherchildren]...)

        #evolve C backwards half a time step
        C, info = exponentiate(x->applyH0(x, F[hn], F[child]), im*dt/2, C; ishermitian=true)
        if verbose
            E = real(dot(C, applyH0(C, F[hn], F[child])))
            println("Sweep L->R: C between sites $hn and $child, energy = $E")
        end

        #contract C with child
        IA = collect(1:ngc+2)
        IB = collect(1:ngc+2)
        IA[1] = -1
        A[child] = tensorcontract(A[child], IA, C, [1,-1], IB)
        #(OC is now on child)
        
        #evolve child forward one full time step
        A, F = tdvp1sweep!(dt, A, M, F, lc, child; verbose=verbose, kwargs...)

        #extract C from child
        AR, C = QR(A[child], 1)
        A[child] = AR
        F[child] = updaterightenv(AR, M[child], F[grandchildren]...)
        
        #evolve C backwards half a time step
        C, info = exponentiate(x->applyH0(x, F[child], F[hn]), im*dt/2, C; ishermitian=true)
        if verbose
            E = real(dot(C, applyH0(C, F[child], F[hn])))
            println("Sweep R->L: C between sites $child and $hn, energy = $E")
        end

        #contract C with headnode
        IA = collect(1:nc+2)
        IB = collect(1:nc+2)
        IA[i+1] = -1
        AC = tensorcontract(AL, IA, C, [i+1,-1], IB)
        #(OC is now on headnode)
    end

    #evolve headnode forwards half a time step
    AC, info = exponentiate(x->applyH1(x, M[hn], F0, F[children]...), -im*dt/2, AC; ishermitian=true)
    if verbose
        E = real(dot(AC, applyH1(AC, M[hn], F0, F[children]...)))
        println("Sweep R->L: AC on site $hn, energy = $E")
    end
    A[hn] = AC
    #the tree has now been evolved by one full time step dt
    return A, F
end
function tdvp1sweep_lc!(dt, A::TreeNetwork, M::TreeNetwork, lc::TreeLightCone, F::Vector, id::Int; verbose=false, kwargs...)

    children = A.tree[id].children
    parent = A.tree[id].parent
    nc = length(children)
    AC = A[id]

    #(OC begins on node)
    #evolve node forward half a time step
    AC, info = exponentiate(x->applyH1(x, M[id], F[parent], F[children]...), -im*dt/2, AC; ishermitian=true)
    if verbose
        E = real(dot(AC, applyH1(AC, M[id], F[parent], F[children]...)))
        println("Sweep L->R: AC on site $id, energy = $E")
    end

    #expand light-cone if necessary
    if in(id, lc.edge)
        v = tensorcontract(conj(AC), collect(1:nc+2), lc.ref[id], collect(1:nc+2))[1]
        if 1-norm(v) > lc.thresh
            filter!(x->x!=id, lc.edge)
            push!(lc.edge, children...)
        else
            #evolve node forward half a time step
            AC, info = exponentiate(x->applyH1(x, M[id], F[parent], F[children]...), -im*dt/2, AC; ishermitian=true)
            if verbose
                E = real(dot(AC, applyH1(AC, M[id], F[parent], F[children]...)))
                println("Sweep R->L: AC on site $id, energy = $E")
            end

            A[id] = AC
            #node id has now been evolved one full time step
            return A, F
        end
    end

    for (i, child) in enumerate(children)
        grandchildren = A.tree[child].children
        otherchildren = filter(x->x!=child, children)
        ngc = length(grandchildren)

        #extract C from node
        AL, C = QR(AC, i+1)
        F[id] = updateleftenv(AL, M[id], i, F[parent], F[otherchildren]...)

        #evolve C backwards half a time step
        C, info = exponentiate(x->applyH0(x, F[id], F[child]), im*dt/2, C; ishermitian=true)
        if verbose
            E = real(dot(C, applyH0(C, F[id], F[child])))
            println("Sweep L->R: C between sites $id and $child, energy = $E")
        end

        #contract C with child
        IA = collect(1:ngc+2)
        IB = collect(1:ngc+2)
        IA[1] = -1
        A[child] = tensorcontract(A[child], IA, C, [1,-1], IB)
        #(OC is now on child)
        
        #evolve child forward one full time step
        A, F = tdvp1sweep!(dt, A, M, F, lc, child; verbose=verbose, kwargs...)

        #extract C from child
        AR, C = QR(A[child], 1)
        A[child] = AR
        F[child] = updaterightenv(AR, M[child], F[grandchildren]...)
        
        #evolve C backwards half a time step
        C, info = exponentiate(x->applyH0(x, F[child], F[id]), im*dt/2, C; ishermitian=true)
        if verbose
            E = real(dot(C, applyH0(C, F[child], F[id])))
            println("Sweep R->L: C between sites $child and $id, energy = $E")
        end

        #contract C with node
        IA = collect(1:nc+2)
        IB = collect(1:nc+2)
        IA[i+1] = -1
        AC = tensorcontract(AL, IA, C, [i+1,-1], IB)
        #(OC is now on node)
    end

    #evolve node forward half a time step
    AC, info = exponentiate(x->applyH1(x, M[id], F[parent], F[children]...), -im*dt/2, AC; ishermitian=true)
    if verbose
        E = real(dot(AC, applyH1(AC, M[id], F[parent], F[children]...)))
        println("Sweep R->L: AC on site $id, energy = $E")
    end

    A[id] = AC
    #node id has now been evolved one full time step
    return A, F
end

function productstatemps(tree_::Tree, physdims::Dims, Dmax::Int=1; state=:Vacuum)
    tree = deepcopy(tree_)
    hn = findheadnode(tree)
    leafnodes = leaves(tree)
    N = length(tree)
    setheadnode!(tree, leafnodes[1])
    bonddims1 = calcbonddims(tree, physdims, Dmax)
    setheadnode!(tree, leafnodes[2])
    bonddims2 = calcbonddims(tree, physdims, Dmax)
    bonddims = min.(bonddims1, bonddims2)
    tree = deepcopy(tree_)
    
    if state==:Vacuum
        statelist = [unitcol(1, physdims[i]) for i in 1:N]
    elseif state==:FullOccupy
        statelist = [unitcol(physdims[i], physdims[i]) for i in 1:N]
    elseif typeof(state)<:Vector
        statelist = state
    else
        throw(ErrorException("state input not recognised"))
    end

    A = Vector{Any}(undef, N)

    for (id, node) in enumerate(tree)
        Dpar = id==hn ? 1 : bonddims[id, node.parent]
        Dchildren = bonddims[id, node.children]
        B = zeros(ComplexF64, Dpar, prod(Dchildren), physdims[id])
        for j in 1:min(Dpar, prod(Dchildren))
            B[j,j,:] = statelist[id]
        end
        B = reshape(B, Dpar, Dchildren..., physdims[id])
        A[id] = B
    end
    net = TreeNetwork(tree, A)
    mpsrightnorm!(net)
    return net
end
productstatemps(tree::Tree, physdims::Int, Dmax::Int; state=:Vacuum) =
    productstatemps(tree, ntuple(i -> physdims, length(tree)), Dmax; state=state)

function mpsembed!(A::TreeNetwork, Dmax::Int)
    tree = deepcopy(A.tree)
    pdims = physdims(A)
    hn = findheadnode(tree)
    leafnodes = leaves(tree)
    setheadnode!(tree, leafnodes[1])
    bonddims1 = calcbonddims(tree, pdims, Dmax)
    setheadnode!(tree, leafnodes[2])
    bonddims2 = calcbonddims(tree, pdims, Dmax)
    bonddims = min.(bonddims1, bonddims2)

    for (id, nd) in enumerate(A.tree.nodes)
        parent = nd.parent
        children = nd.children
        if parent != 0
            A[id] = setbond(A[id], bonddims[id, [parent, children...]]...)
        else
            A[id] = setbond(A[id], 1, bonddims[id, children]...)
        end
    end

    mpsrightnorm!(A)
    return A
end

function bonddims(A::TreeNetwork)
    N = length(A)
    mat = zeros(Int, N, N)
    for bond in bonds(A)
        id1 = bond[1]
        id2 = bond[2]
        dir = findbond(A.tree[id1], id2)
        D = size(A[id1])[dir]
        mat[bond...] = D
        mat[reverse(bond)...] = D
    end
    mat
end
