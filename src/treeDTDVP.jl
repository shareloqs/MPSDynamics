function initenvs_full!(A::TreeNetwork, M::TreeNetwork, Afull::Vector, F::Vector, id::Int; Dplusmax=nothing)
    for child in A.tree[id].children
        F, Afull = initenvs_full!(A, M, Afull, F, child)
    end

    par = A.tree[id].parent

    if par != 0
        dims = size(A[id])
        aleft = dims[1]
        aright = dims[2:end-1]
        d = dims[end]
        
        nc = length(A.tree[par].children)
        dimspar = size(A[par])
        dir = findbond(A.tree[par], id)
        otherlegs = filter(x->x!=dir, collect(1:nc+1))
        Dpar = dimspar[otherlegs]
        dpar = dimspar[end]
        Dmaxleft = prod(Dpar)*dpar

        Dmaxright = Dplusmax != nothing ? min(prod(aright)*d, aleft+Dplusmax) : prod(aright)*d

        Dmax = min(Dmaxright, Dmaxleft)
        Afull[id], C = QR_full(A[id], 1)
        F[id] = updaterightenv(truncAR(Afull[id], Dmax), M[id], truncF.(F[A.tree[id].children], aright)...)
    end
    return F, Afull
end
initenvs_full(A::TreeNetwork, M::TreeNetwork; Dplusmax=nothing) =
    initenvs_full!(A, M, Vector{Any}(undef, length(A)), Vector{Any}(undef, length(A)), findheadnode(A); Dplusmax=Dplusmax)

tdvpdsweep!(dt, A::TreeNetwork, M::TreeNetwork; Dplusmax=nothing, verbose=false, kwargs...) =
    tdvpdsweep!(dt, A, M, initenvs_full(A, M, Dplusmax=Dplusmax)..., findheadnode(A); Dplusmax=Dplusmax, verbose=verbose, kwargs...)

tdvpdsweep!(dt, A::TreeNetwork, M::TreeNetwork, FR::Vector, Afull::Vector; Dplusmax=nothing, verbose=false, kwargs...) =
    tdvpdsweep!(dt, A, M, FR, Afull, findheadnode(A); Dplusmax=Dplusmax, verbose=verbose, kwargs...)
    

function tdvpdsweep!(dt, A::TreeNetwork, M::TreeNetwork, F::Vector, Afull::Vector, id::Int; Dplusmax=nothing, verbose=false, kwargs...)
    parent = A.tree[id].parent
    children = A.tree[id].children
    nc = length(children)
    AC = A[id]
    F0 = parent==0 ? fill!(similar(M[1], (1,1,1)), 1) : F[parent]

    Dnew = getindex.(size.(F[children]), 1)
    Dold = size(AC)[2:end-1]
    AC = setbond(AC, Dnew...)

    AC, info = exponentiate(x->applyH1(x, M[id], F0, F[children]...), -im*dt/2, AC; ishermitian=true, kwargs...)
    if verbose
        E = real(dot(AC, applyH1(AC, M[id], F0, F[children]...)))
        println("Sweep L->R: AC on site $id, energy = $E")
        for (i, child) in enumerate(children)
            Dnew[i]!=Dold[i] && println("*BondDimension $id-$child changed from $(Dold[i]) to $(Dnew[i])")
            Dnew[i]==Dold[i] && println("*BondDimension $id-$child constant at $(Dnew[i])")
        end
    end

    for (i, child) in enumerate(children)
        grandchildren = A.tree[child].children
        otherchildren = filter(x->x!=child, children)

        Dotherchildren = size(AC)[[collect(2:i)..., collect(i+2:nc+1)...]]
        D0 = size(A[id])[1]

        AL, C = QR(AC, i+1)
        F[id] = updateleftenv(AL, M[id], i, truncF(F0, D0), truncF.(F[otherchildren], Dotherchildren)...)

        C, info = exponentiate(x->applyH0(x, F[id], F[child]), im*dt/2, C; ishermitian=true, kwargs...)
        if verbose
            E = real(dot(C, applyH0(C, F[id], F[child])))
            println("Sweep L->R: C between sites $id and $child, energy = $E")
        end

        A[child] = contractC!(truncAR(Afull[child], Dnew[i]), C, 1)

        A, Afull, F = tdvpdsweep!(dt, A, H, F, Afull, child; Dplusmax=Dplusmax, verbose=verbose, kwargs...)

        Dr = Int[getindex.(size.(A[grandchildren]), 1)...]
        d = size(A[child])[end]
        size(A[id])
        ###max from other side
        Dmax = (Dplusmax != nothing ? min(prod(Dr)*d, Dnew[i]+Dplusmax) : prod(Dr)*d)

        Afull[child], C = QR_full(A[child], 1)
        A[child] = truncAR(Afull[child], Dnew[i])
        F[child] = updaterightenv(truncAR(Afull[child], Dmax), M[child], truncF.(F[grandchildren], Dr)...)

        C, info = exponentiate(x->applyH0(x, truncF(F[child], Dnew[i]), F[id]), im*dt/2, C; ishermitian=true, kwargs...)
        if verbose
            E = real(dot(C, applyH0(C, truncF(F[child], Dnew[i]), F[id])))
            println("Sweep R->L: C between sites $child and $id, energy = $E")
        end

        AC = contractC!(AL, C, i+1)
    end

    AC, info = exponentiate(x->applyH1(x, M[id], F0, truncF.(F[children], Dnew)...), -im*dt/2, AC; ishermitian=true, kwargs...)
    if verbose
        E = real(dot(AC, applyH1(AC, M[id], F0, truncF.(F[children], Dnew)...)))
        println("Sweep R->L: AC on site $id, energy = $E")
    end
    A[id] = AC

    return A, Afull, F
end

function tdvp1sweep_dynamic!(dt, A::TreeNetwork, H::TreeNetwork, Afull=nothing, FRs=nothing; obs=Vector{Observable}, prec=10^-3, Dlim=50, verbose=false, error=false, timed=false, Dplusmax=nothing, kwargs...)
    

end

function updatebonds(F, A::TreeNetwork, M::TreeNetwork, prec, Dlim)

    N = length(A)
    olddims = bonddims(A)

    ACs = Vector{Any}(undef, N)
    PAs = Vector{Any}(undef, N)
    PCs = Vector{Any}(undef, N-1)

    LA = fill!(similar(M[1], (1,1,1)), 1)

    hn = findheadnode(A)
    AC = A[hn]
    ACs[hn] = AC

    for (dir, id) in Directional(Traverse(A.tree))
        children = A.tree[id].children
        parent = A.tree[id].parent
        if length(children) != 0
            next = children[dir-1]
            otherchildren = filter(x->x!=next, children)
            Dl = [parent!=0 ? olddims[id, parent] : 1, olddims[id, otherchildren]...]
            Dr = olddims[id, next]
            d = size(A[id])[end]
            Dmax = min(prod(Dl)*d, size(F)[1])

            println("id=$id")
            println("dir=$dir")
            println(otherchildren)
            println(size(LA))
            L = LA_M_A(LA, M[id], A[id], dir, truncF2.(F[otherchildren], Dl[2:end])...)
            PAs[id] = L_FR(L, truncF2(F[next], Dr))
            AL, C = QR_full(AC, dir)
            AC = contractC!(A[next], C, 1)##make default method
            ACs[next] = AC
            LA = L_AL(L[[1:D for D in Dl]...,:,:,:], AL, Dmax, dir)
            PCs[next] = LA_FR(LA, F[next])
        else
            L = LA_M_A(LA, M[id], A[id])
            PAs[id] = L
        end
    end
    
    return FR, ACs, PAs, PCs
end
