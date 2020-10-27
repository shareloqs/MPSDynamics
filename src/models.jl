function xyzmpo(N::Int; Jx=1.0, Jy=Jx, Jz=Jx, hx=0., hz=0.)
    u = unitmat(2)

    D = 5 - iszero(Jx) - iszero(Jy) - iszero(Jz)
    M = zeros(ComplexF64, D, D, 2, 2)
    M[1,1,:,:] = M[D,D,:,:] = u
    M[1,D,:,:] = -hx*sx + -hz*sz
    i = 2
    if !iszero(Jx)
        M[1,i,:,:] = -Jx*sx
        M[i,D,:,:] = sx
        i += 1
    end
    if !iszero(Jy)
        M[1,i,:,:] = -Jy*sy
        M[i,D,:,:] = sy
        i += 1
    end
    if !iszero(Jz)
        M[1,i,:,:] = -Jz*sz
        M[i,D,:,:] = sz
        i += 1
    end
    return Any[M[1:1,:,:,:], fill(M, N-2)..., M[:,D:D,:,:]]
end

isingmpo(N::Int; J=1.0, h=1.0) = xyzmpo(N; Jx=J, Jy=0., Jz=0., hz=h, hx=0.)
heisenbergmpo(N::Int, J=1.0) = xyzmpo(N; Jx=J)
xxzmpo(N::Int, Δ = 1.0, J=1.0) = xyzmpo(N; Jx=J, Jy=J, Jz=J*Δ)

function longrange_xyzmpo(N::Int, α::Float64=0.; Jx=1.0, Jy=Jx, Jz=Jx, hx=0., hz=0.)
    u = unitmat(2)

    l = exp(-α)

    D = 5 - iszero(Jx) - iszero(Jy) - iszero(Jz)
    M = zeros(ComplexF64, D, D, 2, 2)
    M[1,1,:,:] = M[D,D,:,:] = u
    M[1,D,:,:] = -hx*sx + -hz*sz
    i = 2
    if !iszero(Jx)
        M[1,i,:,:] = -l*Jx*sx
        M[i,D,:,:] = sx
        M[i,i,:,:] = u*l
        i += 1
    end
    if !iszero(Jy)
        M[1,i,:,:] = -l*Jy*sy
        M[i,D,:,:] = sy
        M[i,i,:,:] = u*l
        i += 1
    end
    if !iszero(Jz)
        M[1,i,:,:] = -l*Jz*sz
        M[i,D,:,:] = sz
        M[i,i,:,:] = u*l
        i += 1
    end
    return Any[M[1:1,:,:,:], fill(M, N-2)..., M[:,D:D,:,:]]
end
longrange_isingmpo(N::Int, α::Float64=0.; J=1.0, h=1.0) = longrange_xyzmpo(N, α; Jx=J, Jy=0., Jz=0., hz=h, hx=0.)

function spinchainmpo(N::Int; J=1.0, hz=1.0, hx=0.0, i=div(N,2))
    u = unitmat(2)
        
    Mx = zeros(4,4,2,2)
    Mx[1,1,:,:] = u
    Mx[1,2,:,:] = J*sp
    Mx[1,3,:,:] = J*sm
    Mx[1,4,:,:] = hz*sz + hx*sx
    Mx[2,4,:,:] = sm
    Mx[3,4,:,:] = sp
    Mx[4,4,:,:] = u

    M0 = zeros(4,4,2,2)
    M0[1,1,:,:] = u
    M0[1,2,:,:] = J*sp
    M0[1,3,:,:] = J*sm
    M0[1,4,:,:] = hz*sz
    M0[2,4,:,:] = sm
    M0[3,4,:,:] = sp
    M0[4,4,:,:] = u

    return Any[M0[1:1,:,:,:], fill(M0, i-2)..., Mx, fill(M0, N-i-1)..., M0[:,4:4,:,:]]
end

function tightbindingmpo(N::Int, d::Int; J=1.0, e=1.0)
    b = anih(d)
    bd = crea(d)
    u = unitmat(d)    
    n = numb(d)

    M = zeros(4,4,d,d)
    M[1,1,:,:] = u
    M[1,2,:,:] = J*b
    M[1,3,:,:] = J*bd
    M[1,4,:,:] = e*n
    M[2,4,:,:] = bd
    M[3,4,:,:] = b
    M[4,4,:,:] = u

    return Any[M[1:1,:,:,:], fill(M, N-2)..., M[:,4:4,:,:]]
end

function hbathchain(N::Int, d::Int, chainparams; tree=false, reverse=false)
    b = anih(d)
    bd = crea(d)
    n = numb(d)
    u = unitmat(d)
    e = chainparams[1]
    t = chainparams[2]

    up(H, h, hd) = permutedims(cat(H, h, hd, u; dims=3), [3,1,2])
    dn(H, h, hd) = permutedims(cat(u, hd, h, H; dims=3), [3,1,2])

    H=Vector{Any}()
    for i=1:N-1
        M=zeros(4, 4, d, d)
        M[4, :, :, :] = up(e[i]*n, t[i]*b, t[i]*bd)
        M[:, 1, :, :] = dn(e[i]*n, b, bd)
        push!(H, M)
    end
    M=zeros(4, d, d)
    M[:, :, :] = dn(e[N]*n, b, bd)
    push!(H, M)
    if tree
        return TreeNetwork(H)
    else
        H[end] = reshape(H[end],4,1,d,d)
        reverse && reversempo!(H)
        return H
    end
end
hbathchain(N::Int, d::Int, e::Int, t::Int; tree=false) = hbathchain(N, d, (fill(e, N), fill(t, N-1), nothing); tree=tree)

function methylbluempo(e1, e2, δ, N1, N2, N3, d1, d2, d3, chainparams1, chainparams2, chainparams3)
    u = unitmat(3)
    
    up(H, h, hd) = permutedims(cat(H, h, hd, u; dims=3), [3,1,2])

    c1 = chainparams1[3]
    c2 = chainparams2[3]
    c3 = chainparams3[3]

    s2 = unitcol(1, 3)
    s1 = unitcol(2, 3)
       
    #Hs = e1*s1*s1' + e2*s2*s2' + δ*(s1*s2' + s2*s1')
    Hs = (e2-e1)*s2*s2' + δ*(s1*s2' + s2*s1') # e^(-is1*s1't)He^(is1*s1't)

    M=zeros(1,4,4,4,3,3)
    M[1, :, 1, 1, :, :] = up(Hs, c1*s1*s1', c1*s1*s1')
    M[1, 1, :, 1, :, :] = up(Hs, c2*s2*s2', c2*s2*s2')
    M[1, 1, 1, :, :, :] = up(Hs, c3*(s1*s2'+s2*s1'), c3*(s1*s2'+s2*s1'))

    H = TreeNetwork(Any[M])
    addtree!(H, 1, hbathchain(N1, d1, chainparams1, tree=true))
    addtree!(H, 1, hbathchain(N2, d2, chainparams2, tree=true))
    addtree!(H, 1, hbathchain(N3, d3, chainparams3, tree=true))
    return H
end

function methylblue_S1_mpo(e1, N, d, chainparams; tree=false)
    u = unitmat(2)

    up(H, h, hd) = permutedims(cat(H, h, hd, u; dims=3), [3,1,2])

    c = chainparams[3]
    s1 = unitcol(1, 2)
    #Hs = e1*s1*s1'
    Hs = zero(u) # e^(-is1*s1't)He^(is1*s1't)

    M=zeros(1,4,2,2)
    M[1, :, :, :] = up(Hs, c*s1*s1', c*s1*s1')

    if tree
        chain = hbathchain(N, d, chainparams; tree=true)
        H = TreeNetwork(Any[M])
        addtree!(H, 1, chain)
        return H
    else
        chain = hbathchain(N, d, chainparams; tree=false)
        return [M, chain...]
    end
end

function spinbosonmpo(ω0, Δ, d, N, chainparams; rwa=false, tree=false)
    u = unitmat(2)
    
    up(H, h, hd) = permutedims(cat(H, h, hd, u; dims=3), [3,1,2])

    c0 = chainparams[3]

    Hs = (ω0/2)*sz + Δ*sx

    M=zeros(1,4,2,2)
    M[1, :, :, :] = up(Hs, rwa ? c0*sm : c0*sx, rwa ? c0*sp : c0*sx)

    if tree
        chain = hbathchain(N, d, chainparams; tree=true)
        H = TreeNetwork(Vector{AbstractArray}([M]))
        addtree!(H, 1, chain)
        return H
    else
        chain = hbathchain(N, d, chainparams; tree=false)
        return Any[M, chain...]
    end
end

function twobathspinmpo(ω0, Δ, Nl, Nr, dl, dr, chainparamsl=[fill(1.0,N),fill(1.0,N-1), 1.0], chainparamsr=chainparamsl; rwar=false, rwal=false, tree=false)
    u = unitmat(2)

    up(H, h, hd) = permutedims(cat(H, h, hd, u; dims=3), [3,1,2])

    cl = chainparamsl[3]
    cr = chainparamsr[3]

    Hs = (ω0/2)*sz + Δ*sx

    Al = rwal ? cl*sm : cl*sx
    Ar = rwar ? cr*sm : cr*sx

    M=zeros(4,4,2,2)
    M[1,:,:,:] = up(Hs, Ar, Matrix(Ar'))
    M[:,1,:,:] = up(Hs, Al, Matrix(Al'))

    return Any[hbathchain(Nl, dl, chainparamsl; reverse=true)..., M, hbathchain(Nr, dr, chainparamsr)...]
end

"""
    chaincoeffs_ohmic(nummodes, α, s, beta="inf"; wc=1, soft=false)

Generate chain coefficients for an Harmonic bath coupled to a spin-1/2 with spectral density given by: 

soft cutoff: ``J(ω) = 2παω_c (\\frac{ω}{ω_c})^s \\exp(-ω/ω_c)`` \n
hard cutoff: ``J(ω) = 2παω_c (\\frac{ω}{ω_c})^s θ(ω-ω_c)``

The Hamiltonian is given by:

``H = \\frac{ω_0}{2}σ_z + Δσ_x + σ_x\\sum_kg_k(b_k^\\dagger+b_k) + \\sum_kω_kb_k^\\dagger b_k``

And the spectral density is defined by:

``J(ω) ≡ π\\sum_k|g_k|^2δ(ω-ω_k)``
"""
function chaincoeffs_ohmic(nummodes, α, s, beta="inf"; ωc=1, soft=false)
    if beta=="inf"
        if soft
            c0 = ωc*sqrt(2*α*gamma(s+1))
            e = [ωc*(2n + 1 + s) for n in 0:(nummodes-1)]
            t = [ωc*sqrt((n + 1)*(n + s + 1)) for n in 0:(nummodes-2)]
            return e, t, c0
        else    
            c0 = sqrt(2α/(s+1))*ωc
            e = [(ωc/2)*(1 + (s^2)/((s+2n)*(2+s+2n))) for n in 0:(nummodes-1)]
            t = [ωc*(1+n)*(1+s+n)/((s+2+2n)*(s+3+2n))*sqrt((3+s+2n)/(1+s+2n)) for n in 0:(nummodes-2)]
            return e, t, c0
        end
    else
        if soft
            throw(ErrorException("no data for soft cut-off at finite beta"))
        end
        return getchaincoeffs(nummodes, α, s, beta, ωc)
    end
end

function getchaincoeffs(nummodes, α, s, beta, ωc=1)
    matlabdir = MATDIR
    astr = paramstring(α, 2)
    bstr = paramstring(beta, 3)
    datfname = string("jacerg_","a",astr,"s$s","beta",bstr,".csv")
    chaincoeffs = readdlm(string(matlabdir,datfname),',',Float64,'\n')
    es = ωc .* chaincoeffs[:,1]
    ts = ωc .* chaincoeffs[1:end-1,2]
    c0 = ωc * sqrt(chaincoeffs[end,2]/pi)
    Nmax = length(es)
    if Nmax < nummodes
        throw(ErrorException("no data for nummodes > $Nmax"))
    end
    return [es[1:nummodes], ts[1:nummodes-1], c0]
end


function ibmmpo(ω0, d, N, chainparams; tree=false)
    u = unitmat(2)
    
    up(H, h, hd) = permutedims(cat(H, h, hd, u; dims=3), [3,1,2])

    c0 = chainparams[3]

    Hs = (ω0/2)*sz

    M=zeros(1,4,2,2)
    M[1, :, :, :] = up(Hs, c0*sz, c0*sz)

    chain = hbathchain(N, d, chainparams)
    if tree
        H = TreeNetwork(Vector{AbstractArray}([M]))
        addtree!(H, 1, chain)
        return H
    else
        chain = chain.sites
        chain[end] = reshape(chain[end],4,1,d,d)
        return Any[M, chain...]
    end
end

function tunnelingmpo(ϵ, delta, α, s, β, d::Int, nummodes::Int; tree=false, ωc=1)
    cps = chaincoeffs_ohmic(nummodes, α, s, β; ωc=ωc)
    λ = 2*α*ωc/s + delta
    
    u = unitmat(2)
    
    up(H, h, hd) = permutedims(cat(H, h, hd, u; dims=3), [3,1,2])

    c0 = cps[3]

    Hs = (ϵ/2)*sz + λ*(u + sx)/2

    M=zeros(1,4,2,2)
    M[1, :, :, :] = up(Hs, c0*(u + sx)/2, c0*(u + sx)/2)

    if tree
        chain = hbathchain(nummodes, d, cps; tree=true)
        H = TreeNetwork(Vector{AbstractArray}([M]))
        addtree!(H, 1, chain)
        return H
    else
        chain = hbathchain(nummodes, d, cps; tree=false)
        return Any[M, chain...]
    end
end

function nearestneighbourmpo(N::Int, h0, A, Ad = A')
    size(h0) == size(A) || error("physical dimensions don't match")
    size(h0) == size(Ad) || error("physical dimensions don't match")
    d = size(h0)[1]

    D = A == Ad ? 3 : 4

    u=unitmat(d)

    up(H, h, hd) = permutedims(cat(H, h, hd, u; dims=3), [3,1,2])
    dn(H, h, hd) = permutedims(cat(u, hd, h, H; dims=3), [3,1,2])
    up(H, h) = permutedims(cat(H, h, u; dims=3), [3,1,2])
    dn(H, h) = permutedims(cat(u, h, H; dims=3), [3,1,2])

    M=zeros(D, D, d, d)
    M[D, :, :, :] = A == Ad ? up(h0, A) : up(h0, A, Ad)
    M[:, 1, :, :] = A == Ad ? dn(h0, A) : dn(h0, A, Ad)

    return Any[M[D:D,:,:,:], fill(M, N-2)..., M[:,1:1,:,:]]
end

function nearestneighbourmpo(tree_::Tree, h0, A, Ad = A')
    size(h0) == size(A) || error("physical dimensions don't match")
    size(h0) == size(Ad) || error("physical dimensions don't match")
    d = size(h0)[1]

    tree = deepcopy(tree_)

    D = A == Ad ? 3 : 4

    u=unitmat(d)

    up(H, h, hd) = permutedims(cat(H, h, hd, u; dims=3), [3,1,2])
    dn(H, h, hd) = permutedims(cat(u, hd, h, H; dims=3), [3,1,2])
    up(H, h) = permutedims(cat(H, h, u; dims=3), [3,1,2])
    dn(H, h) = permutedims(cat(u, h, H; dims=3), [3,1,2])

    N = length(tree)
    hn = findheadnode(tree)

    Ms = Vector{Any}(undef, N)

    for (id, nd) in enumerate(tree)
        nc = length(nd.children)
        dims = [fill(D, nc+1)..., d, d]
        M = zeros(dims...)
        M[:, fill(1, nc)..., :, :] = A == Ad ? dn(h0, A) : dn(h0, A, Ad)
        for i in 1:nc
            M[D, fill(1, i-1)..., :, fill(1, nc-i)..., :, :] = A == Ad ? up(h0, A) : up(h0, A, Ad)
        end
        Ms[id] = M
    end
    nc = length(tree[hn].children)
    Ms[hn] = Ms[hn][D:D, fill(:,nc+2)...]
    return TreeNetwork(tree, Ms)
end
