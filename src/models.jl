up(ops...) = permutedims(cat(ops...; dims=3), [3,1,2])
dn(ops...) = permutedims(cat(reverse(ops)...; dims=3), [3,1,2])

"""
    xyzmpo(N::Int; Jx=1.0, Jy=Jx, Jz=Jx, hx=0., hz=0.)

Return the MPO representation of the `N`-spin XYZ Hamiltonian with external field ``\\vec{h}=(h_x, 0, h_z)``.

``
H = \\sum_{n=1}^{N-1} -J_x\\sigma_x^{n} \\sigma_x^{n+1} - J_y\\sigma_y^{n} \\sigma_y^{n+1} - J_z\\sigma_z^{n} \\sigma_z^{n+1} + \\sum_{n=1}^{N}(- h_x\\sigma_x^{n} - h_z\\sigma_z^{n})  
``.

"""
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

"""
    isingmpo(N; J=1.0, h=1.0)

Return the MPO representation of a `N`-spin 1D Ising model with external field ``\\vec{h} = (0,0,h)``.
"""
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

"""
    hbathchain(N::Int, d::Int, chainparams, longrangecc...; tree=false, reverse=false, coupletox=false)

Create an MPO representing a tight-binding chain of `N` oscillators with `d` Fock states each. Chain parameters are supplied in the standard form: `chainparams` ``=[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]``. The output does not itself represent a complete MPO but will possess an end which is *open* and should be attached to another tensor site, usually representing the *system*.

# Arguments

* `reverse`: If `reverse=true` create a chain were the last (i.e. Nth) site is the site which couples to the system
* `coupletox`: Used to choose the form of the system coupling. `coupletox=true` gives a non-number conserving coupling of the form ``H_{\\text{I}}=A_{\\text{S}}(b_{0}^\\dagger + b_0)`` where ``A_{\\text{S}}`` is a system operator, while `coupletox=false` gives the number-converving coupling ``H_{\\text{I}}=(A_{\\text{S}} b_{0}^\\dagger + A_{\\text{S}}^\\dagger b_0)``
* `tree`: If `true` the resulting chain will be of type `TreeNetwork`; useful for construcing tree-MPOs 

# Example

One can constuct a system site tensor to couple to a chain by using the function `up` to populate the tensor. For example, to construct a system site with Hamiltonian `Hs` and coupling operator `As`, the system tensor `M` is constructed as follows for a non-number conserving interaction:
```julia
u = one(Hs) # system identity
M = zeros(1,3,2,2)
M[1, :, :, :] = up(Hs, As, u)
```

The full MPO can then be constructed with:
```julia
Hmpo = [M, hbathchain(N, d, chainparams, coupletox=true)...]
```

Similarly for a number conserving interaction the site tensor would look like:
```julia
u = one(Hs) # system identity
M = zeros(1,4,2,2)
M[1, :, :, :] = up(Hs, As, As', u)
```
And the full MPO would be
```julia
Hmpo = [M, hbathchain(N, d, chainparams; coupletox=false)...]
```
    
"""
function hbathchain(N::Int, d::Int, chainparams, longrangecc...; tree=false, reverse=false, coupletox=false)
    b = anih(d)
    bd = crea(d)
    n = numb(d)
    u = unitmat(d)
    e = chainparams[1]
    t = chainparams[2]

    numlong = length(longrangecc)
    cc = longrangecc
    D1 = 2 + (coupletox ? 1 : 2)*(1+numlong)
    D2 = coupletox ? D1+1 : D1 

    H=Vector{Any}()

    if coupletox
        M=zeros(D1, D2, d, d)
        M[D1, :, :, :] = up(e[1]*n, t[1]*b, t[1]*bd, fill(zero(u), numlong)..., u)
        M[:, 1, :, :] = dn(e[1]*n, [cc[j][1]*(b+bd) for j=1:numlong]..., b+bd, u)
        for k=1:numlong
            M[k+2,k+3,:,:] = u
        end
        push!(H, M)
    else
        M=zeros(D1, D2, d, d)
        M[D1, :, :, :] = up(e[1]*n, t[1]*b, t[1]*bd, u)
        M[:, 1, :, :] = dn(e[1]*n, b, bd, u)
        numlong > 0 && error("haven't yet coded case of long range couplings with non-hermitian coupling")
        push!(H, M)
    end
    for i=2:N-1
        M=zeros(D2, D2, d, d)
        M[D2, :, :, :] = up(e[i]*n, t[i]*b, t[i]*bd, fill(zero(u), numlong)..., u)
        M[:, 1, :, :] = dn(e[i]*n, [cc[j][i]*(b+bd) for j=1:numlong]..., b, bd, u)
        for k=1:numlong
            M[k+3,k+3,:,:] = u
        end
        push!(H, M)
    end
    M=zeros(D2, d, d)
    M[:, :, :] = dn(e[N]*n, [cc[j][N]*(b+bd) for j=1:numlong]..., b, bd, u)
    push!(H, M)
    if tree
        return TreeNetwork(H)
    else
        H[end] = reshape(H[end], D2, 1, d, d)
        reverse && reversempo!(H)
        return H
    end
end
hbathchain(N::Int, d::Int, e::Int, t::Int; tree=false) = hbathchain(N, d, (fill(e, N), fill(t, N-1), nothing); tree=tree)

function methylbluempo(e1, e2, δ, N1, N2, N3, d1, d2, d3, cparS1, cparS2, cparS1S2)
    u = unitmat(3)
    
    c1 = only(cparS1[3])
    c2 = only(cparS2[3])
    c3 = only(cparS1S2[3])

    s2 = unitcol(1, 3)
    s1 = unitcol(2, 3)
       
    #Hs = e1*s1*s1' + e2*s2*s2' + δ*(s1*s2' + s2*s1')
    Hs = (e2-e1)*s2*s2' + δ*(s1*s2' + s2*s1') # e^(-is1*s1't)He^(is1*s1't)

    M=zeros(1,3,3,3,3,3)
    M[1, :, 1, 1, :, :] = up(Hs, c1*s1*s1', u)
    M[1, 1, :, 1, :, :] = up(Hs, c2*s2*s2', u)
    M[1, 1, 1, :, :, :] = up(Hs, c3*(s1*s2'+s2*s1'), u)

    H = TreeNetwork(Any[M])
    addtree!(H, 1, hbathchain(N1, d1, cparS1, tree=true, coupletox=true))
    addtree!(H, 1, hbathchain(N2, d2, cparS2, tree=true, coupletox=true))
    addtree!(H, 1, hbathchain(N3, d3, cparS1S2, tree=true, coupletox=true))

    return H
end

function methylbluempo2(e1, e2, δ, N1, N2, N3, d1, d2, d3, S1a1, S2a1, S1a2, S2a2, cparS1S2)
    u = unitmat(3)

    c1 = only(S1a1[3])
    c2 = only(S2a2[3])
    c3 = only(cparS1S2[3])

    s2 = unitcol(1, 3)
    s1 = unitcol(2, 3)

    #Hs = e1*s1*s1' + e2*s2*s2' + δ*(s1*s2' + s2*s1')
    Hs = (e2-e1)*s2*s2' + δ*(s1*s2' + s2*s1') # e^(-is1*s1't)He^(is1*s1't)
    M = zeros(1,4,4,3,3,3)
    M[1,:,1,1,:,:] = up(Hs, c1*s1*s1', s2*s2', u)
    M[1,1,:,1,:,:] = up(Hs, c2*s2*s2', s1*s1', u)
    M[1,1,1,:,:,:] = up(Hs, c3*(s1*s2'+s2*s1'), u)

    H = TreeNetwork(Any[M])
    addtree!(H, 1, hbathchain(N1, d1, S1a1, S2a1; coupletox=true, tree=true))
    addtree!(H, 1, hbathchain(N2, d2, S2a2, S1a2; coupletox=true, tree=true))
    addtree!(H, 1, hbathchain(N3, d3, cparS1S2; coupletox=true, tree=true))
    return H
end

function methylbluempo_correlated(e1, e2, δ, N1, N2, d1, d2, cparS1, ccS2, cparS1S2)
    u = unitmat(3)
    s2 = unitcol(1, 3)
    s1 = unitcol(2, 3)

    c1 = only(cparS1[3])
    c3 = only(cparS1S2[3])

    # Hs = system Hamiltonian 
    Hs = (e2-e1)*s2*s2' + δ*(s1*s2' + s2*s1') # e^(-is1*s1't)He^(is1*s1't)

    M = zeros(3,4,3,3)
    M[:, 1, :, :] = up(Hs, c3*(s1*s2'+s2*s1'), u)
    M[1, :, :, :] = up(Hs, c1*(s1*s1'), s2*s2', u)

    chain1 = hbathchain(N1, d1, cparS1S2; coupletox=true, reverse=true)
    chain2 = hbathchain(N2, d2, cparS1, ccS2; coupletox=true)

    return Any[chain1..., M, chain2...]
end

function methylbluempo_correlated_nocoupling(e1, e2, N, d, cparS1, ccS2)
    u = unitmat(3)
    s2 = unitcol(1, 3)
    s1 = unitcol(2, 3)

    c = only(cparS1[3])

    Hs = (e2-e1)*s2*s2'

    M=zeros(1,4,3,3)
    M[1,:,:,:] = up(Hs, c*s1*s1', s2*s2', u)

    chain = hbathchain(N, d, cparS1, ccS2; coupletox=true)

    return Any[M, chain...]
end
function methylbluempo_nocoupling(e1, e2, N1, N2, d1, d2, cparS1, cparS2)
    u = unitmat(3)
    s2 = unitcol(1, 3)
    s1 = unitcol(2, 3)
    
    c1 = only(cparS1[3])
    c2 = only(cparS2[3])

    Hs = (e2-e1)*s2*s2'

    M=zeros(3,3,3,3)
    M[:,1,:,:] = up(Hs, c1*s1*s1', u)
    M[1,:,:,:] = up(Hs, c2*s2*s2', u)

    chain1 = hbathchain(N1, d1, cparS1; coupletox=true, reverse=true)
    chain2 = hbathchain(N2, d2, cparS2; coupletox=true)

    return Any[chain1..., M, chain2...]
end

function methylblue_S1_mpo(e1, N, d, chainparams; tree=false)
    u = unitmat(2)

    c = only(chainparams[3])
    s1 = unitcol(1, 2)
#    Hs = e1*s1*s1'
    Hs = zero(u) # e^(-is1*s1't)He^(is1*s1't)

    M=zeros(1,3,2,2)
    M[1, :, :, :] = up(Hs, c*s1*s1', u)

    if tree
        chain = hbathchain(N, d, chainparams; coupletox=true, tree=true)
        H = TreeNetwork(Any[M])
        addtree!(H, 1, chain)
        return H
    else
        chain = hbathchain(N, d, chainparams; coupletox=true, tree=false)
        return [M, chain...]
    end
end

"""
    spinbosonmpo(ω0, Δ, d, N, chainparams; rwa=false, tree=false)

Generate MPO for a spin-1/2 coupled to a chain of harmonic oscillators, defined by the Hamiltonian

``
H = \\frac{ω_0}{2}σ_z + Δσ_x + c_0σ_x(b_k^\\dagger+b_k) + \\sum_{i=0}^{N-1} t_i (b_{i+1}^\\dagger b_i +h.c.) + \\sum_{i=0}^{N-1} ϵ_ib_i^\\dagger b_i
``.

The spin is on site 1 of the MPS and the bath modes are to the right.

This Hamiltonain is unitarily equivalent (before the truncation to `N` sites) to the spin-boson Hamiltonian defined by

``
H =  \\frac{ω_0}{2}σ_z + Δσ_x + σ_x\\int_0^∞ dω\\sqrt{J(ω)}(b_ω^\\dagger+b_ω) + \\int_0^∞ωb_ω^\\dagger b_ω
``.

The chain parameters, supplied by `chainparams`=``[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]``, can be chosen to represent any arbitrary spectral density ``J(ω)`` at any temperature.

The rotating wave approximation can be made by setting `rwa=true`.

"""
function spinbosonmpo(ω0, Δ, d, N, chainparams; rwa=false, tree=false)
    u = unitmat(2)
    
    c0 = only(chainparams[3])

    Hs = (ω0/2)*sz + Δ*sx

    if rwa
        M=zeros(1,4,2,2)
        M[1, :, :, :] = up(Hs, c0*sm, c0*sp, u)
    else
        M=zeros(1,3,2,2)
        M[1, :, :, :] = up(Hs, c0*sx, u)
    end

    if tree
        chain = hbathchain(N, d, chainparams; tree=true, coupletox=!rwa)
        H = TreeNetwork(Vector{AbstractArray}([M]))
        addtree!(H, 1, chain)
        return H
    else
        chain = hbathchain(N, d, chainparams; tree=false, coupletox=!rwa)
        return Any[M, chain...]
    end
end

function twobathspinmpo(ω0, Δ, Nl, Nr, dl, dr, chainparamsl=[fill(1.0,N),fill(1.0,N-1), 1.0], chainparamsr=chainparamsl; tree=false)
    u = unitmat(2)

    cl = only(chainparamsl[3])
    cr = only(chainparamsr[3])

    Hs = (ω0/2)*sz + Δ*sx

    M=zeros(3,3,2,2)
    M[1,:,:,:] = up(Hs, cr*sx, u)
    M[:,1,:,:] = up(Hs, cl*sx, u)

    chain1 = hbathchain(Nl, dl, chainparamsl; coupletox=true, reverse=true)
    chain2 = hbathchain(Nr, dr, chainparamsr; coupletox=true)
    return Any[chain1..., M, chain2...]
end

"""
    chaincoeffs_ohmic(N, α, s; ωc=1, soft=false)

Generate chain coefficients ``[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]`` for an Harmonic bath at zero temperature with a power
law spectral density given by:

soft cutoff: ``J(ω) = 2παω_c (\\frac{ω}{ω_c})^s \\exp(-ω/ω_c)`` \n
hard cutoff: ``J(ω) = 2παω_c (\\frac{ω}{ω_c})^s θ(ω-ω_c)``

The coefficients parameterise the chain Hamiltonian

``
H = H_S + c_0 A_S⊗B_0+\\sum_{i=0}^{N-1}t_i (b_{i+1}^\\dagger b_i +h.c.) + \\sum_{i=0}^{N-1} ϵ_ib_i^\\dagger b_i
``

which is unitarily equivalent (before the truncation to `N` sites) to

``
H = H_S + A_S⊗\\int_0^∞dω\\sqrt{\\frac{J(ω)}{π}}B_ω + \\int_0^∞dωωb_ω^\\dagger b_ω 
``

"""
function chaincoeffs_ohmic(nummodes, α, s; ωc=1, soft=false)
    if soft
        c0 = ωc*sqrt(2*α*gamma(s+1))
        e = [ωc*(2n + 1 + s) for n in 0:(nummodes-1)]
        t = [ωc*sqrt((n + 1)*(n + s + 1)) for n in 0:(nummodes-2)]
        return [e, t, c0]
    else    
        c0 = sqrt(2α/(s+1))*ωc
        e = [(ωc/2)*(1 + (s^2)/((s+2n)*(2+s+2n))) for n in 0:(nummodes-1)]
        t = [ωc*(1+n)*(1+s+n)/((s+2+2n)*(s+3+2n))*sqrt((3+s+2n)/(1+s+2n)) for n in 0:(nummodes-2)]
        return [e, t, c0]
    end
end

#Deprecated#
function getchaincoeffs(nummodes, α, s, beta, ωc=1)
    matlabdir = ENV["MATDIR"]
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
##

function readchaincoeffs(fdir, params...)
    n = length(params)
    dat = h5open(fdir, "r") do fid
        g = fid
        for i=1:n
            vars = keys(g)
            par = params[i]
            if typeof(par) <: Number
                ind = only(findall(x->x==par, parse.(Float64, vars)))
                g = g[vars[ind]]
            elseif typeof(par) <: String
                g = g[par]
            else
                throw(ArgumentError("parameter in position $i not valid"))
            end
        end
        return [read(g["e"]), read(g["t"]), read(g["c"])]
    end
    return dat
end

function ibmmpo(ω0, d, N, chainparams; tree=false)
    u = unitmat(2)
    
    c0 = only(chainparams[3])

    Hs = (ω0/2)*sz
#    Hs = (ω0/2)*(sz+u)

    M=zeros(1,3,2,2)
    M[1, :, :, :] = up(Hs, c0*sz, u)
#    M[1, :, :, :] = up(Hs, c0*(sz+u), c0*(sz+u), u)

    chain = hbathchain(N, d, chainparams; tree=tree, coupletox=true)
    if tree
        H = TreeNetwork(Vector{AbstractArray}([M]))
        addtree!(H, 1, chain)
        return H
    else
        return Any[M, chain...]
    end
end

function tunnelingmpo(ϵ, delta, α, s, β, d::Int, nummodes::Int; tree=false, ωc=1)
    cps = chaincoeffs_ohmic(nummodes, α, s, β; ωc=ωc)
    λ = 2*α*ωc/s + delta
    
    u = unitmat(2)
    
    c0 = only(cps[3])

    Hs = (ϵ/2)*sz + λ*(u + sx)/2

    M=zeros(1,4,2,2)
    M[1, :, :, :] = up(Hs, c0*(u + sx)/2, c0*(u + sx)/2, u)

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

    M=zeros(D, D, d, d)
    M[D, :, :, :] = A == Ad ? up(h0, A, u) : up(h0, A, Ad, u)
    M[:, 1, :, :] = A == Ad ? dn(h0, A, u) : dn(h0, A, Ad, u)

    return Any[M[D:D,:,:,:], fill(M, N-2)..., M[:,1:1,:,:]]
end

function nearestneighbourmpo(tree_::Tree, h0, A, Ad = A')
    size(h0) == size(A) || error("physical dimensions don't match")
    size(h0) == size(Ad) || error("physical dimensions don't match")
    d = size(h0)[1]

    tree = deepcopy(tree_)

    D = A == Ad ? 3 : 4

    u=unitmat(d)

    N = length(tree)
    hn = findheadnode(tree)

    Ms = Vector{Any}(undef, N)

    for (id, nd) in enumerate(tree)
        nc = length(nd.children)
        dims = [fill(D, nc+1)..., d, d]
        M = zeros(dims...)
        M[:, fill(1, nc)..., :, :] = A == Ad ? dn(h0, A, u) : dn(h0, A, Ad, u)
        for i in 1:nc
            M[D, fill(1, i-1)..., :, fill(1, nc-i)..., :, :] = A == Ad ? up(h0, A, u) : up(h0, A, Ad, u)
        end
        Ms[id] = M
    end
    nc = length(tree[hn].children)
    Ms[hn] = Ms[hn][D:D, fill(:,nc+2)...]
    return TreeNetwork(tree, Ms)
end
