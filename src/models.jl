up(ops...) = permutedims(cat(ops...; dims=3), [3,1,2])
dn(ops...) = permutedims(cat(reverse(ops)...; dims=3), [3,1,2])

"""
    xyzmpo(N::Int; Jx=1.0, Jy=Jx, Jz=Jx, hx=0., hz=0.)

Generate MPO for the `N`-spin XYZ model with external field ``\\vec{h}=(h_x, 0, h_z)``, , defined by the Hamiltonian

``
H = \\sum_{n=1}^{N-1} -J_x σ_x^{n} σ_x^{n+1} - J_y σ_y^{n} σ_y^{n+1} - J_z σ_z^{n} σ_z^{n+1} + \\sum_{n=1}^{N}(- h_x σ_x^{n} - h_z σ_z^{n})  
``

with ``σ_x^{n}, σ_y^{n}, σ_z^{n}`` the Pauli spin-1/2 matrices of the ``n^\\text{th}`` site.

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

Generate MPO for the `N`-spin 1D Ising model with external field ``\\vec{h} = (0,0,h)``, defined by the Hamiltonian

``
H = \\sum_{n=1}^{N-1} -J_x σ_x^{n} σ_x^{n+1} + \\sum_{n=1}^{N}(- h_z σ_z^{n})  
``

with ``σ_x^{n}, σ_y^{n}, σ_z^{n}`` the Pauli spin-1/2 matrices of the ``n^\\text{th}`` site.

"""
isingmpo(N::Int; J=1.0, h=1.0) = xyzmpo(N; Jx=J, Jy=0., Jz=0., hz=h, hx=0.)
"""
    heisenbergmpo(N::Int, J=1.0) = xyzmpo(N; Jx=J)

Generate MPO for the `N`-spin Heisenberg XXX model, defined by the Hamiltonian

``
H = \\sum_{n=1}^{N-1} -J σ_x^{n} σ_x^{n+1} - J σ_y^{n} σ_y^{n+1} - J σ_z^{n} σ_z^{n+1}   
``

with ``σ_x^{n}, σ_y^{n}, σ_z^{n}`` the Pauli spin-1/2 matrices of the ``n^\\text{th}`` site.


"""
heisenbergmpo(N::Int, J=1.0) = xyzmpo(N; Jx=J)
"""
    xxzmpo(N::Int, Δ = 1.0, J=1.0) = xyzmpo(N; Jx=J, Jy=J, Jz=J*Δ)

Generate MPO for the `N`-spin XXZ model, defined by the Hamiltonian

``
H = \\sum_{n=1}^{N-1} -J σ_x^{n} σ_x^{n+1} - J σ_y^{n} σ_y^{n+1} - \\Delta J σ_z^{n} σ_z^{n+1}   
``

with ``σ_x^{n}, σ_y^{n}, σ_z^{n}`` the Pauli spin-1/2 matrices of the ``n^\\text{th}`` site.

"""
xxzmpo(N::Int, Δ = 1.0, J=1.0) = xyzmpo(N; Jx=J, Jy=J, Jz=J*Δ)

"""
    longrange_xyzmpo(N::Int, α::Float64=0.; Jx=1.0, Jy=Jx, Jz=Jx, hx=0., hz=0.)

Gennerate MPO for the `N`-spin long-range XYZ model with external field ``\\vec{h}=(h_x, 0, h_z)``, , defined by the Hamiltonian


"""
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

"""
    longrange_isingmpo(N::Int, α::Float64=0.; J=1.0, h=1.0) = longrange_xyzmpo(N, α; Jx=J, Jy=0., Jz=0., hz=h, hx=0.)

"""
longrange_isingmpo(N::Int, α::Float64=0.; J=1.0, h=1.0) = longrange_xyzmpo(N, α; Jx=J, Jy=0., Jz=0., hz=h, hx=0.)

"""
    spinchainmpo(N::Int; J=1.0, hz=1.0, hx=0.0, i=div(N,2))


"""
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

"""
    tightbindingmpo(N::Int, d::Int; J=1.0, e=1.0)



"""
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

Generate MPO representing a tight-binding chain of `N` oscillators with `d` Fock states each. Chain parameters are supplied in the standard form: `chainparams` ``=[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]``. The output does not itself represent a complete MPO but will possess an end which is *open* and should be attached to another tensor site, usually representing the *system*.

# Arguments

* `reverse`: If `reverse=true` create a chain were the last (i.e. Nth) site is the site which couples to the system
* `coupletox`: Used to choose the form of the system coupling. `coupletox=true` gives a non-number conserving coupling of the form ``H_{\\text{I}}= A_{\\text{S}}(b_{0}^\\dagger + b_0)`` where ``A_{\\text{S}}`` is a system operator, while `coupletox=false` gives the number-converving coupling ``H_{\\text{I}}=(A_{\\text{S}} b_{0}^\\dagger + A_{\\text{S}}^\\dagger b_0)``
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
H = \\frac{ω_0}{2}σ_z + Δσ_x + c_0σ_x(b_0^\\dagger+b_0) + \\sum_{i=0}^{N-1} t_i (b_{i+1}^\\dagger b_i +h.c.) + \\sum_{i=0}^{N} ϵ_ib_i^\\dagger b_i
``.

The spin is on site 1 of the MPS and the bath modes are to the right.

This Hamiltonain is unitarily equivalent (before the truncation to `N` sites) to the spin-boson Hamiltonian defined by

``
H =  \\frac{ω_0}{2}σ_z + Δσ_x + σ_x\\int_0^∞ dω\\sqrt{J(ω)}(b_ω^\\dagger+b_ω) + \\int_0^∞ dω ωb_ω^\\dagger b_ω
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

"""
    twobathspinmpo(ω0, Δ, Nl, Nr, dl, dr, chainparamsl=[fill(1.0,N),fill(1.0,N-1), 1.0], chainparamsr=chainparamsl; tree=false)

Generate MPO for a spin-1/2 coupled to two chains of harmonic oscillators, defined by the Hamiltonian

``
H = \\frac{ω_0}{2}σ_z + Δσ_x + c_0^rσ_x(b_0^\\dagger+b_0) + \\sum_{i=0}^{N_r-1} t_i^r (b_{i+1}^\\dagger b_i +h.c.) + \\sum_{i=0}^{N_r} ϵ_i^rb_i^\\dagger b_i + c_0^lσ_x(d_0^\\dagger+d_0) + \\sum_{i=0}^{N_l-1} t_i^l (d_{i+1}^\\dagger d_i +h.c.) + \\sum_{i=0}^{N_l} ϵ_i^l d_i^\\dagger d_i
``.

The spin is on site ``N_l + 1`` of the MPS, surrounded by the left chain modes and the right chain modes.

This Hamiltonain is unitarily equivalent (before the truncation to `N` sites) to the spin-boson Hamiltonian defined by

``
H =  \\frac{ω_0}{2}σ_z + Δσ_x + σ_x\\int_0^∞ dω\\sqrt{\\frac{J(ω)}{π}}(b_ω^\\dagger+b_ω) + \\int_0^∞ dω ωb_ω^\\dagger b_ωi + σ_x\\int_0^∞ dω\\sqrt{\\frac{J^l(ω)}{π}}(d_ω^\\dagger+d_ω) + \\int_0^∞ dω ωd_ω^\\dagger d_ω
``.

The chain parameters, supplied by `chainparams`=``[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]``, can be chosen to represent any arbitrary spectral density ``J(ω)`` at any temperature. The two chains can have a different spectral density.

"""
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

soft cutoff: ``J(ω) = 2αω_c (\\frac{ω}{ω_c})^s \\exp(-ω/ω_c)`` \n
hard cutoff: ``J(ω) = 2αω_c (\\frac{ω}{ω_c})^s θ(ω-ω_c)``

The coefficients parameterise the chain Hamiltonian

``
H = H_S + c_0 A_S⊗B_0+\\sum_{i=0}^{N-1}t_i (b_{i+1}^\\dagger b_i +h.c.) + \\sum_{i=0}^{N} ϵ_ib_i^\\dagger b_i
``

which is unitarily equivalent (before the truncation to `N` sites) to

``
H = H_S + A_S⊗\\int_0^∞dω\\sqrt{J(ω)}B_ω + \\int_0^∞dωωb_ω^\\dagger b_ω 
``.

"""
function chaincoeffs_ohmic(nummodes, α, s; ωc=1, soft=false)
    if soft
        c0 = ωc*sqrt(2*α*gamma(s+1))
        e = [ωc*(2n + 1 + s) for n in 0:(nummodes-1)]
        t = [ωc*sqrt((n + 1)*(n + s + 1)) for n in 0:(nummodes-2)]
        return [e, t, c0]
    else
        if s==0
            c0 = sqrt(2α)*ωc
            e = fill(0,nummodes)
            t = [ωc*(n+1)/(2n+1) for n in 0:(nummodes-2)]
        else
            c0 = sqrt(2α/(s+1))*ωc
            e = [(ωc/2)*(1 + (s^2)/((s+2n)*(2+s+2n))) for n in 0:(nummodes-1)]
            t = [ωc*(1+n)*(1+s+n)/((s+2+2n)*(s+3+2n))*sqrt((3+s+2n)/(1+s+2n)) for n in 0:(nummodes-2)]
        end
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

"""
    readchaincoeffs(fdir, params...)




"""
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

"""
    ibmmpo(ω0, d, N, chainparams; tree=false)

Generate MPO for a spin-1/2 coupled to a chain of harmonic oscillators with the interacting boson model (IBM), defined by the Hamiltonian

``
H = \\frac{ω_0}{2}σ_z +  c_0σ_z(b_0^\\dagger+b_0) + \\sum_{i=0}^{N-1} t_i (b_{i+1}^\\dagger b_i +h.c.) + \\sum_{i=0}^{N} ϵ_ib_i^\\dagger b_i
``.

The spin is on site 1 of the MPS and the bath modes are to the right.

This Hamiltonain is unitarily equivalent (before the truncation to `N` sites) to the spin-boson Hamiltonian defined by

``
H =  \\frac{ω_0}{2}σ_z + σ_z\\int_0^∞ dω\\sqrt{\\frac{J(ω)}{π}}(b_ω^\\dagger+b_ω) + \\int_0^∞ dω ωb_ω^\\dagger b_ω
``.

The chain parameters, supplied by `chainparams`=``[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]``, can be chosen to represent any arbitrary spectral density ``J(ω)`` at any temperature.

"""
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

"""
    tunnelingmpo(ϵ, delta, α, s, β, d::Int, nummodes::Int; tree=false, ωc=1)





"""
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

"""
    nearestneighbourmpo(N::Int, h0, A, Ad = A')





"""
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

"""
    nearestneighbourmpo(tree_::Tree, h0, A, Ad = A')





"""
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

"""
    puredephasingmpo(ΔE, dchain, Nchain, chainparams; tree=false)

    Generate MPO for a pure dephasing model, defined by the Hamiltonian
    ``H = \\frac{ΔE}{2} σ_z +  \\frac{σ_z}{2} c_0 (b_0^\\dagger + b_0) + \\sum_{i=0}^{N-1} t_i (b_{i+1}^\\dagger b_i +h.c.) + \\sum_{i=0}^{N-1} ϵ_i b_i^\\dagger b_i  ``

    The spin is on site 1 of the MPS and the bath modes are to the right.

    ### Arguments
    * `ΔE::Real`: energy splitting of the spin
    * `dchain::Int`: physical dimension of the chain sites truncated Hilbert spaces
    * `Nchain::Int`: number of sites in the chain
    * `chainparams::Array{Real,1}`: chain parameters for the bath chain. The chain parameters are given in the standard form: `chainparams` ``=[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]``.
    * `tree::Bool`: if true, return a `TreeNetwork` object, otherwise return a vector of MPO tensors
"""

function puredephasingmpo(ΔE, dchain, Nchain, chainparams; tree=false)
    u = unitmat(2)

    cparams = only(chainparams[3])

    Hs = (ΔE/2)*sz

    M=zeros(1,3,2,2)
    M[1, :, :, :] = up(Hs, (cparams/2)*sz, u)

    chain = hbathchain(Nchain, dchain, chainparams; tree=false, reverse=false, coupletox=true)
    return Any[M, chain...]
end

"""
    tightbinding_mpo(N, ϵd, chainparams1, chainparams2)

    Generate MPO for a tight-binding chain of N fermionic sites with a single impurity site (fermionic as well) 
    of energy ϵd at the center. The impurity is coupled to two leads, each described by a set of chain parameters.
    The interactions are nearest-neighbour, with the first N/2-1 sites corresponding to the first lead,
    the Nth site corresponding to the impurity, and the rest of the sites corresponding to the second
    lead.

    * `N::Int`: number of sites in the chain
    * `ϵd::Real`: energy of the impurity site at the center, as Ed - μ, where μ is the chemical potential
    * `chainparams1::Array{Real,1}`: chain parameters for the first lead
    * `chainparams2::Array{Real,1}`: chain parameters for the second lead

    The chain parameters are given in the standard form: `chainparams` ``=[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]``.
"""

function tightbinding_mpo(N, ϵd, chainparams1, chainparams2)

    e1 = chainparams1[1]
    t1 = chainparams1[2]
    c1 = chainparams1[3]

    e2 = chainparams2[1]
    t2 = chainparams2[2]
    c2 = chainparams2[3]

    d = 2 # physical dimension for fermionic Hilbert space
    D = 4 # bond dimensions of the MPO
    u = unitmat(d)
    n = numb(d)
    b = anih(d)
    bd = crea(d)

    W = Any[]

    # SECOND CHAIN: filled modes

    # the initial matrix for the first filled site:
    Min = zeros(ComplexF64, 1, D, d, d)
    Min[1,1,:,:] = u
    Min[1,D,:,:] = e2[N]*n
    Min[1,2,:,:] = t2[N-1]*bd
    Min[1,3,:,:] = t2[N-1]*b
    push!(W, Min)

    if N!=1

        M2 = zeros(ComplexF64, D, D, d, d)
        for i in 1:(N-2)
            M2[1,1,:,:] = u
            M2[D,D,:,:] = u
            M2[1,D,:,:] = e2[N-i]*n
            M2[1,2,:,:] = t2[N-i-1]*bd
            M2[3,D,:,:] = bd
            M2[1,3,:,:] = t2[N-i-1]*b
            M2[2,D,:,:] = b
            push!(W, M2)
        end
    end

    # MPO for the filled site of the second chain coupled to the system:
    Mcoup2 = zeros(ComplexF64, D, D, d, d)
    Mcoup2[1,1,:,:] = u
    Mcoup2[D,D,:,:] = u
    Mcoup2[1,D,:,:] = e2[1]*n
    Mcoup2[1,2,:,:] = c2*bd
    Mcoup2[3,D,:,:] = bd
    Mcoup2[1,3,:,:] = c2*b
    Mcoup2[2,D,:,:] = b
    push!(W, Mcoup2)


    #system MPO
    Md = zeros(ComplexF64, D, D, d, d)
    Md[1,1,:,:] = u
    Md[D,D,:,:] = u
    Md[1,D,:,:] = ϵd*n
    Md[1,2,:,:] = bd
    Md[3,D,:,:] = bd
    Md[1,3,:,:] = b
    Md[2,D,:,:] = b
    push!(W, Md)

    # FIRST CHAIN: empty modes

    # MPO for the empty site of the first chain coupled to the system:
    Mcoup1 = zeros(ComplexF64, D, D, d, d)
    Mcoup1[1,1,:,:] = u
    Mcoup1[D,D,:,:] = u
    Mcoup1[1,D,:,:] = e1[1]*n
    Mcoup1[1,2,:,:] = bd
    Mcoup1[3,D,:,:] = c1*bd
    Mcoup1[1,3,:,:] = b
    Mcoup1[2,D,:,:] = c1*b
    push!(W, Mcoup1)


    if N!=1
        M1 = zeros(ComplexF64, D, D, d, d)
        for i in 2:(N-1)
            M1[1,1,:,:] = u
            M1[D,D,:,:] = u
            M1[1,D,:,:] = e1[i]*n
            M1[1,2,:,:] = bd
            M1[3,D,:,:] = t1[i-1]*bd
            M1[1,3,:,:] = b
            M1[2,D,:,:] = t1[i-1]*b
            push!(W, M1)
        end
    end
    # the final matrix for the last empty site:
    Mfin = zeros(ComplexF64, D, 1, d, d)
    Mfin[1,1,:,:] = e1[N]*n
    Mfin[D,1,:,:] = u
    Mfin[3,1,:,:] = t1[N-1]*bd
    Mfin[2,1,:,:] = t1[N-1]*b
    push!(W, Mfin)
    return W
end

"""
    interleaved_tightbinding_mpo(N, ϵd, chainparams1, chainparams2)

    Generate MPO for a tight-binding chain of N fermionic sites with a single impurity site (fermionic as well)
    of energy ϵd. The impurity is coupled to two leads, each described by a set of chain parameters. 
    The interactions are next-nearest-neighbour, with the first site corresponding to the impurity, and the
    two chains organised in an interleaved fashion.

    # Arguments

    * `N::Int`: number of sites in the chain
    * `ϵd::Real`: energy of the impurity site at the first site, as Ed - μ, where μ is the chemical potential
    * chainparams1::Array{Real,1}: chain parameters for the first lead
    * chainparams2::Array{Real,1}: chain parameters for the second lead

    The chain parameters are given in the standard form: `chainparams` ``=[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]``.
"""

function interleaved_tightbinding_mpo(N, ϵd, chainparams1, chainparams2)

    e1 = chainparams1[1]
    t1 = chainparams1[2]
    c1 = chainparams1[3]
#    c1 = 0.

    e2 = chainparams2[1]
    t2 = chainparams2[2]
    c2 = chainparams2[3]
#    c2 = 0.

    d = 2 # physical dimension for fermionic Hilbert space
    D = 6 # bond dimensions of the MPO
    u = unitmat(d)
    n = numb(d)
    b = anih(d)
    bd = crea(d)
    F = u .- 2.0.*n

    W = Any[]


    # the initial matrix for the system:
    Min = zeros(ComplexF64, 1, D, d, d)
    Min[1,1,:,:] = u
    Min[1,D,:,:] = ϵd*n
    Min[1,2,:,:] = b
    Min[1,3,:,:] = bd

    push!(W, Min)

    if N!=1
        # SECOND CHAIN: filled modes
        M2in = zeros(ComplexF64, D, D, d, d)
        # FIRST CHAIN: empty modes
        M1in = zeros(ComplexF64, D, D, d, d)
        # first site, COUPLED TO THE SYSTEM
        M2in[1,1,:,:] = u
        M2in[1,2,:,:] = b
        M2in[1,3,:,:] = bd
        M2in[1,D,:,:] = e2[1]*n
        M2in[2,4,:,:] = F
        M2in[3,5,:,:] = F
        M2in[2,D,:,:] = c2*bd
        M2in[3,D,:,:] = c2*b
        M2in[D,D,:,:] = u
        push!(W, M2in)

        M1in[1,1,:,:] = u
        M1in[1,2,:,:] = b
        M1in[1,3,:,:] = bd
        M1in[1,D,:,:] = e1[1]*n
        M1in[2,4,:,:] = F
        M1in[3,5,:,:] = F
        M1in[4,D,:,:] = c1*bd
        M1in[5,D,:,:] = c1*b
        M1in[D,D,:,:] = u
        push!(W, M1in)

        for i in 2:(N-1)
            # SECOND CHAIN: filled modes
            M2 = zeros(ComplexF64, D, D, d, d)
            # FIRST CHAIN: empty modes
            M1 = zeros(ComplexF64, D, D, d, d)

            M2[1,1,:,:] = u
            M2[1,2,:,:] = b
            M2[1,3,:,:] = bd
            M2[1,D,:,:] = e2[i]*n
            M2[2,4,:,:] = F
            M2[3,5,:,:] = F
            M2[4,D,:,:] = t2[i-1]*bd
            M2[5,D,:,:] = t2[i-1]*b
            M2[D,D,:,:] = u
            push!(W, M2)

            M1[1,1,:,:] = u
            M1[1,2,:,:] = b
            M1[1,3,:,:] = bd
            M1[1,D,:,:] = e1[i]*n
            M1[2,4,:,:] = F
            M1[3,5,:,:] = F
            M1[4,D,:,:] = t1[i-1]*bd
            M1[5,D,:,:] = t1[i-1]*b
            M1[D,D,:,:] = u
            push!(W, M1)

        end
    end
    # The final matrix for the SECOND chain:
    M2fin = zeros(ComplexF64, D, D, d, d)
    M2fin[1,1,:,:] = u
    M2fin[1,2,:,:] = b
    M2fin[1,3,:,:] = bd
    M2fin[1,D,:,:] = e2[N]*n
    M2fin[2,4,:,:] = F
    M2fin[3,5,:,:] = F
    M2fin[4,D,:,:] = t2[N-1]*bd
    M2fin[5,D,:,:] = t2[N-1]*b
    M2fin[D,D,:,:] = u
    push!(W, M2fin)

    # The final matrix for the FIRST chain:
    M1fin = zeros(ComplexF64, D, 1, d, d)
    M1fin[1,1,:,:] = e1[N]*n
    M1fin[4,1,:,:] = t1[N-1]*bd
    M1fin[5,1,:,:] = t1[N-1]*b
    M1fin[D,1,:,:] = u
    push!(W, M1fin)
    return W
end

"""
    correlatedenvironmentmpo(R::Vector, Nm::Int, d::Int; chainparams, fnamecc::String, s=1, α=1, ωc=1, c_phonon=1, β="inf", issoft=false)

Generate a MPO for a one-dimensional bosonic bath spatially correlated to a multi-component system 

``
H_B + H_int = \\int_{-∞}^{+∞} dk ω_k b_k^\\dagger b_k + ∑_j \\int_{-∞}^{+∞}dk \\sqrt{J(k)}(A_j b_k e^{i k R_j} + h.c.)
``.

The interactions between the system and the chain-mapped bath are long range, i.e. each site interacts with all the chain modes. The spectral density is assumed to be Ohmic ``J(ω) = 2αωc(ω/ωc)^s``.

# Arguments

* `R`: List of system's components positions
* `Nm`: Number of chain modes. The actual number of mode will be doubled to account for the left and right moving excitations.
* `d`: Local Hilbert space dimension of the bath modes
* `chainparams`: chain parameters, of the form `chainparams`=``[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]``, can be chosen to represent any arbitrary spectral density ``J(ω)`` at any temperature.
* `fnamecc`: Path to a file containing pre-computed long-range coupling coefficient. If not provided, the coupling coefficients will be computed and stored.
* `s`: Ohmicity
* `α`: Kondo parameter
* `ωc`: Bath cut-off frequency
* `c_phonon`: Speed of sound in the bath
* `β`: Inverse temperature 
* `issoft`: Is the cut-off of the Ohmic SD soft or hard?
"""
function correlatedenvironmentmpo(R::Vector, Nm::Int, d::Int; chainparams, fnamecc::String, s=1, α=1, ωc=1, c_phonon=1, β="inf", issoft=false)

    function polybeta(t::Float64, n::Int, a::Array, b::Array, temp::Array)
    """
        polybeta recursively constructs the polynomials used to compute the coupling coefficients given the coefficients a and b
        This function is useful when working at finite temperature (β != inf)
    """
        if n==-1
            return 0
        elseif n==0
            if length(temp)>=2
                temp[2] = 1
            else
                push!(temp,1)
            end

            return 1
        elseif n==1
            pn = (t - a[n])
            if length(temp) == 2
                push!(temp,pn)
            elseif length(temp) == 1
                push!(temp, 1)
                push!(temp,pn)
            end

            return pn
        else
            if length(temp)<n+1 && temp[1] == t
                pn = (t - a[n])*polybeta(t,n-1,a,b,temp) - b[n-1]*polybeta(t,n-2,a,b,temp) #P_{n}(t) = (t-a_{n-1})P_{n-1} - b_{n-1}P_{n-2}
                push!(temp, pn)

                return pn
            elseif length(temp)<n+1 && temp[1] != t
                temp = [t]
                pn = (t - a[n])*polybeta(t,n-1,a,b,temp) - b[n-1]*polybeta(t,n-2,a,b,temp) #P_{n}(t) = (t-a_{n-1})P_{n-1} - b_{n-1}P_{n-2}
                push!(temp, pn)

                return pn
            elseif length(temp) == n+1 && temp[1] == t
                pn = (t - a[n])*temp[n+1] - b[n-1]*temp[n]
                push!(temp,pn)

                return pn
            elseif length(temp) == n+1 && temp[1] != t
                temp = [t]
                pn = (t - a[n])*polybeta(t,n-1,a,b,temp) - b[n-1]*polybeta(t,n-2,a,b,temp) #P_{n}(t) = (t-a_{n-1})P_{n-1} - b_{n-1}P_{n-2}
                push!(temp, pn)

                return pn
            elseif length(temp) > n+1 && temp[1] == t
                pn = temp[n+2]

                return pn
            else
                temp = [t]
                pn = (t - a[n])*polybeta(t,n-1,a,b,temp) - b[n-1]*polybeta(t,n-2,a,b,temp) #P_{n}(t) = (t-a_{n-1})P_{n-1} - b_{n-1}P_{n-2}
                push!(temp, pn)

                return pn
            end
        end
    end
    
    function SDOhmic(t)
    """
        Bath Ohmic Spectral Density for zero temperature chain mapping of the bath
    """
        if t==0
            return 0
        elseif t>-1 && t<1
            return 2*α*abs(t)*ωc
        elseif abs(t)==1
            return 2*α*ωc
        else
            return 0
        end
    end
    
        
    function SDTOhmic(t)
    """
        Bath Ohmic Spectral Density after the finite temperature chain mapping of the bath
    """
        if t==0
            return 2*α/β
        elseif t>-1 && t<1
            return α*t*ωc*(1+coth(β*t*ωc*0.5))
        elseif abs(t)==1
            return α*t*ωc*(1+coth(β*t*ωc*0.5))
        else
            return 0
        end
    end

    a_chain = chainparams[1]
    b_chain = chainparams[2].^2

    Norm = zeros(Nm)

    function γ(x::Int, n::Int, issoft::Bool; β="inf", temp=[1.])
    """
        Definition of the coupling coefficient between the site x and the mode n for a Ohmic spectral density with a hard cut-off (Jacobi Polynomials) or a soft cut-off Laguerre Polynomials
    """
        if β=="inf"
            if issoft==true
                polynomial0(t) = sf_laguerre_n(n,s,t)*exp(-im*t*R[x]*ωc/c_phonon)*t^s*exp(-s)
                return sqrt(2*α*gamma(n+s + 1)/gamma(n+1))*ωc*quadgk(polynomial0, 0, 1)[1]
            else
                polynomial(t) = jacobi(2*t-1,n-1, 0, s)*exp(-im*t*R[x]*ωc/c_phonon)*t^s
		        return sqrt(2*α*(2*(n-1) + s + 1))*ωc*quadgk(polynomial, 0, 1)[1]
            end
        elseif β!="inf"
           polynomial1(t) = polybeta(t,n-1,a_chain,b_chain,[t])
           integrand(t) = polynomial1(t)*exp(im*t*R[x]*ωc/c_phonon)*SDTOhmic(t)
           N2(t) = polynomial1(t)^2*SDTOhmic(t)
           if Norm[n]==0
               Norm[n] = sqrt(quadgk(N2,-1,1)[1])
           end
           return (ωc/Norm[n])*quadgk(integrand, -1, 1)[1]

        end
    end

    # Construction of the MPO
    W = Any[] # list of the MPO's tensors

    ### Construction of the bosonic bath MPO ###
    e = chainparams[1] # list of the energy of each modes
    t = chainparams[2] # list of the hopping energy between modes
    u = unitmat(d)
    # modes creation, anihilation and number operators
    bd = crea(d)
    b = anih(d)
    n = numb(d)

    N = length(R) # Number of system's sites

    coupling_stored = zeros(ComplexF64,N,Nm) # just need NxNm because the two chains have the same coupling coeff up to complex conjugation
    arestored = 1 # are the coupling coefficient stored
    Nstored, Nmstored = N, Nm # number of stored coeff.
    couplinglist = []
    try
        couplinglist = readdlm(fnamecc,',',ComplexF64,'\n')
    catch stored_error
        if isa(stored_error, ArgumentError)
            print("Coupling coefficient not found. They will be computed and stored.\n")
            arestored = 0
        end
    end
    if couplinglist != []
        Nstored, Nmstored = size(couplinglist)
        if Nmstored>=Nm && Nstored>=N
            coupling_stored = couplinglist[1:N,1:Nm]
        else
            coupling_stored[1:min(N,Nstored),1:min(Nm,Nmstored)] = couplinglist[1:min(N,Nstored),1:min(Nm,Nmstored)]
            print("Less coupling coefficient stored than needed. Available ones will be used and missing one will be computed and stored.\n")
        end
    end

    ## First chain MPO
    D = 2*(N + 2) #Bond dimension
    M = zeros(ComplexF64,D-2, D, d, d)
    M[1,1,:,:] = M[D-2,D,:,:] = u
    M[1,D,:,:] = e[1]*n
    i = 2
    M[1,i,:,:] = t[1]*bd
    M[1,i+1,:,:] = t[1]*b

    a = 0 #site counter
    while i<D-2
        a+=1
        if arestored==1
            couplingcoeff = coupling_stored[a,1]
        else
            couplingcoeff = γ(a,1,issoft, β=β)
            coupling_stored[a,1] = couplingcoeff
        end
        couplingcoeff = couplingcoeff
        M[i,D,:,:] = couplingcoeff*b
        M[i,i+2,:,:] = u
        i+=1
        M[i,D,:,:] = conj(couplingcoeff)*bd
        M[i,i+2,:,:] = u
        i+=1
    end
    M = reshape(M,D-2,D,d,d)
    push!(W, M)

    for m = 2:Nm-1
        D = 2*(N + 2)
        M = zeros(ComplexF64,D, D, d, d)
        M[1,1,:,:] = M[D,D,:,:] = u
        M[1,D,:,:] = e[m]*n
        i = 2
        M[1,i,:,:] = t[m]*bd
        M[1,i+1,:,:] = t[m]*b
        M[i,D,:,:] = b
        M[i+1,D,:,:] = bd
        i += 2

        a = 0 #site counter
        while i<D
            a+=1
            if arestored==1 && m<=Nmstored && a<=Nstored
                couplingcoeff = coupling_stored[a,m]
            else
                couplingcoeff = γ(a, m, issoft,β=β)
                coupling_stored[a,m] = couplingcoeff
            end
            M[i,D,:,:] = couplingcoeff*b
            if i<D-1
                M[i,i,:,:] = u
            end
            i+=1
            M[i,D,:,:] = conj(couplingcoeff)*bd
            if i<D
                M[i,i,:,:] = u
            end
            i+=1
        end

        M = reshape(M,D,D,d,d)
        push!(W, M)
    end

    # Last Mode of the First Chain
    D = 2*(N + 2)
    M = zeros(ComplexF64,D, D, d, d)
    M[1,1,:,:] = M[D,D,:,:] = u
    M[1,D,:,:] = e[Nm]*n
    i = 2
    M[i,D,:,:] = b
    M[i+1,D,:,:] = bd
    i += 2

    a = 0 #site counter
    while i<D
        a+=1
        if arestored==1 && Nm<=Nmstored && a<=Nstored
            couplingcoeff = coupling_stored[a,Nm]
        else
            couplingcoeff = γ(a, Nm, issoft,β=β)
            coupling_stored[a,Nm] = couplingcoeff
        end
        M[i,D,:,:] = couplingcoeff*b
        if i<D-1
            M[i,i,:,:] = u
        end
        i+=1
        M[i,D,:,:] = conj(couplingcoeff)*bd
        if i<D
            M[i,i,:,:] = u
        end
        i+=1
    end

    M = reshape(M,D,D,d,d)
    push!(W, M)

    # Second chain
    for m = 1:Nm-1
        D = 2*(N + 2)
        M = zeros(ComplexF64,D, D, d, d)
        M[1,1,:,:] = M[D,D,:,:] = u
        M[1,D,:,:] = e[m]*n
        i = 2
        M[1,i,:,:] = t[m]*bd
        M[1,i+1,:,:] = t[m]*b
        M[i,D,:,:] = b
        M[i+1,D,:,:] = bd
        i += 2

        a = 0 #site counter
        while i<D
            a+=1
            couplingcoeff = coupling_stored[a,m]
            M[i,D,:,:] = couplingcoeff*bd
            if i<D-1
                M[i,i,:,:] = u
            end
            i+=1
            M[i,D,:,:] = conj(couplingcoeff)*b
            if i<D
                M[i,i,:,:] = u
            end
            i+=1
        end

        M = reshape(M,D,D,d,d)
        push!(W, M)
    end

    # Last mode
    WNm = zeros(ComplexF64,D, 1, d, d)
    WNm[1,1,:,:] = e[Nm]*n
    WNm[2,1,:,:] = b
    WNm[3,1,:,:] = bd
    a = 0 #site counter
    i = 4 #row index
    while i<D
        a+=1
        couplingcoeff = coupling_stored[a, Nm]
        WNm[i,1,:,:] = couplingcoeff*bd
        i+=1
        WNm[i,1,:,:] = conj(couplingcoeff)*b
        i+=1
    end
    WNm[D,1,:,:] = u
    WNm = reshape(WNm,D,1,d,d)

    push!(W, WNm)

    if arestored==0 || Nmstored<Nm || Nstored<N
        writedlm(fnamecc, coupling_stored, ',')
    end

    return W
end

"""
    protontransfermpo(ω0e,ω0k,x0e,x0k, Δ, dRC, d, N, chainparams, RCparams, λreorg)

Generate a MPO for a system described in space with a reaction coordinate (RC) tensor. The RC tensor is coupled to a bosonic bath, taking into account the induced reorganization energy. 

``
H_S + H_RC + H_int^{S-RC} = \\omega^0_{e} |e\\rangle \\langle e| + \\omega^0_{k} |k\\rangle \\langle k| + \\Delta (|e\\rangle \\langle k| + |k\\rangle \\langle e|) + \\omega_{RC} (d^{\\dagger}d + \\frac{1}{2}) + g_{e} |e\\rangle \\langle e|( d + d^{\\dagger})+ g_{k} |k \\rangle \\langle k|( d + d^{\\dagger})
``
``
H_B + H_int^{RC-B} = \\int_{-∞}^{+∞} dk ω_k b_k^\\dagger b_k - (d + d^{\\dagger})\\int_0^∞ dω\\sqrt{J(ω)}(b_ω^\\dagger+b_ω) + \\lambda_{reorg}(d + d^{\\dagger})^2
``.
``
\\lambda_{reorg} = \\int \\frac{J(\\omega)}{\\omega}d\\omega
``.


# Arguments

* `ω0e`: enol energy at x=0 
* `ω0k`: keto energy at x=0
* `x0e`: enol equilibrium displacement
* `x0k`: keto equilibrium displacement 
* `Δ`: direct coupling between enol and keto
* `dRC`: fock space of the RC tensor 
* `d`: number of Fock states of the chain modes
* `N`: length of the chain
* `chainparams`: chain parameters, of the form `chainparams`=``[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]``, can be chosen to represent any arbitrary spectral density ``J(ω)`` at any temperature. 
* `RCparams`: RC tensor parameter, of the form `RCparams`=``[ωRC,-g/x]`` 
* `λreorg`: reorganization energy
"""
function protontransfermpo(ω0e,ω0k,x0e,x0k, Δ, dRC, d, N, chainparams, RCparams, λreorg)
    u = unitmat(2)

    b = anih(dRC)
    bd = crea(dRC)
    n = numb(dRC)
    e = RCparams[1]
    uRC = unitmat(dRC)

    spos = [-x0e 0.; 0. -x0k]

    c0 = only(chainparams[3])
    cRC = only(RCparams[2])

    Hs = (ω0e)*[1. 0.; 0. 0.] + (ω0k)*[0. 0.; 0. 1.] + Δ*sx

    M=zeros(1,3,2,2)
    M[1,:,:,:] = up(Hs, cRC*spos, u)

    MRC = zeros(3,3,dRC,dRC)

    MRC[:, 1, :, :] = dn(e[1]*(n.+diagm([1/2 for i=1:dRC]))+(λreorg)*(b+bd)^2,(b+bd), uRC)
    MRC[3, :, :, :] = up(e[1]*(n.+diagm([1/2 for i=1:dRC]))+(λreorg)*(b+bd)^2, -c0*(b+bd), uRC)

    chain = hbathchain(N, d, chainparams; coupletox=true)

    return Any[M,MRC,chain...]
end

