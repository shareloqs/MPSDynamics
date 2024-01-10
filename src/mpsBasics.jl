"""
    orthcentersmps(A)

Compute the orthoganality centres of MPS `A`.

Return value is a list in which each element is the corresponding site tensor of `A` with the
orthoganility centre on that site. Assumes `A` is right normalised.

"""
function orthcentersmps(A::Vector)
    B = deepcopy(A)
    N = length(B)
    for i in 2:N
        AL, C = QR(B[i-1], 2)
        @tensor AC[:] := C[-1,1] * B[i][1,-2,-3]
        B[i] = AC
    end
    return B
end

"""
    normmps(A::Vector; mpsorthog=:None)

Calculate norm of MPS `A`.

Setting `mpsorthog`=`:Right`/`:Left` will calculate the norm assuming right/left canonical form.
Setting `mpsorthog=OC::Int` will cause the norm to be calculated assuming the orthoganility center is on site
`OC`. If mpsorthog is `:None` the norm will be calculated as an MPS-MPS product.

"""
function normmps(A::Vector; mpsorthog=:None)
    N = length(A)
    if mpsorthog==:Right
        @tensor n[a',b'] := A[1][a',c,s]*conj(A[1][b',c,s])
        return sqrt(real(n[1,1]))
    elseif mpsorthog==:Left
        @tensor n[a',b'] := A[N][a,a',s]*conj(A[N][a,b',s])
        return sqrt(real(n[1,1]))
    elseif typeof(mpsorthog)<:Int
        OC = mpsorthog
        @tensor n = A[OC][a,b,s]*conj(A[OC][a,b,s])
        return sqrt(real(n))
    elseif mpsorthog==:None
        ρ = ones(eltype(A[1]), 1, 1)
        for k=1:N
            @tensor ρ[a,b] := ρ[a',b']*A[k][b',b,s]*conj(A[k][a',a,s])
        end
        return sqrt(real(ρ[1,1]))
    end
end

"""
    mpsrightnorm!(A::Vector, jq::Int=1)

Right orthoganalise MPS `A` up to site `jq`.

"""
function mpsrightnorm!(A::Vector, jq::Int=1)
    N = length(A)
    for i=N:-1:jq+1
        aleft, aright, aup = size(A[i])
        C, AR = lq(reshape(A[i], aleft, aright*aup))
        A[i] = reshape(Matrix(AR), aleft, aright, aup)
        @tensor AC[:] := A[i-1][-1,1,-3] * C[1,-2]
        A[i-1] = AC
    end
end

"""
    mpsleftnorm!(A::Vector, jq::Int=length(A))

Left orthoganalise MPS `A` up to site `jq`.

"""
function mpsleftnorm!(A::Vector, jq::Int=length(A))
    for i=1:jq-1
        AL, C = QR(A[i], 2)
        A[i] = AL
        @tensor AC[:] := C[-1,1] * A[i+1][1,-2,-3]
        A[i+1] = AC
    end
end

"""
    mpsmixednorm!(A::Vector, OC::Int)

Put MPS `A` into mixed canonical form with orthogonality centre on site `OC`.

"""                   
function mpsmixednorm!(A::Vector, OC::Int)
    mpsleftnorm!(A, OC)
    mpsrightnorm!(A, OC)
end

"""
    randmps(physdims::Dims{N}, Dmax::Int, T::Type{<:Number} = Float64) where {N}

Construct a random, right-normalised MPS with local Hilbert space dimensions given by `physdims` and max
bond-dimension given by `Dmax`. 

`T` specifies the element type, eg. use `T=ComplexF64` for a complex valued MPS.

"""
function randmps(physdims::Dims{N}, Dmax::Int, T::Type{<:Number} = Float64) where {N}
    bonddims = Vector{Int}(undef, N+1)
    bonddims[1] = 1
    bonddims[N+1] = 1
    Nhalf = div(N,2)
    for n = 2:N
        bonddims[n] = min(Dmax, bonddims[n-1]*physdims[n-1])
    end
    for n = N:-1:1
        bonddims[n] = min(bonddims[n], bonddims[n+1]*physdims[n])
    end

    As = Vector{Any}(undef, N)
    for n = 1:N
        d = physdims[n]
        Dl = bonddims[n]
        Dr = bonddims[n+1]
        As[n] = reshape(randisometry(T, Dl, Dr*d), (Dl, Dr, d))
    end
    return As
end

"""
    randmps(N::Int, d::Int, Dmax::Int, T=Float64)

Construct a random, `N`-site, right-normalised MPS with all local Hilbert space dimensions given by `d`.

"""
randmps(N::Int, d::Int, Dmax::Int, T=Float64) = randmps(ntuple(n->d, N), Dmax, T)

"""
    productstatemps(physdims::Dims, Dmax=1; state=:Vacuum, mpsorthog=:Right)

Return an MPS representing a product state with local Hilbert space dimensions given by `physdims`.

By default all bond-dimensions will be 1 since the state is a product state. However, to
embed the product state in a manifold of greater bond-dimension, `Dmax` can be set accordingly.

The indvidual states of the MPS sites can be provdided by setting `state` to a list of
column vectors. Setting `state=:Vacuum` will produce an MPS in the vacuum state (where the
state of each site is represented by a column vector with a 1 in the first row and zeros
elsewhere). Setting `state=:FullOccupy` will produce an MPS in which each site is fully
occupied (ie. a column vector with a 1 in the last row and zeros elsewhere).

The argument `mpsorthog` can be used to set the gauge of the resulting MPS.

"""
function productstatemps(physdims::Dims, Dmax=1; state=:Vacuum, mpsorthog=:Right)

    N = length(physdims)

    if typeof(Dmax) <: Int
        Dmax=fill(Dmax, N-1)
    end
    
    bonddims = Vector{Int}(undef, N+1)
    bonddims[1] = 1
    bonddims[N+1] = 1

    for i=2:N
        bonddims[i] = min(Dmax[i-1], bonddims[i-1] * physdims[i-1])
    end
    for i=N:-1:2
        bonddims[i] = min(bonddims[i], bonddims[i+1] * physdims[i])
    end

    if state == :Vacuum
        statelist = [unitcol(1, physdims[i]) for i in 1:N]
    elseif state == :FullOccupy
        statelist = [unitcol(physdims[i], physdims[i]) for i in 1:N]
    elseif typeof(state)<:Vector
        statelist = state
        length(state)==N || throw(ErrorException("state list has length $(length(state)) while MPS has $N sites"))
    else
        throw(ErrorException("state input not recognised"))
    end
        
    B0 = Vector{Any}(undef, N)
    for i=1:N
        A = zeros(eltype(statelist[i]), bonddims[i], bonddims[i+1], physdims[i])
        for j=1:min(bonddims[i], bonddims[i+1])
            A[j,j,:] = statelist[i]
        end
        B0[i] = A
    end

    if mpsorthog==:Right
        mpsrightnorm!(B0)
    elseif mpsorthog==:Left
        mpsleftnorm!(B0)
    elseif typeof(mpsorthog)<:Int
        OC = mpsorthog
        mpsmixednorm!(B0, OC)
    end
    return B0
end

"""
    productstatemps(N::Int, d::Int, Dmax=1; state=:Vacuum, mpsorthog=:Right)

Return an `N`-site MPS with all local Hilbert space dimensions given by `d`. 

"""
productstatemps(N::Int, d::Int, Dmax=1; state=:Vacuum, mpsorthog=:Right) =
    productstatemps(ntuple(i->d, N), Dmax; state=state, mpsorthog=mpsorthog)

"""
    chainmps(N::Int, site::Int, numex::Int)

Generate an MPS with `numex` excitations on `site`

The returned MPS will have bond-dimensions and physical dimensions `numex+1`

"""
function chainmps(N::Int, site::Int, numex::Int)
    1 <= site <= N || throw(ErrorException("tried to excite site $site of $N-site chain"))
    D = numex + 1
    A = Vector{Any}(undef, N)
    for n=1:N
        M = zeros(D, D, D)
        M[:,:,1] = unitmat(D)
        for i=2:D
            j=i-1
            M[:,:,i] = diagm(j => fill(n==site ? 1/sqrt(factorial(j)) : 0.0, D-j))
        end
        A[n] = M
    end
    A[1] = A[1][1:1,:,:] * sqrt(factorial(numex))
    A[N] = A[N][:,D:D,:]
    mpsrightnorm!(A)
    return A
end

"""
    chainmps(N::Int, sites::Vector{Int}, numex::Int)

Generate an MPS with `numex` excitations of an equal super-position over `sites`

"""
function chainmps(N::Int, sites::Vector{Int}, numex::Int)
    all(x-> 1<=x<=N, sites) || throw(ErrorException("tried to excite non-existant site of $N-site chain"))
    D = numex + 1
    numsites = length(sites)
    A = Vector{Any}(undef, N)
    for n=1:N
        M = zeros(D, D, D)
        M[:,:,1] = unitmat(D)
        for i=2:D
            j=i-1
            M[:,:,i] = diagm(j => fill(in(n, sites) ? 1/sqrt(factorial(j)) : 0.0, D-j))
        end
        A[n] = M
    end
    A[1] = A[1][1:1,:,:] * sqrt(factorial(numex)) * (1/sqrt(numsites))^numex
    A[N] = A[N][:,D:D,:]
    mpsrightnorm!(A)
    return A
end

"""
    modemps(N::Int, k::Int, numex::Int, chainparams=[fill(1.0,N), fill(1.0,N-1)])

Generate an MPS with `numex` excitations of mode `k` of a bosonic tight-binding chain. 

`chainparams` takes the form `[e::Vector, t::Vector]` where `e` are the on-site energies and `t` are the hoppping
parameters.

The returned MPS will have bond-dimensions and physical dimensions `numex+1`

"""
function modemps(N::Int, k::Int, numex::Int, chainparams=[fill(1.0,N), fill(1.0,N-1)])
    1 <= k <= N || throw(ErrorException("tried to excite mode $k of an $N-site chain"))
    
    es=chainparams[1]
    ts=chainparams[2]
    hmat = diagm(0=>es, 1=>ts, -1=>ts)
    F = eigen(hmat)
    U = F.vectors
    
    D = numex + 1
    A = Vector{Any}(undef, N)
    for n=1:N
        M = zeros(eltype(U), D, D, D)
        M[:,:,1] = unitmat(D)
        for i=2:D
            j=i-1
            M[:,:,i] = diagm(j => fill(1/sqrt(factorial(j))*U[n,k]^j, D-j))
        end
        A[n] = M
    end
    A[1] = A[1][1:1,:,:] * sqrt(factorial(numex))
    A[N] = A[N][:,D:D,:]
    mpsrightnorm!(A)
    return A
end

"""
    modemps(N::Int, k::Vector{Int}, numex::Int, chainparams=[fill(1.0,N), fill(1.0,N-1)])

Generate an MPS with `numex` excitations of an equal superposition of modes `k` of a bosonic tight-binding chain.

"""
function modemps(N::Int, k::Vector{Int}, numex::Int, chainparams=[fill(1.0,N), fill(1.0,N-1)])
    all(x -> 1<=x<= N, k) || throw(ErrorException("tried to excite non-existant mode of an $N-site chain"))
    numks = length(k)
    
    es=chainparams[1]
    ts=chainparams[2]
    hmat = diagm(0=>es, 1=>ts, -1=>ts)
    F = eigen(hmat)
    U = F.vectors

    D = numex + 1
    A = Vector{Any}(undef, N)
    for n=1:N
        M = zeros(eltype(U), D, D, D)
        M[:,:,1] = unitmat(D)
        for i=2:D
            j=i-1
            M[:,:,i] = diagm(j => fill(1/sqrt(factorial(j))*sum(U[n,k])^j, D-j))
        end
        A[n] = M
    end
    A[1] = A[1][1:1,:,:] * sqrt(factorial(numex)) * (1/sqrt(numks))^numex
    A[N] = A[N][:,D:D,:]
    mpsrightnorm!(A)
    return A
end

"""
    electronkmps(N::Int, k::Int, spin=:Up, chainparams=[fill(1.0,N), fill(1.0,N-1)])

Generate an MPS for an electron with momentum `k`.

"""
function electronkmps(N::Int, k::Int, spin=:Up, chainparams=[fill(1.0,N), fill(1.0,N-1)])
    1 <= k <= N || throw(ErrorException("tried to excite mode $k of an $N-site chain"))
    
    es=chainparams[1]
    ts=chainparams[2]
    hmat = diagm(0=>es, 1=>ts, -1=>ts)
    F = eigen(hmat)
    U = F.vectors

    if spin==:Up
        i=2
    elseif spin==:Down
        i=3
    elseif spin==:UpDown
        i=4
    end
    
    A = Vector{Any}(undef, N)
    for n=1:N
        M = zeros(eltype(U), 2, 2, 4)
        M[:,:,1] = unitmat(2)
        M[:,:,i] = diagm(1 => [U[n,k]])
        A[n] = M
    end
    A[1] = A[1][1:1,:,:]
    A[N] = A[N][:,2:2,:]
    mpsrightnorm!(A)
    return A
end

"""
    electronkmps(N::Int, k::Vector{Int}, spin=:Up, chainparams=[fill(1.0,N), fill(1.0,N-1)])

Generate an MPS with 2 electrons in k-states `k1` and `k2`.

"""
function electron2kmps(N::Int, k1::Int, k2::Int, spin=:Up, chainparams=[fill(1.0,N), fill(1.0,N-1)])
    all(x -> 1<=x<= N, [k1,k2]) || throw(ErrorException("tried to excite non-existant mode of an $N-site chain"))

    es=chainparams[1]
    ts=chainparams[2]
    hmat = diagm(0=>es, 1=>ts, -1=>ts)
    F = eigen(hmat)
    U = F.vectors

    if spin==:Up
        i=2
    elseif spin==:Down
        i=3
    elseif spin==:UpDown
        i=4
    end
    
    A = Vector{Any}(undef, N)

    for n=1:N
        M = zeros(eltype(U), 4, 4, 4)
        M[:,:,1] = unitmat(4)
        M[:,:,i] = [0.0 U[n,k1] -U[n,k2] 0.0    ;
                    0.0 0.0      0.0     U[n,k2];
                    0.0 0.0      0.0     U[n,k1];
                    0.0 0.0      0.0     0.0    ]
        A[n] = M
    end
    A[1] = A[1][1:1,:,:]
    A[N] = A[N][:,4:4,:]
    mpsrightnorm!(A)
    return A
end

"""
    entanglemententropy(A)

For a list of tensors `A` representing a right orthonormalized MPS, compute the entanglement
entropy for a bipartite cut for every bond.

"""
function entanglemententropy(A)
    N = length(A)
    entropy = Vector{Float64}(undef, N-1)

    A1 = A[1]
    aleft, aright, aup = size(A1)
    U, S, V = svdtrunc(reshape(permutedims(A1, [1,3,2]), aleft*aup, :))
    schmidt = diag(S)
    entropy[1] = sum([schmidt[i]==0 ? 0 : -schmidt[i]^2 * log(schmidt[i]^2) for i=1:length(schmidt)])
    for k = 2:N-1
        Ak = A[k]
        aleft, aright, aup = size(Ak)
        Ak = reshape(S*V*reshape(Ak, aleft, :), aleft, aright, aup)
        U, S, V = svdtrunc(reshape(permutedims(Ak, [1,3,2]), aleft*aup, :))
        schmidt = diag(S)
        entropy[k] = sum([schmidt[i]==0 ? 0 : -schmidt[i]^2 * log(schmidt[i]^2) for i=1:length(schmidt)])
    end
    return entropy
end
#If you want the entanglements at every point during an evolution it would be more efficient to
#replace the QRs with SVDs in the tdvp routine

"""
    svdmps(A)

For a right normalised mps `A` compute the full svd spectrum for a bipartition at every bond.
"""
function svdmps(A)
    N = length(A)
    spec = Vector{Any}(undef, N-1)    
    A1 = A[1]
    aleft, aright, aup = size(A1)
    U, S, V = svdtrunc(reshape(permutedims(A1, [1,3,2]), aleft*aup, :))
    spec[1] = diag(S)
    for k = 2:N-1
        Ak = A[k]
        aleft, aright, aup = size(Ak)
        Ak = reshape(S*V*reshape(Ak, aleft, :), aleft, aright, aup)
        U, S, V = svdtrunc(reshape(permutedims(Ak, [1,3,2]), aleft*aup, :))
        spec[k] = diag(S)
    end
    return spec
end

"""
    elementmps(A, el...)

Return the element of the MPS `A` for the set of physical states `el...`

# Examples

```julia-repl
julia> A = chainmps(6, [2,4], 1);

julia> elementmps(A, 1, 2, 1, 1, 1, 1)
0.7071067811865475

julia> elementmps(A, 1, 1, 1, 2, 1, 1)
0.7071067811865475

julia> elementmps(A, 1, 2, 1, 2, 1, 1)
0.0

julia> elementmps(A, 1, 1, 1, 1, 1, 1)
0.0
```
"""
function elementmps(A, el...)
    nsites = length(A)
    length(el) == nsites || throw(ArgumentError("indices do not match MPS"))
    c = A[1][:,:,el[1]]
    for i=2:nsites
        c = c * A[i][:,:,el[i]]
    end
    c[1,1]
end

"""
    elementmpo(M, el...)

Return the element of the MPO `M` for the set of physical states `el...`
"""
function elementmpo(M, el...)
    nsites = length(M)
    length(el) == nsites || throw(ArgumentError("indices do not match MPO"))
    c = M[1][:,:,el[1][1],el[1][2]]
    for i=2:nsites
        c *= M[i][:,:,el[i][1],el[i][2]]
    end
    c[1,1]
end

function apply1siteoperator!(A, O, sites::Vector{Int})
    for i in sites
        @tensor R[a,b,s] := O[s,s']*A[i][a,b,s']
        A[i] = R
    end
end
apply1siteoperator!(A, O, site::Int) = apply1siteoperator!(A, O, [site])

function applympo!(A, H)
    N = length(H)
    N == length(A) || throw(ArgumentError("MPO has $N site while MPS has $(length(A)) sites"))
    for i=1:N
        Al, Ar, d = size(A[i])
        Hl, Hr, d, d = size(H[i])
        @tensor X[a',a,b',b,s] := H[i][a',b',s,s'] * A[i][a,b,s']
        A[i] = reshape(X, Al*Hl, Ar*Hr, d)
    end
end

function reversemps!(A)
    N = length(A)
    reverse!(A)
    for i=1:N
        A[i] = permutedims(A[i], [2,1,3])
    end
end
function reversemps(A)
    N = length(A)
    Ar = Vector{Any}(undef, N)
    for i=1:N
        Ar[N-i+1] = permutedims(A[i], [2,1,3])
    end
    return Ar
end

function reversempo!(M)
    N = length(M)
    reverse!(M)
    for i=1:N
        M[i] = permutedims(M[i], [2,1,3,4])
    end
end
function reversempo(M)
    N = length(M)
    Mr = Vector{Any}(undef, N)
    for i=1:N
        Mr[N-i+1] = permutedims(M[i], [2,1,3,4])
    end
    return Mr
end

"""
    physdims(M)

Return the physical dimensions of an MPS or MPO `M`.
"""
function physdims(M::Vector)
    N = length(M)
    res = Vector{Int}(undef, N)
    for (i, site) in enumerate(M)
        res[i] = size(site)[end]
    end
    return Dims(res)
end

function rtrunc(A, D::Int)
    A[1:D,:,:]
end
function ltrunc(A, D::Int)
    A[:,1:D,:]
end
#storing two copys of mps is not so memory efficient
function mpsform(A_::Vector)
    A = deepcopy(A_)
    N = length(A)
    for i=1:N-1
        D = min(size(A[i])[2], size(A[i+1])[1])
        A[i] = ltrunc(A[i], D)
        A[i+1] = rtrunc(A[i+1], D)
    end
    return A
end
function mpsform!(A::Vector)
    N = length(A)
    for i=1:N-1
        D = min(size(A[i])[2], size(A[i+1])[1])
        A[i] = ltrunc(A[i], D)
        A[i+1] = rtrunc(A[i+1], D)
    end
end

function bonddimsmps(A::Vector)
    N = length(A)
    res = Vector{Int}(undef, N+1)
    res[1] = res[N+1] = 1
    for i=1:N-1
        res[i+1] = min(size(A[i])[2], size(A[i+1])[1])
    end
    return Dims(res)
end
function bonddims(M::Vector)
    N = length(M)
    res = Vector{Int}(undef, N+1)
    res[1] = res[N+1] = 1
    for i=1:N-1
        res[i+1] = size(M[i])[2]
    end
    return Dims(res)
end

#calculates M1*M2 where M1 and M2 are MPOs
function multiply(M1::Vector, M2::Vector)
    N = length(M1)
    length(M2) == N || throw(ArgumentError("MPOs do not have the same length!"))
    res = Vector{Any}(undef, N)
    for k=1:N
        Dl1, Dr1, d, d = size(M1[k])
        Dl2, Dr2, d, d = size(M2[k])
        @tensor M12[a1,a2,b1,b2,s1,s2] := M1[k][a1,b1,s1,s] * M2[k][a2,b2,s,s2]
        res[k] = reshape(M12, Dl1*Dl2, Dr1*Dr2, d, d)
    end
    return res
end

"""
    mpsembed(A::Vector, Dmax::Int)

Embed MPS `A` in manifold of max bond-dimension `Dmax`

"""
function mpsembed!(A::Vector, Dmax::Int)

    N = length(A)
    pdims = physdims(A)
    
    bonddims = Vector{Int}(undef, N+1)
    bonddims[1] = 1
    bonddims[N+1] = 1

    for i=2:N
        bonddims[i] = min(Dmax, bonddims[i-1] * pdims[i-1])
    end
    for i=N:-1:2
        bonddims[i] = min(bonddims[i], bonddims[i+1] * pdims[i])
    end

    for i=1:N
        A[i] = setleftbond(A[i], bonddims[i])
        A[i] = setrightbond(A[i], bonddims[i+1])
    end
    mpsrightnorm!(A)
    return A
end

