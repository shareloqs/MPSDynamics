abstract type Observable end

struct OneSiteObservable <: Observable
    name::String
    op::Array{<:Number, 2}
    sites::Union{Int, Tuple{Int,Int}, Vector{Int}, Nothing}
    hermitian::Bool
    allsites::Bool
end

"""
    OneSiteObservable(name,op,sites)

Computes the local expectation value of the one-site operator `op` on the specified sites. Used to define
one-site observables that are obs and convobs parameters for the `runsim` function.

"""
OneSiteObservable(name, op, sites) = OneSiteObservable(name, op, sites, ishermitian(op), false)

"""
    OneSiteObservable(name,op)

Computes the local expectation value of the one-site operator `op` on the every site. Used to define
one-site observables that are obs and convobs parameters for the `runsim` function.

"""

OneSiteObservable(name, op) = OneSiteObservable(name, op, nothing, ishermitian(op), true)

struct TwoSiteObservable <: Observable
    name::String
    op1::Array{<:Number, 2}
    op2::Array{<:Number, 2}
    sites1::Union{Int, Tuple{Int, Int}, Vector{Int}, Nothing}
    sites2::Union{Int, Tuple{Int, Int}, Vector{Int}, Nothing}
    allsites::Bool
end

struct CdagCup <: Observable
    name::String
    sites::Tuple{Int,Int}
end
CdagCup(sites::Tuple{Int,Int}) = CdagCup("CdagCup", sites)
CdagCup(i1::Int, i2::Int) = CdagCup("CdagCup", (i1,i2))

struct CdagCdn <: Observable
    name::String
    sites::Tuple{Int,Int}
end
CdagCdn(sites::Tuple{Int,Int}) = CdagCdn("CdagCup", sites)
CdagCdn(i1::Int, i2::Int) = CdagCdn("CdagCdn", (i1,i2))

"""
    TwoSiteObservable(name,op1,op2,sites1=nothing,sites2=nothing)

Computes the local expectation value of operators `op1` and `op2` where `op1` acts on sites1 and `op2` acts on sites2. Used to define
several-site observables that are obs and convobs parameters for the `runsim` function.

"""

function TwoSiteObservable(name, op1, op2, sites1=nothing, sites2=nothing)
    return TwoSiteObservable(name, op1, op2, sites1, sites2, sites1==nothing && sites2==nothing)
end

function Base.ndims(ob::OneSiteObservable)
    if typeof(ob.sites) <: Int
        return 0
    else
        return 1
    end
end
function Base.ndims(ob::TwoSiteObservable)
    if typeof(ob.sites1) <: Int && typeof(ob.sites2) <: Int
        return 0
    elseif typeof(ob.sites1) <: Int || typeof(ob.sites2) <: Int
        return 1
    else
        return 2
    end
end
Base.ndims(::CdagCup) = 2
Base.ndims(::CdagCdn) = 2

span(A, ob::OneSiteObservable) = span(A, ob.sites)
span(A, sites::Tuple) = collect(sites[1]:(sites[1] <= sites[2] ? 1 : -1):sites[2])
span(A, site::Int) = Int[site]
span(A, sites::Vector) = sites
span(A, sites::Nothing) = collect(1:length(A))

reach(A, ob::OneSiteObservable) = max(span(A, ob)...)
reach(A, ob::TwoSiteObservable) = min(max(span(A, ob.sites1)...), max(span(A, ob.sites2)...))
reach(A, ob::CdagCup) = max(ob.sites...)
reach(A, ob::CdagCdn) = max(ob.sites...)
reach(A, ob::Observable) = 1

"""
    measurempo(A::Vector, M::Vector)

For a list of tensors `A` representing a right orthonormalized MPS, compute the local expectation
value of the MPO M on every site.

"""

function measurempo(A::Vector, M::Vector)
    N = length(M)
    N == length(A) || throw(ArgumentError("MPO has $N site while MPS has $(length(A)) sites"))
    F = fill!(similar(M[1], (1,1,1)), 1)
    for k=1:N
        F = updateleftenv(A[k], M[k], F)
    end
    real(only(F))
end
"""
    measurempo(A::Vector, M::Vector, sites::Tuples{Int,Int})

For a list of tensors `A` representing a right orthonormalized MPS, compute the local expectation
value of the MPO M on specified sites.

"""

function measurempo(A::Vector, M::Vector, sites::Tuple{Int, Int})
    N = sites[2] - sites[1] + 1
    F = fill!(similar(M[1], (1,1,1)), 1)
    for k=1:sites[1]-1
        d = size(A[k])[3]
        id = reshape(unitmat(d), 1, 1, d, d)
        F = updateleftenv(A[k], id, F)
    end
    for k=sites[1]:sites[2]
        F = updateleftenv(A[k], M[k], F)
    end
    F = tensortrace(F, [1,2,1], [2])
    real(only(F))
end

"""
    measure(A, O; kwargs...)

Measure observable `O` on mps state `A`

"""
measure(A, O::OneSiteObservable; acs=nothing, ρ=nothing, kwargs...) = measure(A, O, acs, ρ)
measure(A, O::OneSiteObservable, ::Nothing, ::Nothing) = measure1siteoperator(A, O.op, O.sites)
#ignore ρ if acs is supplied:
measure(A, O::OneSiteObservable, acs::Vector, ρ::Vector) = measure(A, O, acs, nothing)
function measure(A, O::OneSiteObservable, acs::Vector, ::Nothing)
    T = O.hermitian ? Float64 : ComplexF64
    if typeof(O.sites) <: Int
        k = O.sites
        v = ACOAC(acs[k], O.op)
        T<:Real && (v=real(v))
        return v
    else
        sites = span(A, O)
    end
    N = length(sites)
    expval = zeros(T, N)
    for (i, k) in enumerate(sites)
        v = ACOAC(acs[k], O.op)
        T<:Real && (v=real(v))
        expval[i] = v
    end
    return expval
end
function measure(A, O::OneSiteObservable, ::Nothing, ρ::Vector)
    T = O.hermitian ? Float64 : ComplexF64
    if typeof(O.sites) <: Int
        k = O.sites
        v = rhoAOAstar(ρ[k], A[k], O.op, nothing)
        T<:Real && (v=real(v))
        return v
    else
        sites = span(A, O)
    end
    N = length(sites)
    expval = zeros(T, N)
    for (i, k) in enumerate(sites)
        v = rhoAOAstar(ρ[k], A[k], O.op, nothing)
        T<:Real && (v=real(v))
        expval[i] = v
    end
    return expval
end

measure(A, O::TwoSiteObservable; ρ=nothing, kwargs...) = measure(A, O, ρ)
measure(A, O::TwoSiteObservable, ::Nothing) =
    measure2siteoperator(A, O.op1, O.op2, O.sites1, O.sites2)
measure(A, O::TwoSiteObservable, ρ::Vector) =
    measure2siteoperator(A, O.op1, O.op2, O.sites1, O.sites2, ρ)

"""
    measure1siteoperator(A::Vector, O, sites::Vector{Int})

For a list of tensors `A` representing a right orthonormalized MPS, compute the local expectation
value of a one-site operator O for every site or just one if it is specified.

For calculating operators on single sites this will be more efficient if the site is on the left of the mps.

"""
function measure1siteoperator(A::Vector, O, sites::Vector{Int})
    N = length(A)
    ρ = ones(ComplexF64, 1, 1)
    T = ishermitian(O) ? Float64 : ComplexF64
    expval = zeros(T, N)
    for i=1:N
        if in(i, sites)
            v = rhoAOAstar(ρ, A[i], O, nothing)
            T<:Real && (v=real(v))
            expval[i] = v
        end
        ρ = rhoAAstar(ρ, A[i])
    end
    return expval[sites]
end

"""
    measure1siteoperator(A::Vector, O, chainsection::Tuple{Int64,Int64})
For a list of tensors `A` representing a right orthonormalized MPS, compute the local expectation
value of a one-site operator O for a chainsection.

"""

function measure1siteoperator(A::Vector, O, chainsection::Tuple{Int64,Int64})
    ρ = ones(ComplexF64, 1, 1)

    T = ishermitian(O) ? Float64 : ComplexF64

    l=min(chainsection...)#leftmost site in section
    r=max(chainsection...)#rightmost site in section
    rev = l != chainsection[1]
    N = r-l+1
    expval = zeros(T, N)

    for i=1:l-1
        ρ = rhoAAstar(ρ, A[i])
    end
    for i=l:r-1
        v = rhoAOAstar(ρ, A[i], O, nothing)
        ρ = rhoAAstar(ρ, A[i])
        T<:Real && (v=real(v))
        expval[i-l+1] = v
    end
    v = rhoAOAstar(ρ, A[r], O, nothing)
    T<:Real && (v=real(v))
    expval[N] = v
    if rev
        reverse!(expval)
    end
    return expval
end

"""
    measure1siteoperator(A::Vector, O)
For a list of tensors `A` representing a right orthonormalized MPS, compute the local expectation
value of a one-site operator O for every site.

"""
function measure1siteoperator(A::Vector, O)
    N = length(A)
    ρ = ones(ComplexF64, 1, 1)
    T = ishermitian(O) ? Float64 : ComplexF64
    expval = zeros(T, N)
    for i=1:N
        v = rhoAOAstar(ρ, A[i], O, nothing)
        T<:Real && (v=real(v))
        expval[i] = v
        ρ = rhoAAstar(ρ, A[i])
    end
    return expval
end
measure1siteoperator(A::Vector, O, ::Nothing) = measure1siteoperator(A, O)

"""
    measure1siteoperator(A::Vector, O, site::Int)
For a list of tensors `A` representing a right orthonormalized MPS, compute the local expectation
value of a one-site operator O for a single site.

"""

function measure1siteoperator(A::Vector, O, site::Int)
    ρ = ones(ComplexF64, 1, 1)
    T = ishermitian(O) ? Float64 : ComplexF64
    for i=1:site-1
        ρ = rhoAAstar(ρ, A[i])
    end
    v = rhoAOAstar(ρ, A[site], O, nothing)
    T<:Real && (v=real(v))
    return v
end

"""
     measure2siteoperator(A::Vector, M1, M2, j1, j2)

Caculate expectation of M1*M2 where M1 acts on site j1 and M2 acts on site j2, assumes A is right normalised.

"""
function measure2siteoperator(A::Vector, M1, M2, j1::Int64, j2::Int64)
    herm_trans = ishermitian(M1) && ishermitian(M2)
    ρ = ones(ComplexF64, 1, 1)
    i1=min(j1,j2)
    i2=max(j1,j2)
    m1 = j1<j2 ? M1 : M2
    m2 = j1<j2 ? M2 : M1
    if j1==j2
        return measure1siteoperator(A, M1*M2, j1)
    else
        T = herm_trans ? Float64 : ComplexF64
        for k=1:i1-1
            ρ = rhoAAstar(ρ, A[k])
        end
        ρ = rhoAOAstar(ρ, A[i1], m1)
        for k=i1+1:i2-1
            ρ = rhoAAstar(ρ, A[k])
        end
        v = rhoAOAstar(ρ, A[i2], m2, nothing)
        T<:Real && (v=real(v))
        return v
    end
end
function measure2siteoperator(A::Vector, M1, M2, j1::Int64, j2::Int64, ρ::Vector)
    herm_trans = ishermitian(M1) && ishermitian(M2)
    T = herm_trans ? Float64 : ComplexF64
    i1=min(j1,j2)
    i2=max(j1,j2)
    m1 = j1<j2 ? M1 : M2
    m2 = j1<j2 ? M2 : M1

    if j1==j2
        herm_cis = ishermitian(M1*M2)
        v = rhoAOAstar(ρ[j1], A[j1], M1*M2, nothing)
        return herm_cis ? real(v) : v
    else
        ρ1 = rhoAOAstar(ρ[i1], A[i1], m1)
        for k=i1+1:i2-1
            ρ1 = rhoAAstar(ρ1, A[k])
        end
        v = rhoAOAstar(ρ1, A[i2], m2, nothing)
        T<:Real && (v=real(v))
        return v
    end
end

function measure2siteoperator(A::Vector, M1, M2)
    herm_cis = ishermitian(M1*M2)
    herm_trans = ishermitian(M1) && ishermitian(M2)
    pair = M1 == M2
    cpair = M1 == M2'
    if pair || cpair
        return measure2siteoperator_pair(A, M1, conjugate=!pair)
    else
        N = length(A)
        ρ = ones(ComplexF64, 1, 1)
        
        T = (herm_cis && herm_trans) ? Float64 : ComplexF64
        
        expval = zeros(T, N, N)

        for i in 1:N
            v = rhoAOAstar(ρ, A[i], M1*M2, nothing)
            herm_cis && (v=real(v))
            expval[i,i] = v
            ρ12 = rhoAOAstar(ρ, A[i], M1)
            ρ21 = rhoAOAstar(ρ, A[i], M2)
            for j in i+1:N
                v = rhoAOAstar(ρ12, A[j], M2, nothing)
                herm_trans && (v=real(v))
                expval[i,j] = v
                v = rhoAOAstar(ρ21, A[j], M1, nothing)
                herm_trans && (v=real(v))
                expval[j,i] = v
                ρ12 = rhoAAstar(ρ12, A[j])
                ρ21 = rhoAAstar(ρ21, A[j])
            end
            ρ = rhoAAstar(ρ, A[i])
        end
        return expval
    end
end
function measure2siteoperator(A::Vector, M1, M2, ρ::Vector)
    herm_cis = ishermitian(M1*M2)
    herm_trans = ishermitian(M1) && ishermitian(M2)
    pair = M1 == M2
    cpair = M1 == M2'
    if pair || cpair
        return measure2siteoperator_pair(A, M1, ρ, conjugate=!pair)
    else
        N = length(A)
        
        T = (herm_cis && herm_trans) ? Float64 : ComplexF64
        
        expval = zeros(T, N, N)

        for i in 1:N
            v = rhoAOAstar(ρ[i], A[i], M1*M2, nothing)
            herm_cis && (v=real(v))
            expval[i,i] = v
            ρ12 = rhoAOAstar(ρ[i], A[i], M1)
            ρ21 = rhoAOAstar(ρ[i], A[i], M2)
            for j in i+1:N
                v = rhoAOAstar(ρ12, A[j], M2, nothing)
                herm_trans && (v=real(v))
                expval[i,j] = v
                v = rhoAOAstar(ρ21, A[j], M1, nothing)
                herm_trans && (v=real(v))
                expval[j,i] = v
                ρ12 = rhoAAstar(ρ12, A[j])
                ρ21 = rhoAAstar(ρ21, A[j])
            end
        end
        return expval
    end
end
measure2siteoperator(A::Vector, M1, M2, ::Nothing, ::Nothing) = measure2siteoperator(A, M1, M2)
measure2siteoperator(A::Vector, M1, M2, ::Nothing, ::Nothing, ρ::Vector) =
    measure2siteoperator(A, M1, M2, ρ)

function measure2siteoperator_pair(A::Vector, M1; conjugate=false)
    M2 = conjugate ? Matrix(M1') : M1    
    herm_cis = ishermitian(M1*M2)
    herm_trans = ishermitian(M1) && ishermitian(M2)

    N = length(A)
    ρ = ones(ComplexF64, 1, 1)

    T = (herm_cis && herm_trans) ? Float64 : ComplexF64
    
    expval = zeros(T, N, N)
    for i in 1:N
        v = rhoAOAstar(ρ, A[i], M1*M2, nothing)
        herm_cis && (v=real(v))
        expval[i,i] = v
        ρ12 = rhoAOAstar(ρ, A[i], M1)
        for j in i+1:N
            v = rhoAOAstar(ρ12, A[j], M2, nothing)
            herm_trans && (v=real(v))
            expval[i,j] = v
            ρ12 = rhoAAstar(ρ12, A[j])
        end
        ρ = rhoAAstar(ρ, A[i])
    end
    dia = diagm(0 => diag(expval))
    return expval + (conjugate ? expval' : transpose(expval)) - dia
end
function measure2siteoperator_pair(A::Vector, M1, ρ::Vector; conjugate=false)
    M2 = conjugate ? Matrix(M1') : M1
    herm_cis = ishermitian(M1*M2)
    herm_trans = ishermitian(M1) && ishermitian(M2)

    N = length(A)

    T = (herm_cis && herm_trans) ? Float64 : ComplexF64
    
    expval = zeros(T, N, N)
    for i in 1:N
        v = rhoAOAstar(ρ[i], A[i], M1*M2, nothing)
        herm_cis && (v=real(v))
        expval[i,i] = v
        ρ12 = rhoAOAstar(ρ[i], A[i], M1)
        for j in i+1:N
            v = rhoAOAstar(ρ12, A[j], M2, nothing)
            herm_trans && (v=real(v))
            expval[i,j] = v
            ρ12 = rhoAAstar(ρ12, A[j])
        end
    end
    dia = diagm(0 => diag(expval))
    return expval + (conjugate ? expval' : transpose(expval)) - dia
end

"""
     measure2siteoperator(A::Vector, M1, M2, sites1::Vector{Int}, sites2::Vector{Int})

Caculate expectation of M1*M2 where M1 acts on sites1 and M2 acts on sites2, assumes A is right normalised.

"""

function measure2siteoperator(A::Vector, M1, M2, sites1::Vector{Int}, sites2::Vector{Int})
    if size(M1) == size(M2)
        herm_cis = ishermitian(M1*M2)
    else
        herm_cis = false
    end
    herm_trans = ishermitian(M1) && ishermitian(M2)

    maxsites1 = max(sites1...)
    maxsites2 = max(sites2...)
    maxsites = max(maxsites1, maxsites2)
    
    N = length(A)
    ρ = ones(ComplexF64, 1, 1)

    T = herm_trans ? Float64 : ComplexF64

    expval = zeros(T, N, N)
    for i in 1:maxsites
        if in(i, sites1)
            if in(i, sites2)
                v = rhoAOAstar(ρ, A[i], M1*M2, nothing)
                expval[i,i] = herm_cis ? real(v) : v
            end
            if any(x->x>i, sites2)
                ρ12 = rhoAOAstar(ρ, A[i], M1)
                for j in i+1:maxsites2
                    if in(j, sites2)
                        v = rhoAOAstar(ρ12, A[j], M2, nothing)
                        expval[i,j] = herm_trans ? real(v) : v
                    end
                    ρ12 = rhoAAstar(ρ12, A[j])
                end
            end
        end
        
        if in(i, sites2)
            if any(x->x>i,sites1)
                ρ21 = rhoAOAstar(ρ, A[i], M2)
                for j in i+1:maxsites1
                    if in(j, sites1)
                        v = rhoAOAstar(ρ21, A[j], M1, nothing)
                        expval[j,i] = herm_trans ? real(v) : v
                    end
                    ρ21 = rhoAAstar(ρ21, A[j])
                end
            end
        end
        ρ = rhoAAstar(ρ, A[i])
    end
    return expval[sites1,sites2]
end
function measure2siteoperator(A::Vector, M1, M2, sites1::Vector{Int}, sites2::Vector{Int}, ρ::Vector)
    if size(M1) == size(M2)
        herm_cis = ishermitian(M1*M2)
    else
        herm_cis = false
    end
    herm_trans = ishermitian(M1) && ishermitian(M2)

    maxsites1 = max(sites1...)
    maxsites2 = max(sites2...)
    maxsites = max(maxsites1, maxsites2)

    N = length(A)

    T = (herm_cis && herm_trans) ? Float64 : ComplexF64

    expval = zeros(T, N, N)
    for i in 1:maxsites
        if in(i, sites1)
            if in(i, sites2)
                v = rhoAOAstar(ρ[i], A[i], M1*M2, nothing)
                expval[i,i] = herm_cis ? real(v) : v
            end
            if any(x->x>i, sites2)
                ρ12 = rhoAOAstar(ρ[i], A[i], M1)
                for j in i+1:maxsites2
                    if in(j, sites2)
                        v = rhoAOAstar(ρ12, A[j], M2, nothing)
                        expval[i,j] = herm_trans ? real(v) : v
                    end
                    ρ12 = rhoAAstar(ρ12, A[j])
                end
            end
        end
        
        if in(i, sites2)
            if any(x->x>i,sites1)
                ρ21 = rhoAOAstar(ρ[i], A[i], M2)
                for j in i+1:maxsites1
                    if in(j, sites1)
                        v = rhoAOAstar(ρ21, A[j], M1, nothing)
                        expval[j,i] = herm_trans ? real(v) : v
                    end
                    ρ21 = rhoAAstar(ρ21, A[j])
                end
            end
        end
    end
    return expval[sites1,sites2]
end

function measure2siteoperator(A::Vector, M1, M2, sites::Vector{Int})
    pair = M1 == M2
    cpair = M1 == M2'
    if pair || cpair
        return measure2siteoperator_pair(A, M1, sites, conjugate=!pair)
    else
        return measure2siteoperator(A, M1, M2, sites, sites)
    end
end
function measure2siteoperator(A::Vector, M1, M2, sites::Vector{Int}, ρ::Vector)
    pair = M1 == M2
    cpair = M1 == M2'
    if pair || cpair
        return measure2siteoperator_pair(A, M1, sites, ρ, conjugate=!pair)
    else
        return measure2siteoperator(A, M1, M2, sites, sites, ρ)
    end
end
measure2siteoperator(A::Vector, M1, M2, sites::Vector{Int}, ::Nothing) =
    measure2siteoperator(A, M1, M2, sites)
measure2siteoperator(A::Vector, M1, M2, sites::Vector{Int}, ::Nothing, ρ::Vector) =
    measure2siteoperator(A, M1, M2, sites, ρ)

function measure2siteoperator_pair(A::Vector, M1, sites::Vector{Int}; conjugate=false)
    M2 = conjugate ? Matrix(M1') : M1
    herm_cis = ishermitian(M1*M2)
    herm_trans = ishermitian(M1) && ishermitian(M2)
    
    N = length(A)
    ρ = ones(ComplexF64, 1, 1)

    T = (herm_cis && herm_trans) ? Float64 : ComplexF64

    expval = zeros(T, N, N)
    for i in 1:N
        if in(i, sites)
            v = rhoAOAstar(ρ, A[i], M1*M2, nothing)
            herm_cis && (v=real(v))
            expval[i,i] = v

            ρ12 = rhoAOAstar(ρ, A[i], M1)
            for j in i+1:N
                if in(j, sites)
                    v = rhoAOAstar(ρ12, A[j], M2, nothing)
                    herm_trans && (v=real(v))
                    expval[i,j] = v
                end
                ρ12 = rhoAAstar(ρ12, A[j])
            end
        end
        ρ = rhoAAstar(ρ, A[i])
    end
    expval = expval[sites, sites]
    dia = diagm(0 => diag(expval))
    return expval + (conjugate ? expval' : transpose(expval)) - dia
end
function measure2siteoperator_pair(A::Vector, M1, sites::Vector{Int}, ρ::Vector; conjugate=false)
    M2 = conjugate ? Matrix(M1') : M1
    herm_cis = ishermitian(M1*M2)
    herm_trans = ishermitian(M1) && ishermitian(M2)

    N = length(A)

    T = (herm_cis && herm_trans) ? Float64 : ComplexF64

    expval = zeros(T, N, N)
    for i in 1:N
        if in(i, sites)
            v = rhoAOAstar(ρ[i], A[i], M1*M2, nothing)
            herm_cis && (v=real(v))
            expval[i,i] = v

            ρ12 = rhoAOAstar(ρ[i], A[i], M1)
            for j in i+1:N
                if in(j, sites)
                    v = rhoAOAstar(ρ12, A[j], M2, nothing)
                    herm_trans && (v=real(v))
                    expval[i,j] = v
                end
                ρ12 = rhoAAstar(ρ12, A[j])
            end
        end
    end
    expval = expval[sites, sites]
    dia = diagm(0 => diag(expval))
    return expval + (conjugate ? expval' : transpose(expval)) - dia
end

function measure2siteoperator(A::Vector, M1, M2, chainsection::Tuple{Int,Int})
    herm_cis = ishermitian(M1*M2)
    herm_trans = ishermitian(M1) && ishermitian(M2)
    pair = M1 == M2
    cpair = M1 == M2'
    if pair || cpair
        return measure2siteoperator_pair(A, M1, chainsection; conjugate=!pair)
    end
    N = length(A)
    ρ = ones(ComplexF64, 1, 1)
    
    l=min(chainsection...)#leftmost site in section
    r=max(chainsection...)#rightmost site in section
    rev = l != chainsection[1]
    N = r-l+1

    T = (herm_cis && herm_trans) ? Float64 : ComplexF64

    expval = zeros(T, N, N)

    for i in 1:l-1
        ρ = rhoAAstar(ρ, A[i])
    end
    for i in l:r
        v = rhoAOAstar(ρ, A[i], M1*M2, nothing)
        herm_cis && (v=real(v))
        expval[i-l+1,i-l+1] = v
        ρ12 = rhoAOAstar(ρ, A[i], M1)
        ρ21 = rhoAOAstar(ρ, A[i], M2)
        for j in i+1:r
            v = rhoAOAstar(ρ12, A[j], M2, nothing)
            herm_trans && (v=real(v))
            expval[i-l+1,j-l+1] = v
            v = rhoAOAstar(ρ21, A[j], M1, nothing)
            herm_trans && (v=real(v))
            expval[j-l+1,i-l+1] = v
            ρ12 = rhoAAstar(ρ12, A[j])
            ρ21 = rhoAAstar(ρ21, A[j])
        end
        ρ = rhoAAstar(ρ, A[i])
    end
    if rev
        expval = reverse(reverse(expval, dims=1), dims=2)
    end
    return expval
end
function measure2siteoperator(A::Vector, M1, M2, chainsection::Tuple{Int,Int}, ρ::Vector)
    herm_cis = ishermitian(M1*M2)
    herm_trans = ishermitian(M1) && ishermitian(M2)
    pair = M1 == M2
    cpair = M1 == M2'
    if pair || cpair
        return measure2siteoperator_pair(A, M1, chainsection, ρ, conjugate=!pair)
    end
    N = length(A)
    
    l=min(chainsection...)#leftmost site in section
    r=max(chainsection...)#rightmost site in section
    rev = l != chainsection[1]
    N = r-l+1

    T = (herm_cis && herm_trans) ? Float64 : ComplexF64

    expval = zeros(T, N, N)

    for i in l:r
        v = rhoAOAstar(ρ[i], A[i], M1*M2, nothing)
        herm_cis && (v=real(v))
                expval[i-l+1,i-l+1] = v
        ρ12 = rhoAOAstar(ρ[i], A[i], M1)
        ρ21 = rhoAOAstar(ρ[i], A[i], M2)
        for j in i+1:r
            v = rhoAOAstar(ρ12, A[j], M2, nothing)
            herm_trans && (v=real(v))
            expval[i-l+1,j-l+1] = v
            v = rhoAOAstar(ρ21, A[j], M1, nothing)
            herm_trans && (v=real(v))
            expval[j-l+1,i-l+1] = v
            ρ12 = rhoAAstar(ρ12, A[j])
            ρ21 = rhoAAstar(ρ21, A[j])
        end
    end
    if rev
        expval = reverse(reverse(expval, dims=1), dims=2)
    end
    return expval
end
measure2siteoperator(A, M1, M2, chainsection::Tuple{Int,Int}, ::Nothing) =
    measure2siteoperator(A, M1, M2, chainsection)

function measure2siteoperator_pair(A::Vector, M1, chainsection::Tuple{Int,Int}; conjugate=false)
    M2 = conjugate ? Matrix(M1') : M1
    herm_cis = ishermitian(M1*M2)
    herm_trans = ishermitian(M1) && ishermitian(M2)

    N = length(A)
    ρ = ones(ComplexF64, 1, 1)
    
    l=min(chainsection...)#leftmost site in section
    r=max(chainsection...)#rightmost site in section
    rev = l != chainsection[1]
    N = r-l+1

    T = (herm_cis && herm_trans) ? Float64 : ComplexF64

    expval = zeros(T, N, N)

    for i in 1:l-1
        ρ = rhoAAstar(ρ, A[i])
    end
    for i in l:r
        v = rhoAOAstar(ρ, A[i], M1*M2, nothing)
        herm_cis && (v=real(v))
        expval[i-l+1,i-l+1] = v
        ρ12 = rhoAOAstar(ρ, A[i], M1)
        for j in i+1:r
            v = rhoAOAstar(ρ12, A[j], M2, nothing)
            herm_trans && (v=real(v))
            expval[i-l+1,j-l+1] = v
            ρ12 = rhoAAstar(ρ12, A[j])
        end
        ρ = rhoAAstar(ρ, A[i])
    end
    dia = diagm(0 => diag(expval))
    expval = expval + (conjugate ? expval' : transpose(expval)) - dia
    if rev
        expval = reverse(reverse(expval, dims=1), dims=2)
    end
    return expval
end
function measure2siteoperator_pair(A::Vector, M1, chainsection::Tuple{Int,Int}, ρ::Vector; conjugate=false)
    M2 = conjugate ? Matrix(M1') : M1
    herm_cis = ishermitian(M1*M2)
    herm_trans = ishermitian(M1) && ishermitian(M2)

    N = length(A)
    
    l=min(chainsection...)#leftmost site in section
    r=max(chainsection...)#rightmost site in section
    rev = l != chainsection[1]
    N = r-l+1

    T = (herm_cis && herm_trans) ? Float64 : ComplexF64

    expval = zeros(T, N, N)

    for i in l:r
        v = rhoAOAstar(ρ[i], A[i], M1*M2, nothing)
        herm_cis && (v=real(v))
        expval[i-l+1,i-l+1] = v
        ρ12 = rhoAOAstar(ρ[i], A[i], M1)
        for j in i+1:r
            v = rhoAOAstar(ρ12, A[j], M2, nothing)
            herm_trans && (v=real(v))
            expval[i-l+1,j-l+1] = v
            ρ12 = rhoAAstar(ρ12, A[j])
        end
    end
    dia = diagm(0 => diag(expval))
    expval = expval + (conjugate ? expval' : transpose(expval)) - dia
    if rev
        expval = reverse(reverse(expval, dims=1), dims=2)
    end
    return expval
end
measure2siteoperator(A, M1, M2, chainsection::Tuple{Int,Int}, ::Nothing, ::Nothing) =
    measure2siteoperator(A, M1, M2, chainsection)
measure2siteoperator(A, M1, M2, chainsection::Tuple{Int,Int}, ::Nothing, ρ::Vector) =
    measure2siteoperator(A, M1, M2, chainsection, ρ)

measure(A::Vector, O::CdagCup; ρ=nothing, kwargs...) = measure(A, O, ρ)
function measure(A::Vector, O::CdagCup, ::Nothing)
    first = O.sites[1]
    last = O.sites[2]

    if first < last
        nearest = first
         farthest = last
    else
        nearest = last
        farthest = first
    end
    
    N = farthest - nearest + 1 

    ρ1 = ones(ComplexF64, 1, 1)
    ρ = ones(ComplexF64, 1, 1)
    
    expval = zeros(ComplexF64, N, N)

    for i=1:nearest-1
        ρ = rhoAAstar(ρ, A[i])
    end

    for i=nearest:farthest
        v = rhoAOAstar(ρ, A[i], Adagup*Aup, nothing)
        expval[i-nearest+1,i-nearest+1] = v
        ρ1 = rhoAOAstar(ρ, A[i], Adagup*parity)
        for j=i+1:farthest
            v = rhoAOAstar(ρ1, A[j], Aup, nothing)
            expval[i-nearest+1, j-nearest+1] = v
            expval[j-nearest+1, i-nearest+1] = conj(v)
            ρ1 = rhoAOAstar(ρ1, A[j], parity)
        end
        ρ = rhoAAstar(ρ, A[i])
    end
    if nearest != first
        expval = reverse(reverse(expval, dims=1), dims=2)
    end
    return expval
end
function measure(A::Vector, O::CdagCup, ρ::Vector)
    first = O.sites[1]
    last = O.sites[2]

    if first < last
        nearest = first
         farthest = last
    else
        nearest = last
        farthest = first
    end
    
    N = farthest - nearest + 1 

    ρ1 = ones(ComplexF64, 1, 1)
    
    expval = zeros(ComplexF64, N, N)

    for i=nearest:farthest
        v = rhoAOAstar(ρ[i], A[i], Adagup*Aup, nothing)
        expval[i-nearest+1,i-nearest+1] = v
        ρ1 = rhoAOAstar(ρ[i], A[i], Adagup*parity)
        for j=i+1:farthest
            v = rhoAOAstar(ρ1, A[j], Aup, nothing)
            expval[i-nearest+1, j-nearest+1] = v
            expval[j-nearest+1, i-nearest+1] = conj(v)
            ρ1 = rhoAOAstar(ρ1, A[j], parity)
        end
    end
    if nearest != first
        expval = reverse(reverse(expval, dims=1), dims=2)
    end
    return expval
end

measure(A::Vector, O::CdagCdn; ρ=nothing) = measure(A, O, ρ)
function measure(A::Vector, O::CdagCdn, ::Nothing)
    first = O.sites[1]
    last = O.sites[2]

    if first < last
        nearest = first
         farthest = last
    else
        nearest = last
        farthest = first
    end
    
    N = farthest - nearest + 1 

    ρ1 = ones(ComplexF64, 1, 1)
    ρ = ones(ComplexF64, 1, 1)
    
    expval = zeros(ComplexF64, N, N)

    for i=1:nearest-1
        ρ = rhoAAstar(ρ, A[i])
    end

    for i=nearest:farthest
        v = rhoAOAstar(ρ, A[i], Adagdn*Adn, nothing)
        expval[i-nearest+1,i-nearest+1] = v
        ρ1 = rhoAOAstar(ρ, A[i], Adagdn)
        for j=i+1:farthest
            v = rhoAOAstar(ρ1, A[j], parity*Adn, nothing)
            expval[i-nearest+1, j-nearest+1] = v
            expval[j-nearest+1, i-nearest+1] = conj(v)
            ρ1 = rhoAOAstar(ρ1, A[j], parity)
        end
        ρ = rhoAAstar(ρ, A[i])
    end
    if nearest != first
        expval = reverse(reverse(expval, dims=1), dims=2)
    end
    return expval
end
function measure(A::Vector, O::CdagCdn, ρ::Vector)
    first = O.sites[1]
    last = O.sites[2]

    if first < last
        nearest = first
         farthest = last
    else
        nearest = last
        farthest = first
    end
    
    N = farthest - nearest + 1 

    ρ1 = ones(ComplexF64, 1, 1)
    
    expval = zeros(ComplexF64, N, N)

    for i=nearest:farthest
        v = rhoAOAstar(ρ[i], A[i], Adagdn*Adn, nothing)
        expval[i-nearest+1,i-nearest+1] = v
        ρ1 = rhoAOAstar(ρ[i], A[i], Adagdn)
        for j=i+1:farthest
            v = rhoAOAstar(ρ1, A[j], parity*Adn, nothing)
            expval[i-nearest+1, j-nearest+1] = v
            expval[j-nearest+1, i-nearest+1] = conj(v)
            ρ1 = rhoAOAstar(ρ1, A[j], parity)
        end
    end
    if nearest != first
        expval = reverse(reverse(expval, dims=1), dims=2)
    end
    return expval
end

"""
     measure(A::Vector, Os::Vector; kwargs...)

Caculate expectation of Os on MPS A.

"""


function measure(A::Vector, Os::Vector; kwargs...)
    numobs = length(Os)
    numobs==0 && return Any[]
    N = max(reach.((A,), Os)...)
    res = Vector{Any}(undef, numobs)
    ρ = leftcontractmps(A, N)
    for (k, obs) in enumerate(Os)
        res[k] = measure(A, obs; ρ=ρ, kwargs...)
    end
    return res
end

function leftcontractmps(A, N::Int=length(A))
    ρ = Vector{Any}(undef, N)
    ρ[1] = ones(ComplexF64, 1, 1)
    for i=2:N
        ρ[i] = rhoAAstar(ρ[i-1], A[i-1])
    end
    return ρ
end
function leftcontractmps(A, O::Vector, N::Int=length(A))
    ρ = Vector{Any}(undef, N)
    ρ[1] = ones(ComplexF64, 1, 1)
    numops = length(O)
    for i=2:numops
        ρ[i] = rhoAOAstar(ρ[i-1], A[i-1], O[i-1])
    end
    for i=numops+1:N
        ρ[i] = rhoAAstar(ρ[i-1], A[i-1])
    end
    return ρ
end
