function ACOAC(AC::Array{T1, 3}, O::Array{T2, 2}) where {T1,T2}
    @tensor v = scalar(conj(AC[a,b,s']) * O[s',s] * AC[a,b,s])
end

function contractC!(A::Array{T1,2}, C::Array{T2,2}, dir::Int) where {T1,T2}
    return @tensor A[a0,s] := A[a0',s] * C[a0,a0']
end
function contractC!(A::Array{T1,3}, C::Array{T2,2}, dir::Int) where {T1,T2}
    dir==1 && return @tensor A[a0,a1,s] := A[a0',a1,s] * C[a0,a0']
    dir==2 && return @tensor A[a0,a1,s] := A[a0,a1',s] * C[a1,a1']
    throw("invalid argument dir=$dir")
end
function contractC!(A::Array{T1,4}, C::Array{T2,2}, dir::Int) where {T1,T2}
    dir==1 && return @tensor A[a0,a1,a2,s] := A[a0',a1,a2,s] * C[a0,a0']
    dir==2 && return @tensor A[a0,a1,a2,s] := A[a0,a1',a2,s] * C[a1,a1']
    dir==3 && return @tensor A[a0,a1,a2,s] := A[a0,a1,a2',s] * C[a2,a2']
    throw("invalid argument dir=$dir")
end
function contractC!(A::Array{T1,5}, C::Array{T2,2}, dir::Int) where {T1,T2}
    dir==1 && return @tensor A[a0,a1,a2,a3,s] := A[a0',a1,a2,a3,s] * C[a0,a0']
    dir==2 && return @tensor A[a0,a1,a2,a3,s] := A[a0,a1',a2,a3,s] * C[a1,a1']
    dir==3 && return @tensor A[a0,a1,a2,a3,s] := A[a0,a1,a2',a3,s] * C[a2,a2']
    dir==4 && return @tensor A[a0,a1,a2,a3,s] := A[a0,a1,a2,a3',s] * C[a3,a3']
    throw("invalid argument dir=$dir")
end
function contractC!(A::Array{T1,6}, C::Array{T2,2}, dir::Int) where {T1,T2}
    dir==1 && return @tensor A[a0,a1,a2,a3,a4,s] := A[a0',a1,a2,a3,a4,s] * C[a0,a0']
    dir==2 && return @tensor A[a0,a1,a2,a3,a4,s] := A[a0,a1',a2,a3,a4,s] * C[a1,a1']
    dir==3 && return @tensor A[a0,a1,a2,a3,a4,s] := A[a0,a1,a2',a3,a4,s] * C[a2,a2']
    dir==4 && return @tensor A[a0,a1,a2,a3,a4,s] := A[a0,a1,a2,a3',a4,s] * C[a3,a3']
    dir==5 && return @tensor A[a0,a1,a2,a3,a4,s] := A[a0,a1,a2,a3,a4',s] * C[a4,a4']
    throw("invalid argument dir=$dir")
end

function rhoAAstar(ρ::Array{T1,2}, A::Array{T2,2}, indir::Int) where {T1,T2}
    @tensoropt ρO = scalar(ρ[a0,b0] * conj(A[a0,s]) * A[b0,s])
end
function rhoAAstar(ρ::Array{T1,2}, A::Array{T2,2}) where {T1,T2}
    @tensoropt ρO = scalar(ρ[a0,b0] * conj(A[a0,s]) * A[b0,s])
end
function rhoAAstar(ρ::Array{T1,2}, A::Array{T2,3}, indir::Int, outdir::Int) where {T1,T2}
    indir==1 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a0,a,s]) * A[b0,b,s]
    indir==2 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a,a0,s]) * A[b,b0,s]
    throw("invalid arguments indir=$indir, outdir=$outdir")
end
function rhoAAstar(ρ::Array{T1,2}, A::Array{T2,3}) where {T1,T2}
    return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a0,a,s]) * A[b0,b,s]
end

function rhoAAstar(ρ::Array{T1,2}, A::Array{T2,4}, indir::Int, outdir::Int) where {T1,T2}
    indir==1 && outdir==2 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a0,a,c0,s]) * A[b0,b,c0,s]
    indir==1 && outdir==3 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a0,c0,a,s]) * A[b0,c0,b,s]

    indir==2 && outdir==1 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a,a0,c0,s]) * A[b,b0,c0,s]
    indir==2 && outdir==3 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,a0,a,s]) * A[c0,b0,b,s]

    indir==3 && outdir==1 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a,c0,a0,s]) * A[b,c0,b0,s]
    indir==3 && outdir==2 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,a,a0,s]) * A[c0,b,b0,s]
    throw("invalid arguments indir=$indir, outdir=$outdir")
end
function rhoAAstar(ρ::Array{T1,2}, A::Array{T2,5}, indir::Int, outdir::Int) where {T1,T2}
    indir==1 && outdir==2 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a0,a,c0,c1,s]) * A[b0,b,c0,c1,s]
    indir==1 && outdir==3 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a0,c0,a,c1,s]) * A[b0,c0,b,c1,s]
    indir==1 && outdir==4 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a0,c0,c1,a,s]) * A[b0,c0,c1,b,s]

    indir==2 && outdir==1 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a,a0,c0,c1,s]) * A[b,b0,c0,c1,s]
    indir==2 && outdir==3 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,a0,a,c1,s]) * A[c0,b0,b,c1,s]
    indir==2 && outdir==4 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,a0,c1,a,s]) * A[c0,b0,c1,b,s]

    indir==3 && outdir==1 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a,c0,a0,c1,s]) * A[b,c0,b0,c1,s]
    indir==3 && outdir==2 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,a,a0,c1,s]) * A[c0,b,b0,c1,s]
    indir==3 && outdir==4 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,c1,a0,a,s]) * A[c0,c1,b0,b,s]

    indir==4 && outdir==1 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a,c0,c1,a0,s]) * A[b,c0,c1,b0,s]
    indir==4 && outdir==2 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,a,c1,a0,s]) * A[c0,b,c1,b0,s]
    indir==4 && outdir==3 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,c1,a,a0,s]) * A[c0,c1,b,b0,s]
    throw("invalid arguments indir=$indir, outdir=$outdir")
end

function rhoAOAstar(ρ::Array{T1,2}, A::Array{T2,2}, O::Array{T3,2}, indir::Int, ::Nothing) where {T1,T2,T3}
    @tensoropt ρO = scalar(ρ[a0,b0] * conj(A[a0,s']) * O[s',s] * A[b0,s])
end
function rhoAOAstar(ρ::Array{T1,2}, A::Array{T2,2}, O::Array{T3,2}, ::Nothing) where {T1,T2,T3}
    @tensoropt ρO = scalar(ρ[a0,b0] * conj(A[a0,s']) * O[s',s] * A[b0,s])
end
function rhoAOAstar(ρ::Array{T1,2}, A::Array{T2,3}, O::Array{T3,2}, indir::Int, ::Nothing) where {T1,T2,T3}
    indir==1 && return @tensoropt ρO = scalar(ρ[a0,b0] * conj(A[a0,c0,s']) * O[s',s] * A[b0,c0,s])
    indir==2 && return @tensoropt ρO = scalar(ρ[a0,b0] * conj(A[c0,a0,s']) * O[s',s] * A[c0,b0,s])
end
function rhoAOAstar(ρ::Array{T1,2}, A::Array{T2,3}, O::Array{T3,2}, ::Nothing) where {T1,T2,T3}
    return @tensoropt ρO = scalar(ρ[a0,b0] * conj(A[a0,c0,s']) * O[s',s] * A[b0,c0,s])
end

function rhoAOAstar(ρ::Array{T1,2}, A::Array{T2,4}, O::Array{T3,2}, indir::Int, ::Nothing) where {T1,T2,T3}
    indir==1 && return @tensoropt ρO = scalar(ρ[a0,b0] * conj(A[a0,c0,c1,s']) * O[s',s] * A[b0,c0,c1,s])
    indir==2 && return @tensoropt ρO = scalar(ρ[a0,b0] * conj(A[c0,a0,c1,s']) * O[s',s] * A[c0,b0,c1,s])
    indir==3 && return @tensoropt ρO = scalar(ρ[a0,b0] * conj(A[c0,c1,a0,s']) * O[s',s] * A[c0,c1,b0,s])
end
function rhoAOAstar(ρ::Array{T1,2}, A::Array{T2,5}, O::Array{T3,2}, indir::Int, ::Nothing) where {T1,T2,T3}
    indir==1 && return @tensoropt ρO = scalar(ρ[a0,b0] * conj(A[a0,c0,c1,c2,s']) * O[s',s] * A[b0,c0,c1,c2,s])
    indir==2 && return @tensoropt ρO = scalar(ρ[a0,b0] * conj(A[c0,a0,c1,c2,s']) * O[s',s] * A[c0,b0,c1,c2,s])
    indir==3 && return @tensoropt ρO = scalar(ρ[a0,b0] * conj(A[c0,c1,a0,c2,s']) * O[s',s] * A[c0,c1,b0,c2,s])
    indir==4 && return @tensoropt ρO = scalar(ρ[a0,b0] * conj(A[c0,c1,c2,a0,s']) * O[s',s] * A[c0,c1,c2,b0,s])
end

function rhoAOAstar(ρ::Array{T1,2}, A::Array{T2,3}, O::Array{T3,2}, indir::Int, outdir::Int) where {T1,T2,T3}
    indir==1 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a0,a,s']) * O[s',s] * A[b0,b,s]
    indir==2 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a,a0,s']) * O[s',s] * A[b,b0,s]
    throw("invalid arguments indir=$indir, outdir=$outdir")
end
function rhoAOAstar(ρ::Array{T1,2}, A::Array{T2,3}, O::Array{T3,2}) where {T1,T2,T3}
    return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a0,a,s']) * O[s',s] * A[b0,b,s]
end

function rhoAOAstar(ρ::Array{T1,2}, A::Array{T2,4}, O::Array{T3,2}, indir::Int, outdir::Int) where {T1,T2,T3}
    indir==1 && outdir==2 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a0,a,c0,s']) * O[s',s] * A[b0,b,c0,s]
    indir==1 && outdir==3 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a0,c0,a,s']) * O[s',s] * A[b0,c0,b,s]

    indir==2 && outdir==1 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a,a0,c0,s']) * O[s',s] * A[b,b0,c0,s]
    indir==2 && outdir==3 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,a0,a,s']) * O[s',s] * A[c0,b0,b,s]

    indir==3 && outdir==1 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a,c0,a0,s']) * O[s',s] * A[b,c0,b0,s]
    indir==3 && outdir==2 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,a,a0,s']) * O[s',s] * A[c0,b,b0,s]
    throw("invalid arguments indir=$indir, outdir=$outdir")
end
function rhoAOAstar(ρ::Array{T1,2}, A::Array{T2,5}, O::Array{T3,2}, indir::Int, outdir::Int) where {T1,T2,T3}
    indir==1 && outdir==2 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a0,a,c0,c1,s']) * O[s',s] * A[b0,b,c0,c1,s]
    indir==1 && outdir==3 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a0,c0,a,c1,s']) * O[s',s] * A[b0,c0,b,c1,s]
    indir==1 && outdir==4 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a0,c0,c1,a,s']) * O[s',s] * A[b0,c0,c1,b,s]

    indir==2 && outdir==1 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a,a0,c0,c1,s']) * O[s',s] * A[b,b0,c0,c1,s]
    indir==2 && outdir==3 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,a0,a,c1,s']) * O[s',s] * A[c0,b0,b,c1,s]
    indir==2 && outdir==4 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,a0,c1,a,s']) * O[s',s] * A[c0,b0,c1,b,s]

    indir==3 && outdir==1 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a,c0,a0,c1,s']) * O[s',s] * A[b,c0,b0,c1,s]
    indir==3 && outdir==2 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,a,a0,c1,s']) * O[s',s] * A[c0,b,b0,c1,s]
    indir==3 && outdir==4 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,c1,a0,a,s']) * O[s',s] * A[c0,c1,b0,b,s]

    indir==4 && outdir==1 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[a,c0,c1,a0,s']) * O[s',s] * A[b,c0,c1,b0,s]
    indir==4 && outdir==2 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,a,c1,a0,s']) * O[s',s] * A[c0,b,c1,b0,s]
    indir==4 && outdir==3 && return @tensoropt ρO[a,b] := ρ[a0,b0] * conj(A[c0,c1,a,a0,s']) * O[s',s] * A[c0,c1,b,b0,s]
    throw("invalid arguments indir=$indir, outdir=$outdir")
end

#dir gives the open index, where dir=1 is the first child
#Fs must be given in the correct order
#left = direction of parent
#right = direction of children
function updateleftenv(A::Array{T1,2}, M::Array{T2,3}, dir::Int) where {T1,T2}
    @tensor F[a,b,c] := conj(A[a,s'])*M[b,s',s]*A[c,s]
end
function updaterightenv(A::Array{T1,2}, M::Array{T2,3}) where {T1,T2}
    @tensor F[a,b,c] := conj(A[a,s'])*M[b,s',s]*A[c,s]
end
function updateleftenv(A::Array{T1,3}, M::Array{T2,4}, FL) where {T1,T2}
    @tensor F[a,b,c] := FL[a0,b0,c0]*conj(A[a0,a,s'])*M[b0,b,s',s]*A[c0,c,s]
end
function updateleftenv(A::Array{T1,3}, M::Array{T2,4}, dir::Int, F0) where {T1,T2}
    @tensor F[a,b,c] := F0[a0,b0,c0]*conj(A[a0,a,s'])*M[b0,b,s',s]*A[c0,c,s]
end
function updaterightenv(A::Array{T1,3}, M::Array{T2,4}, FR) where {T1,T2}
    @tensor F[a,b,c] := FR[a0,b0,c0]*conj(A[a,a0,s'])*M[b,b0,s',s]*A[c,c0,s]
end

function updateleftenv(A::Array{T1,4}, M::Array{T2,5}, dir::Int, F0, F1) where {T1,T2} # dir in {1,2}
    dir==1 && (Aperm = [1,3,2,4]; Mperm = [1,3,2,4,5])
    dir==2 && (Aperm = [1,2,3,4]; Mperm = [1,2,3,4,5])
    At = permutedims(A, Aperm)
    Mt = permutedims(M, Mperm)
    @tensoropt !(b,b0,b1) F[a,b,c] := F0[a0,b0,c0]*F1[a1,b1,c1]*conj(At[a0,a1,a,s'])*Mt[b0,b1,b,s',s]*At[c0,c1,c,s]
end
function updaterightenv(A::Array{T1,4}, M::Array{T2,5}, F1, F2) where {T1,T2}
    @tensoropt !(b,b1,b2) F[a,b,c] := F1[a1,b1,c1]*F2[a2,b2,c2]*conj(A[a,a1,a2,s'])*M[b,b1,b2,s',s]*A[c,c1,c2,s]
end
function updateleftenv(A::Array{T1,5}, M::Array{T2,6}, dir::Int, F0, F1, F2) where {T1,T2} # dir in {1,2,3}
    dir==1 && (Aperm = [1,3,4,2,5]; Mperm = [1,3,4,2,5,6])
    dir==2 && (Aperm = [1,2,4,3,5]; Mperm = [1,2,4,3,5,6])
    dir==3 && (Aperm = [1,2,3,4,5]; Mperm = [1,2,3,4,5,6])
    At = permutedims(A, Aperm)
    Mt = permutedims(M, Mperm)
    @tensoropt !(b,b0,b1,b2) F[a,b,c] := F0[a0,b0,c0]*F1[a1,b1,c1]*F2[a2,b2,c2]*conj(At[a0,a1,a2,a,s'])*Mt[b0,b1,b2,b,s',s]*At[c0,c1,c2,c,s]
end
function updaterightenv(A::Array{T1,5}, M::Array{T2,6}, F1, F2, F3) where {T1,T2}
    @tensoropt !(b,b1,b2,b3) F[a,b,c] := F1[a1,b1,c1]*F2[a2,b2,c2]*F3[a3,b3,c3]*conj(A[a,a1,a2,a3,s'])*M[b,b1,b2,b3,s',s]*A[c,c1,c2,c3,s]
end

#dir gives the open index, where dir=1 is the parent
function updateenv(A, M, dir, F1)
    @tensor F[a,b,c] := conj(A[a,s'])*M[b,s',s]*A[c,s]
end
function updateenv(A, M, dir, F1, F2)
    if dir==1
        updateenv1(A, M, F2)
    else
        updateenv2(A, M, F1)
    end
end
function updateenv(A, M, dir, F1, F2, F3)
    if dir==1
        updateenv1(A, M, F2, F3)
    elseif dir==2
        updateenv2(A, M, F1, F3)
    else
        updateenv3(A, M, F1, F2)
    end
end
function updateenv(A, M, dir, F1, F2, F3, F4)
    if dir==1
        updateenv1(A, M, F2, F3, F4)
    elseif dir==2
        updateenv2(A, M, F1, F3, F4)
    elseif dir==3
        updateenv3(A, M, F1, F2, F4)
    else
        updateenv4(A, M, F1, F2, F3)
    end
end
function updateenv1(A, M, F1)
    @tensor F[a,b,c] := F1[a1,b1,c1]*conj(A[a,a1,s'])*M[b,b1,s',s]*A[c,c1,s]
end
function updateenv2(A, M, F1)
    @tensor F[a,b,c] := F1[a1,b1,c1]*conj(A[a1,a,s'])*M[b1,b,s',s]*A[c1,c,s]
end
function updateenv1(A, M, F1, F2)
    @tensor F[a,b,c] := F1[a1,b1,c1]*F2[a2,b2,c2]*conj(A[a,a1,a2,s'])*M[b,b1,b2,s',s]*A[c,c1,c2,s]
end
function updateenv2(A, M, F1, F2)
    @tensor F[a,b,c] := F1[a1,b1,c1]*F2[a2,b2,c2]*conj(A[a1,a,a2,s'])*M[b1,b,b2,s',s]*A[c1,c,c2,s]
end
function updateenv3(A, M, F1, F2)
    @tensor F[a,b,c] := F1[a1,b1,c1]*F2[a2,b2,c2]*conj(A[a1,a2,a,s'])*M[b1,b2,b,s',s]*A[c1,c2,c,s]
end
function updateenv1(A, M, F1, F2, F3)
    @tensor F[a,b,c] := F1[a1,b1,c1]*F2[a2,b2,c2]*F3[a3,b3,c3]*conj(A[a,a1,a2,a3,s'])*M[b,b1,b2,b3,s',s]*A[c,c1,c2,c3,s]
end
function updateenv2(A, M, F1, F2, F3)
    @tensor F[a,b,c] := F1[a1,b1,c1]*F2[a2,b2,c2]*F3[a3,b3,c3]*conj(A[a1,a,a2,a3,s'])*M[b1,b,b2,b3,s',s]*A[c1,c,c2,c3,s]
end
function updateenv3(A, M, F1, F2, F3)
    @tensor F[a,b,c] := F1[a1,b1,c1]*F2[a2,b2,c2]*F3[a3,b3,c3]*conj(A[a1,a2,a,a3,s'])*M[b1,b2,b,b3,s',s]*A[c1,c2,c,c3,s]
end
function updateenv4(A, M, F1, F2, F3)
    @tensor F[a,b,c] := F1[a1,b1,c1]*F2[a2,b2,c2]*F3[a3,b3,c3]*conj(A[a1,a2,a3,a,s'])*M[b1,b2,b3,b,s',s]*A[c1,c2,c3,c,s]
end

function applyH2(AA, H1, H2, F1, F2)
    @tensoropt !(b1,b,b2) HAA[a1,s1,a2,s2] := F1[a1,b1,c1]*AA[c1,s1',c2,s2']*H1[b1,b,s1,s1']*H2[b,b2,s2,s2']*F2[a2,b2,c2]
end
function applyH1(AC, M, F)
    @tensoropt !(b) HAC[a,s'] := F[a,b,c]*AC[c,s]*M[b,s',s]
end
function applyH1(AC, M, F0, F1)
    @tensoropt !(b0,b1) HAC[a0,a1,s'] := F0[a0,b0,c0]*F1[a1,b1,c1]*AC[c0,c1,s]*M[b0,b1,s',s]
end
function applyH1(AC, M, F0, F1, F2)
    @tensoropt !(b0,b1,b2) HAC[a0,a1,a2,s'] := F0[a0,b0,c0]*F1[a1,b1,c1]*F2[a2,b2,c2]*AC[c0,c1,c2,s]*M[b0,b1,b2,s',s]
end
function applyH1(AC, M, F0, F1, F2, F3)
    @tensoropt !(b0,b1,b2,b3) HAC[a0,a1,a2,a3,s'] := F0[a0,b0,c0]*F1[a1,b1,c1]*F2[a2,b2,c2]*F3[a3,b3,c3]*AC[c0,c1,c2,c3,s]*M[b0,b1,b2,b3,s',s]
end
function applyH0(C, F0, F1)
    @tensor HC[α,β] := F0[α,a,α']*C[α',β']*F1[β,a,β']
end

#sets the right/left bond dimension of A, ie will truncate if Dnew is smaller than the current bond dimension and zero pad if it's larger
function setrightbond(A, Dnew::Int)
    Dl, Dr, d  = size(A)
    a = fill!(similar(A, Dl, Dnew, d), 0)
    if Dr > Dnew
        a = A[:,1:Dnew,:]
    else
        a[:,1:Dr,:] = A
    end
    return a
end
function setleftbond(A, Dnew::Int)
    Dl, Dr, d  = size(A)
    a = fill!(similar(A, Dnew, Dr, d), 0)
    if Dl > Dnew
        a = A[1:Dnew,:,:]
    else
        a[1:Dl,:,:] = A
    end
    return a
end

function setbond(A::Array{T, 2}, Dnew::Int) where T
    Dold, d = size(A)
    a = fill!(similar(A, Dnew, d), 0)
    Dold > Dnew ? (D = Dnew) : (D = Dold)
    a[1:D,:] = A[1:D,:]
    return a
end
function setbond(A::Array{T, 2}) where T
    return A
end
function setbond(A::Array{T, 3}, D1new::Int, D2new::Int) where T
    D1old, D2old, d = size(A)
    a = fill!(similar(A, D1new, D2new, d), 0)
    D1old > D1new ? (D1 = D1new) : (D1 = D1old)
    D2old > D2new ? (D2 = D2new) : (D2 = D2old)
    a[1:D1,1:D2,:] = A[1:D1,1:D2,:]
    return a
end
function setbond(A::Array{T, 3}, D2new::Int) where T
    D1old, D2old, d = size(A)
    a = fill!(similar(A, D1old, D2new, d), 0)
    D2old > D2new ? (D2 = D2new) : (D2 = D2old)
    a[:,1:D2,:] = A[:,1:D2,:]
    return a
end
function setbond(A::Array{T, 4}, D1new::Int, D2new::Int, D3new::Int) where T
    D1old, D2old, D3old, d = size(A)
    a = fill!(similar(A, D1new, D2new, D3new, d), 0)
    D1old > D1new ? (D1 = D1new) : (D1 = D1old)
    D2old > D2new ? (D2 = D2new) : (D2 = D2old)
    D3old > D3new ? (D3 = D3new) : (D3 = D3old)
    a[1:D1,1:D2,1:D3,:] = A[1:D1,1:D2,1:D3,:]
    return a
end
function setbond(A::Array{T, 4}, D2new::Int, D3new::Int) where T
    D1old, D2old, D3old, d = size(A)
    a = fill!(similar(A, D1old, D2new, D3new, d), 0)
    D2old > D2new ? (D2 = D2new) : (D2 = D2old)
    D3old > D3new ? (D3 = D3new) : (D3 = D3old)
    a[:,1:D2,1:D3,:] = A[:,1:D2,1:D3,:]
    return a
end
function setbond(A::Array{T, 5}, D1new::Int, D2new::Int, D3new::Int, D4new::Int) where T
    D1old, D2old, D3old, D4old, d = size(A)
    a = fill!(similar(A, D1new, D2new, D3new, D4new, d), 0)
    D1old > D1new ? (D1 = D1new) : (D1 = D1old)
    D2old > D2new ? (D2 = D2new) : (D2 = D2old)
    D3old > D3new ? (D3 = D3new) : (D3 = D3old)
    D4old > D4new ? (D4 = D4new) : (D4 = D4old)
    a[1:D1,1:D2,1:D3,1:D4,:] = A[1:D1,1:D2,1:D3,1:D4,:]
    return a
end
function setbond(A::Array{T, 5}, D2new::Int, D3new::Int, D4new::Int) where T
    D1old, D2old, D3old, D4old, d = size(A)
    a = fill!(similar(A, D1old, D2new, D3new, D4new, d), 0)
    D2old > D2new ? (D2 = D2new) : (D2 = D2old)
    D3old > D3new ? (D3 = D3new) : (D3 = D3old)
    D4old > D4new ? (D4 = D4new) : (D4 = D4old)
    a[:,1:D2,1:D3,1:D4,:] = A[:,1:D2,1:D3,1:D4,:]
    return a
end

truncAR(A::Array{T, 2}, Dnew) where T = A[1:Dnew,:]
truncAR(A::Array{T, 3}, Dnew) where T = A[1:Dnew,:,:]
truncAR(A::Array{T, 4}, Dnew) where T = A[1:Dnew,:,:,:]
truncAR(A::Array{T, 5}, Dnew) where T = A[1:Dnew,:,:,:,:]
truncF(F, D) = F[1:D,:,1:D]
truncF2(F, D) = F[:,:,1:D]

function evolveAC2(dt::Float64, A1, A2, M1, M2, FL, FR, energy=false; kwargs...)
    @tensor AA[a,sa,b,sb] := A1[a,c,sa] * A2[c,b,sb]
    AAnew, info = exponentiate(x->applyH2(x, M1, M2, FL, FR), -im*dt, AA; ishermitian = true, kwargs...)

    if energy
        E = real(dot(AAnew, applyH2(AAnew, M1, M2, FL, FR)))
        return AAnew, (E, info)
    end
    return AAnew, info
end

function evolveAC(dt::Float64, AC, M, FL, FR, energy=false; kwargs...)
    Dlnew, w, Dl = size(FL)
    Drnew, w, Dr = size(FR)

    AC = setrightbond(AC, Drnew)
    AC = setleftbond(AC, Dlnew)

    AC, info = exponentiate(x->applyH1(x, M, FL, FR), -im*dt, AC; ishermitian = true, kwargs...)

    if energy
        E = real(dot(AC, applyH1(AC, M, FL, FR)))
        return AC, (E, info)
    end
    return AC, info
end

function evolveC(dt::Float64, C, FL, FR, energy=false; kwargs...)
        
    C, info = exponentiate(x->applyH0(x, FL, FR), im*dt, C; ishermitian = true, kwargs...)

    if energy
        E = real(dot(C, applyH0(C, FL, FR)))
        return C, (E, info)
    end
    return C, info
end

# @optimalcontractiontree (a0=>x,a1=>x,c0=>x,c1=>x,s=>x,s'=>x,b0=>x,b1=>x) HAC[a0,a1,s'] := F0[a0,b0,c0]*F1[a1,b1,c1]*AC[c0,c1,s]*M[b0,b1,s',s]

# @tensoropt (a0=>x,a1=>x,c0=>x,c1=>x,s=>x,s'=>x,b0=>x,b1=>x) HAC[a0,a1,s'] := F0[a0,b0,c0]*F1[a1,b1,c1]*AC[c0,c1,s]*M[b0,b1,s',s]

# @optimalcontractiontree (a0=,a1,c0,c1,s,s',b0,b1) HAC[a0,a1,s'] := F0[a0,b0,c0]*F1[a1,b1,c1]*AC[c0,c1,s]*M[b0,b1,s',s]
# @optimalcontractiontree (α, β, α', β') HC[α,β] := F0[α,a,α']*C[α',β']*F1[β,a,β']

function LA_M_A(LA, M::Array{T1, 3}, A::Array{T2, 2}) where {T1,T2}
    @tensoropt !(b,b0) L[a0,s] := LA[a0,b0,c0] * M[b0,s,s'] * A[c0,s']
end
function LA_M_A(LA, M::Array{T1, 4}, A::Array{T2, 3}) where {T1,T2}
    @tensoropt !(b,b0) L[a0,s,b,c] := LA[a0,b0,c0] * M[b0,b,s,s'] * A[c0,c,s']
end
LA_M_A(LA, M::Array{T1, 4}, A::Array{T2, 3}, dir::Int) where {T1,T2} = LA_M_A(LA, M, A)
function LA_M_A(LA, M::Array{T1, 5}, A::Array{T2, 4}, dir::Int, F1) where {T1,T2}
    dir==2 && (Aperm = [1,3,2,4]; Mperm = [1,3,2,4,5])
    dir==3 && (Aperm = [1,2,3,4]; Mperm = [1,2,3,4,5])
    At = permutedims(A, Aperm)
    Mt = permutedims(M, Mperm)
    @tensoropt !(b,b0,b1) L[a0,a1,s,b,c] := LA[a0,b0,c0] * Mt[b0,b1,b,s,s'] * At[c0,c1,c,s'] * F1[a1,b1,c1]
end
function LA_M_A(LA, M::Array{T1, 6}, A::Array{T2, 5}, dir::Int, F1, F2) where {T1,T2}
    dir==2 && (Aperm = [1,3,4,2,5]; Mperm = [1,3,4,2,5,6])
    dir==3 && (Aperm = [1,2,4,3,5]; Mperm = [1,2,4,3,5,6])
    dir==4 && (Aperm = [1,2,3,4,5]; Mperm = [1,2,3,4,5,6])
    At = permutedims(A, Aperm)
    Mt = permutedims(M, Mperm)
    @tensoropt !(b,b0,b1,b2) L[a0,a1,a2,s,b,c] := LA[a0,b0,c0] * Mt[b0,b1,b2,b,s,s'] * At[c0,c1,c2,c,s'] * F1[a1,b1,c1] * F2[a2,b2,c2]
end

function L_FR(L::Array{T, 4}, FR) where T
    @tensor PA[a1,s,a] := L[a1,s,b,c] * FR[a,b,c]    
end
function L_FR(L::Array{T, 5}, FR) where T
    @tensor PA[a1,a2,s,a] := L[a1,a2,s,b,c] * FR[a,b,c]
end
function L_FR(L::Array{T, 6}, FR) where T
    @tensor PA[a1,a2,a3,s,a] := L[a1,a2,a3,s,b,c] * FR[a,b,c]
end

function L_AL(L::Array{T1, 4}, AL::Array{T2, 3}, Dmax::Int, dir::Int) where {T1,T2}
    @tensor LA[a,b,c] := L[a0,s,b,c] * (AL[:,1:Dmax,:])[a0,a,s]
end
function L_AL(L::Array{T1, 5}, AL::Array{T2, 4}, Dmax::Int, dir::Int) where {T1,T2}
    dir==2 && (Aperm = [1,3,2,4])
    dir==3 && (Aperm = [1,2,3,4])
    ALt = permutedims(AL, Aperm)
    @tensor LA[a,b,c] := L[a0,a1,s,b,c] * (ALt[:,:,1:Dmax,:])[a0,a1,a,s]
end
function L_AL(L::Array{T1, 6}, AL::Array{T2, 5}, Dmax::Int, dir::Int) where {T1,T2}
    dir==2 && (Aperm = [1,3,4,2,5])
    dir==3 && (Aperm = [1,2,4,3,5])
    dir==4 && (Aperm = [1,2,3,4,5])
    ALt = permutedims(AL, Aperm) 
    @tensor LA[a,b,c] := L[a0,a1,a2,s,b,c] * (ALt[:,:,:,1:Dmax,:])[a0,a1,a2,a,s]
end

function LA_FR(LA, FR)
    @tensor PC[a,a'] := LA[a,b',c'] * FR[a',b',c']
end

import LinearAlgebra: transpose
function transpose(A::AbstractArray, dim1::Int, dim2::Int)
    nd=ndims(A)
    perm=collect(1:nd)
    perm[dim1]=dim2
    perm[dim2]=dim1
    permutedims(A, perm)
end
function QR(A::AbstractArray, i::Int)
    dims = [size(A)...]
    nd = length(dims)
    ds = collect(1:nd)
    AL, C = qr(reshape(permutedims(A, circshift(ds, -i)), :, dims[i]))
    AL = permutedims(reshape(Matrix(AL), circshift(dims, -i)...), circshift(ds, i))
    return AL, C
end
function QR_full(A::AbstractArray, i::Int)
    dims = [size(A)...]
    nd = length(dims)
    ds = collect(1:nd)
    AL, C = qr(reshape(permutedims(A, circshift(ds, -i)), :, dims[i]))
    AL = permutedims(reshape(AL*Matrix(I,size(AL)...), circshift(dims, -i)[1:end-1]..., :), circshift(ds, i))
    return AL, C
end

"""
    U, S, Vd = svdtrunc(A; truncdim = max(size(A)...), truncerr = 0.)
Perform a truncated SVD, with maximum number of singular values to keep equal to `truncdim`
or truncating any singular values smaller than `truncerr`. If both options are provided, the
smallest number of singular values will be kept.
Unlike the SVD in Julia, this returns matrix U, a diagonal matrix (not a vector) S, and
Vt such that A ≈ U * S * Vt
"""
function svdtrunc(A; truncdim = max(size(A)...), truncerr = 0., kwargs...)
    F = svd(A)
    d = min(truncdim, count(F.S .>= truncerr))
    return F.U[:,1:d], diagm(0=>F.S[1:d]), F.Vt[1:d, :]
end

