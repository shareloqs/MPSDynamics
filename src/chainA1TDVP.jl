function initrightenvs_a1tdvp(A, M)
    N = length(A)
    Acomps = Vector{Any}(undef, N)
    F = Vector{Any}(undef, N+2)
    F[1] = fill!(similar(M[1], (1,1,1)), 1)
    F[N+2] = fill!(similar(M[1], (1,1,1)), 1)

    for k=N:-1:1
        AR, ARcomp, C = tensor_lq_full(A[k])
        Acomps[k] = ARcomp
        F[k+1] = updaterightenv(AR, M[k], F[k+2])
    end
    return F, Acomps
end

function tdvpa1sweep!(dt2, A, M, Acomp=nothing, F=nothing; verbose=false, projerr=true, prec=10^-4, Dlim=5, kwargs...)
    dt = dt2/2
    N = length(A)

    if Acomp==nothing || F==nothing
        F, Acomp = initrightenvs_a1tdvp(A, M)
    end

    svs = Vector{Any}(undef, N-1)
    projerr && (h2 = measurempo(A, multiply(M,M)))

    projnorm = 0
    data = []

    AC = A[1]
    for k = 1:N-1
        ARaug, sv = expandAR(F[k], F[k+3], AC, A[k+1], M[k], M[k+1], Acomp[k+1], prec, Dlim)
        svs[k] = sv
        FR = updaterightenv(ARaug, M[k+1], F[k+3])
        
        Dold = size(AC,2)

        AC, info = evolveAC(dt, AC, M[k], F[k], FR, verbose; projerr=projerr, kwargs...)
        Dnew = size(AC,2)
        verbose && println("Sweep L->R: AC site $k, energy $(info[1])")
        verbose && Dnew!=Dold && println("*BondDimension $k-$(k+1) changed from $Dold to $Dnew")
        verbose && Dnew==Dold && println("*BondDimension $k-$(k+1) constant at $Dnew")

        projerr && (projnorm += info[2]^2)

        AL, C = tensor_qr(AC)
        A[k] = AL
        F[k+1] = updateleftenv(A[k], M[k], F[k])

        C, info = evolveC(dt, C, F[k+1], FR, verbose; projerr=projerr, kwargs...)
        verbose && println("Sweep L->R: C between site $k and $(k+1), energy $(info[1])")

        projerr && (projnorm -= info[2]^2)

        @tensor AC[a,c,s] := C[a,b] * ARaug[b,c,s]
    end
    k = N

    AC, info = evolveAC(dt2, AC, M[k], F[k], F[k+2], verbose; projerr=projerr, kwargs...)
    verbose && println("Sweep L->R: AC site $k, energy $(info[1])")

    projerr && (projnorm += info[2]^2)

    for k = N-1:-1:1
        AR, ARcomp, C = tensor_lq_full(AC)
        A[k+1] = AR
        Acomp[k+1] = ARcomp
        F[k+2] = updaterightenv(A[k+1], M[k+1], F[k+3])

        C, info = evolveC(dt, C, F[k+1], F[k+2], verbose; kwargs...)
        verbose && println("Sweep R->L: C between site $k and $(k+1), energy $(info[1])")

        @tensor AC[a,c,s] := A[k][a,b,s] * C[b,c]

        AC, info = evolveAC(dt, AC, M[k], F[k], F[k+2], verbose; kwargs...)
        verbose && println("Sweep R->L: AC site $k, energy $(info[1])")
    end
    A[1] = AC

    projerr && push!(data, ("err", h2 - projnorm))
    projerr && push!(data, ("projerr", projnorm))
    push!(data, ("sigvals",svs))
    push!(data, ("dims", [bonddims(A)...]))
    return A, Acomp, F, Dict(data)
end

function expandAR(F1, F2, AC, AR, M1, M2, Acomp, prec, Dlim)
    @tensoropt G[a1,s1,aa] := F1[a1,b1,c1]*M1[b1,b0,s1,s1']*M2[b0,b2,s2,s2']*AC[c1,c0,s1']*AR[c0,c2,s2']*Acomp[aa,a2,s2]*F2[a2,b2,c2]

    Dl, d1, Dcomp = size(G)
    Dold, Dr, d2 = size(AR)

    Dmax = Dl*d1
    
    F = svd(reshape(G, Dl*d1, Dcomp))

    # Dadd = min(count((F.S ./max(F.S...)) .> thresh), Dmax-Dold)
#    Dadd = min(count(F.S .> prec), Dmax-Dold, Dlim-Dold)
    Dadd = min(1, size(F.Vt,1), Dmax-Dold, Dlim-Dold)

    randVt = randisometry(ComplexF64, size(F.Vt)...)
    @tensor Aadd[a,c,s] := randVt[1:Dadd,:][a,b] * Acomp[b,c,s]
#    @tensor Aadd[a,c,s] := F.Vt[1:Dadd,:][a,b] * Acomp[b,c,s]
    #@tensor Aadd[a,c,s] := F.Vt[end-Dadd+1:end,:][a,b] * Acomp[b,c,s]

    Aaug = fill!(similar(Aadd, Dold+Dadd, Dr, d2), 0)
    Aaug[1:Dold, :, :] = AR
    Aaug[Dold+1:end, :, :] = Aadd
#    Aaug[Dold+1:end, :, :] = Acomp[1:Dadd,:,:]
    return Aaug, F.S
end
