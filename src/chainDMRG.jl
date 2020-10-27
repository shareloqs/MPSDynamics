function dmrg1sweep!(A::Vector, M::Vector, F=nothing; verbose=false, kwargs...)
    N = length(A)

    if F == nothing
        F = Vector{Any}(undef, N+2)
        F[1] = fill!(similar(M[1], (1,1,1)), 1)
        F[N+2] = fill!(similar(M[1], (1,1,1)), 1)
        for k = N:-1:1
            F[k+1] = updaterightenv(A[k], M[k], F[k+2])
        end
    end

    AC = A[1]
    for k = 1:N-1
        Es, ACs, info = eigsolve(x->applyH1(x, M[k], F[k], F[k+2]), AC, 1, :SR; ishermitian = true, kwargs...)
        AC = ACs[1]
        E = Es[1]

        verbose && println("Sweep L->R: site $k -> energy $E")

        AL, C = QR(AC, 2)
        A[k] = AL
        F[k+1] = updateleftenv(A[k], M[k], F[k])

        @tensor AC[-1,-2,-3] := C[-1,1] * A[k+1][1,-2,-3]
    end
    k = N
    Es, ACs, info = eigsolve(x->applyH1(x, M[k], F[k], F[k+2]), AC, 1, :SR; ishermitian = true, kwargs...)
    AC = ACs[1]
    E = Es[1]
    verbose && println("Sweep L->R: site $k -> energy $E")
    for k = N-1:-1:1
        AR, C = QR(AC, 1)
        A[k+1] = AR
        F[k+2] = updaterightenv(A[k+1], M[k+1], F[k+3])

        @tensor AC[-1,-2,-3] := C[-2,2] * A[k][-1,2,-3]

        Es, ACs, info = eigsolve(x->applyH1(x, M[k], F[k], F[k+2]), AC, 1, :SR; ishermitian = true, kwargs...)
        AC = ACs[1]
        E = Es[1]

        verbose && println("Sweep R->L: site $k -> energy $E")
    end
    A[1] = AC
    return E, A, F
end
