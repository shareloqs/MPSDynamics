function tdvp2sweep!(dt2, A::Vector, M::Vector, F=nothing; verbose=false, kwargs...)
    N = length(A)
    dt = dt2/2
    F = initenvs(A, M, F)
    AC = A[1]
    for k=1:N-1
        AR = A[k+1]
        Dold = size(AC)[2]
        AA, info = evolveAC2(dt, AC, AR, M[k], M[k+1], F[k], F[k+3], verbose; kwargs...)

        verbose && println("Sweep L->R: updating sites $k and $(k+1), energy $(info[1])")

        Dl, dl, Dr, dr = size(AA)
        U, S, Vt = svdtrunc(reshape(AA, Dl*dl, Dr*dr); kwargs...)
        Dnew = size(S)[1]

        verbose && Dnew!=Dold && println("*BondDimension $k-$(k+1) changed from $Dold to $Dnew")
        
        AL = permutedims(reshape(U, Dl, dl, Dnew), [1,3,2])
        F[k+1] = updateleftenv(AL, M[k], F[k])
        A[k] = AL
        AC = reshape(S*Vt, Dnew, Dr, dr)
    end
    for k=N-1:-1:1
        AL = A[k]
        Dold = size(AC)[1]
        AA, info = evolveAC2(dt, AL, AC, M[k], M[k+1], F[k], F[k+3], verbose; kwargs...)

        verbose && println("Sweep R->L: updating sites $(k+1) and $k, energy $(info[1])")
        
        Dl, dl, Dr, dr = size(AA)
        U, S, Vt = svdtrunc(reshape(AA, Dl*dl, Dr*dr); kwargs...)
        Dnew = size(S)[1]

        verbose && Dnew!=Dold && println("*BondDimension $(k+1)-$k changed from $Dold to $Dnew")

        AR = reshape(Vt, Dnew, Dr, dr)
        F[k+2] = updaterightenv(AR, M[k+1], F[k+3])
        A[k+1] = AR
        AC = permutedims(reshape(U*S, Dl, dl, Dnew), [1,3,2])
    end
    A[1] = AC
    return A, F
end
