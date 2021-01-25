
function tdvp2sweep!(dt2, A::Vector, M::Vector, F=nothing; verbose=false, kwargs...)
    N = length(A)
    dt = dt2/2
    F = initenvs(A, M, F)
    AC = A[1]
    for k=1:N-1
        AR = A[k+1]
        AA, info = evolveAC2(dt, AC, AR, M[k], M[k+1], F[k], F[k+3])
        Dl, dl, Dr, dr = size(AA)
        U, S, Vt = svdtrunc(reshape(AA, Dl*dl, Dr*dr))
        Dnew = size(S)[1]
        AL = permutedims(reshape(U, Dl, dl, Dnew), [1,3,2])
        F[k+1] = updateleftenv(AL, M[k], F[k])
        A[k] = AL
        AC = reshape(S*Vt, Dnew, Dr, dr)
    end
    A[N] = AC
#    return A
    for k=N-1:-1:1
        AL = A[k]
        AA, info = evolveAC2(dt, AL, AC, M[k], M[k+1], F[k], F[k+3])
        Dl, dl, Dr, dr = size(AA)
        U, S, Vt = svdtrunc(reshape(AA, Dl*dl, Dr*dr))
        Dnew = size(S)[1]
        AR = reshape(Vt, Dnew, Dr, dr)
        F[k+2] = updaterightenv(AR, M[k+1], F[k+3])
        A[k+1] = AR
        AC = permutedims(reshape(U*S, Dl, dl, Dnew), [1,3,2])
#        println(size(AC))
    end
    A[1] = AC
    return A, F
end
