using LinearAlgebra

""" 
    gauss(N,ab)

Gauss quadrature rule for `N` sites on an interval `ab`. Given a weight function w encoded by the `N`x2 array `ab` of the first `N` recurrence coefficients for the associated orthogonal
polynomials, the first column of `ab` containing the `N` alpha-coefficients and the second column the `N` beta-coefficients, the call xw = gauss(N,ab) generates the nodes and weights xw of the `N`-point Gauss quadrature rule for the weight function w.
The nodes, in increasing order, are stored in the first column, the n corresponding weights in the second column, of the `N`x2 array xw.
"""
function gauss(N,ab)
    N0 = size(ab,1)
    if N0<N
        error("input array ab too short")
    end
    J = zeros(N,N)
    for n=1:N
        J[n,n] = ab[n,1]
    end
    for n=2:N
      J[n,n-1] = sqrt(ab[n,2])
      J[n-1,n] = J[n,n-1]
    end
    F = eigen(J)
    V = F.vectors
    D = F.values
    I = sortperm(D)
    D = D[I]
    V = V[:,I]
    xw = hcat(D, ab[1,2]*V[1,:].^2)

    return xw
end
