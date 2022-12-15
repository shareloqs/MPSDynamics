using SpecialFunctions # for Euler Gamma function

""" R_JACOBI Recurrence coefficients for monic Jacobi polynomials.

    ab=R_JACOBI(n,a,b) generates the first n recurrence
    coefficients for monic Jacobi polynomials with parameters
    a and b. These are orthogonal on [-1,1] relative to the
    weight function w(t)=(1-t)^a(1+t)^b. The n alpha-coefficients
    are stored in the first column, the n beta-coefficients in
    the second column, of the nx2 array ab. The call ab=
    R_JACOBI(n,a) is the same as ab=R_JACOBI(n,a,a) and
    ab=R_JACOBI(n) the same as ab=R_JACOBI(n,0,0).

    Supplied by Dirk Laurie, 6-22-1998; edited by Walter
    Gautschi, 4-4-2002.

    Written in Julia by Thibaut Lacroix 10/2020.
"""

function r_jacobi(N::Int64, a=0,b=a)

    if ((N<=0)||(a<=-1)||(b<=-1))
        error("parameter(s) out of range")
    end

    nu = (b-a)/(a+b+2)

    mu = 2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2)

    if N==1
        ab = [nu mu]
        return ab
    end

    N = N-1
    n = [1:N...]
    nab = 2*n .+ a .+ b
    A = (b^2-a^2)*ones(N)./(nab.*(nab.+2))
    pushfirst!(A, nu)
    n = [2:N...]
    nab = nab[n]
    B1 = 4*(a+1)*(b+1)/((a+b+2)^2*(a + b+ 3))
    B = 4*(n.+a).*(n.+b).*n.*(n.+a.+b)./((nab.^2).*(nab.+1).*(nab.-1))
    pushfirst!(B,B1)
    pushfirst!(B,mu)
    ab = hcat(A,B)

    return ab
end
