""" LANCZOS Lanczos algorithm.

    Given the discrete inner product whose nodes are contained
    in the first column, and whose weights are contained in the
    second column, of the nx2 array xw, the call ab=LANCZOS(n,xw)
    generates the first n recurrence coefficients ab of the
    corresponding discrete orthogonal polynomials. The n alpha-
    coefficients are stored in the first column, the n beta-
    coefficients in the second column, of the nx2 array ab.

    The script is adapted from the routine RKPW in
    W.B. Gragg and W.J. Harrod, ``The numerically stable
    reconstruction of Jacobi matrices from spectral data'',
    Numer. Math. 44 (1984), 317-335.
"""
function lanczos(N,xw) #return ab
    Ncap = size(xw,1)
    if (N<=0 || N>Ncap)
        error("N out of range")
    end

    p0 = xw[:,1]
    p1 = zeros(Ncap,1)
    p1[1] = xw[1,2]

    for n=1:Ncap-1
      pn = xw[n+1,2]
      gam = 1
      sig = 0
      t = 0
      xlam = xw[n+1,1]

      for k=1:n+1
        rho = p1[k] + pn
        tmp = gam*rho
        tsig = sig

        if rho<=0
          gam = 1
          sig = 0
        else
          gam = p1[k]/rho
          sig = pn/rho
        end

        tk = sig*(p0[k]-xlam)-gam*t
        p0[k] = p0[k] - (tk-t)
        t = tk

        if sig<=0
          pn = tsig*p1[k]
        else
          pn = (t^2)/sig
        end

        tsig = sig
        p1[k] = tmp
      end
    end

    ab = [p0[1:N] p1[1:N]]

    return ab
end
