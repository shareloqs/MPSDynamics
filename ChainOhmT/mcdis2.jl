include("gauss.jl")
include("r_jacobi.jl")
include("lanczos.jl")
include("stieltjes.jl")

"""
MCDIS Multiple-component discretization procedure.

    This routine performs a sequence of discretizations of the
    given weight function (or measure), each discretization being
    followed by an application of the Stieltjes, or Lanczos,
    procedure to produce approximations to the desired recurrence
    coefficients. The fineness of the discretization is
    characterized by a discretization parameter N. The support of
    the continuous part of the weight function is decomposed into
    a given number mc of subintervals (some or all of which may
    be identical). The routine then applies to each subinterval
    an N-point quadrature rule to discretize the weight function
    on that subinterval. The discrete part of the weight function
    (if there is any) is added on to the discretized continuous
    weight function. The sequence of discretizations, if chosen
    judiciously, leads to convergence of the recurrence
    coefficients for the discretized measures to those of the
    given measure. If convergence to within a prescribed accuracy
    eps0 occurs before N reaches its maximum allowed value Nmax,
    then the value of N that yields convergence is output as
    Ncap, and so is the number of iterations, kount. If there is
    no convergence, the routine displays the message "Ncap
    exceeds Nmax in mcdis" prior to exiting.

    The choice between the Stieltjes and the Lanczos procedure is
    made by setting the parameter irout equal to 1 in the former,
    and different from 1, in the latter case.

    The details of the discretization are to be specified prior
    to calling the procedure. They are embodied in the following
    global parameters:

    mc     = the number of component intervals
    mp     = the number of points in the discrete part of the
             measure (mp=0 if there is none)
    iq     = a parameter to be set equal to 1, if the user
             provides his or her own quadrature routine, and
             different from 1 otherwise
    idelta = a parameter whose default value is 1, but is
             preferably set equal to 2, if iq=1 and the user
             provides Gauss-type quadrature routines

    The component intervals have to be specified (in the order
    left to right) by a global mcx2 array AB=[[a1 b1];[a2 b2];
    ...;[amc bmc]],  where for infinite extreme intervals a1=-Inf
    resp. bmc=Inf. The discrete spectrum (if mp>0) is similarly
    specified by a global mpx2 array DM=[[x1 y1];[x2 y2];...;
    ;[xmp ymp]] containing the abscissae and jumps.

    If the user provides his or her own quadrature routine
    "quadown", the routine mcdis must be called with the input
    parameter "quad" replaced by "@quadown", otherwise with
    "quad" replaced by "@quadgp", a general-purpose routine
    provided in the package. The quadrature routine must have
    the form

                             function xw=quad(N,i)

    where N is the number of nodes and i identifies the interval
    to which the routine is to be applied.

    The routine mcdis also applies to measures given originally
    in multi-component form.
"""

function mcdis(N,eps0,quad::Function,Nmax,idelta,mc,AB,wf,mp,irout)

    suc = true
    f = "Ncap exceeds Nmax in mcdis with irout = "
    if N<1
        error("Input variable N out of range")
    end
    iNcap = 1
    kount = -1
    ab = Array{Float64}(undef, N, 2)
    ab[:,2] = zeros(N,1)
    b = ones(N)
    Ncap = floor((2*N-1)/idelta)
    uv = []

    while any( abs.(ab[:,2]-b) > eps0*abs.(ab[:,2]) )
      b = ab[:,2]
      kount = kount+1;
      if kount>1
           iNcap = 2^(floor(kount/5))*N
      end
      Ncap = Int64(Ncap + iNcap)
      if Ncap>Nmax
        print(f,irout,"\n")
        suc=false
        return [ab,Ncap,kount,suc,uv]
      end
      mtNcap = mc*Ncap

      jacs = r_jacobi(Ncap,0,0.0)
      uv = gauss(Ncap,jacs)

      xwm = Array{Float64}(undef, mc*Ncap, 2)
      for i=1:mc
        im1tn = (i-1)*Ncap
        xw = quad(Ncap,i,uv,mc,AB,wf)
        xwm[im1tn+1:im1tn+Ncap,:] = xw
      end

      if mp!=0
          xwm[mtNcap+1:mtNcap+mp,:] = DM
      end
      if irout==1
        ab = stieltjes(N,xwm)
      else
        ab = lanczos(N,xwm)
      end
    end

    return [ab,Ncap,kount,suc, uv]
end
