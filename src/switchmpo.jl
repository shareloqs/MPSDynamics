import QuadGK: quadgk
import Jacobi: jacobi

function polybeta(t::Float64, n::Int, a::Array, b::Array, temp::Array)
    """
    polybeta recursively constructs the polynomials used to compute the coupling coefficients given the coefficients a and b
    this function is useful when working at finite temperature (beta != inf)
    """
          if n==-1
              return 0
          elseif n==0
              if length(temp)>=2
                  temp[2] = 1
              else
                  push!(temp,1)
              end

              return 1
          elseif n==1
              pn = (t - a[n])
              if length(temp) == 2
                  push!(temp,pn)
              elseif length(temp) == 1
                  push!(temp, 1)
                  push!(temp,pn)
              end

              return pn
          else
              if length(temp)<n+1 && temp[1] == t
                  pn = (t - a[n])*polybeta(t,n-1,a,b,temp) - b[n-1]*polybeta(t,n-2,a,b,temp) #P_{n}(t) = (t-a_{n-1})P_{n-1} - b_{n-1}P_{n-2}
                  push!(temp, pn)

                  return pn
              elseif length(temp)<n+1 && temp[1] != t
                  temp = [t]
                  pn = (t - a[n])*polybeta(t,n-1,a,b,temp) - b[n-1]*polybeta(t,n-2,a,b,temp) #P_{n}(t) = (t-a_{n-1})P_{n-1} - b_{n-1}P_{n-2}
                  push!(temp, pn)

                  return pn
              elseif length(temp) == n+1 && temp[1] == t
                   pn = (t - a[n])*temp[n+1] - b[n-1]*temp[n]
                   push!(temp,pn)

                   return pn
               elseif length(temp) == n+1 && temp[1] != t
                   temp = [t]
                   pn = (t - a[n])*polybeta(t,n-1,a,b,temp) - b[n-1]*polybeta(t,n-2,a,b,temp) #P_{n}(t) = (t-a_{n-1})P_{n-1} - b_{n-1}P_{n-2}
                   push!(temp, pn)

                   return pn
               elseif length(temp) > n+1 && temp[1] == t
                   pn = temp[n+2]

                   return pn
               else
                   temp = [t]
                   pn = (t - a[n])*polybeta(t,n-1,a,b,temp) - b[n-1]*polybeta(t,n-2,a,b,temp) #P_{n}(t) = (t-a_{n-1})P_{n-1} - b_{n-1}P_{n-2}
                   push!(temp, pn)

                   return pn
               end
          end
      end

function SDOhmic(t)
    """
    Bath Ohmic Spectral Density for zero temperature chain mapping of the bath
    """
    if t==0
        return 0
    elseif t>-1 && t<1
        return 2*α*abs(t)*ωc
    elseif abs(t)==1
        return 2*α*ωc
    else
        return 0
    end
end

    
function SDTOhmic(t)
    """
    Bath Ohmic Spectral Density after the finite temperature chain mapping of the bath
    """
    if t==0
        return 2*α/beta
    elseif t>-1 && t<1
        return α*t*ωc*(1+coth(beta*t*ωc*0.5))
    elseif abs(t)==1
        return α*t*ωc*(1+coth(beta*t*ωc*0.5))
    else
        return 0
    end
end

function switchmpo(N::Int, Nm::Int, dhilbert::Int; E=[], J=0.2, chainparams, s=1, α=0.02, ωc=1, Ra=20, R=1, c_phonon=1, beta ="inf", κ=1, signature = [], issoft=false)
    """
    This function construct a MPO for a system made of a TLS called switch (at x = 0) that can be excited or not, and two sites a (x = Ra) and b (x = Ra + R) coherently coupled in the single exciation subspace.
    The switch does not interact with sites a and b.  The switch and the two sites are interacting with a bosonic bath with Nm modes.
    In the chain representation of the bath, the interactions between the system and the bath are long range, i.e. each site interacts with several modes.
    """

    datadir = "./"
    if beta!="inf"
        datfname = "chaincoeffs_ohmic_a$(α)wc$(ωc)xc$(ωc/c_phonon)beta$(beta).csv"
        chaincoeffs = readdlm(string(datadir,datfname),',',Float64,'\n',skipstart=1)
    else 
        chaincoeffs = zeros(Nm,2)
    end

    a_chain = chaincoeffs[:,1]
    b_chain = chaincoeffs[:,2].^2

    Norm = zeros(Nm)

    # Definition of the coupling coefficient between the site x and the mode n for a Ohmic spectral density with a hard cut-off (Jacobi Polynomials) or a soft cut-off Laguerre Polynomials
    function γ(x::Int, n::Int, issoft::Bool; beta="inf", temp=[1.])
        if beta=="inf"
            if issoft==true
                polynomial0(t) = sf_laguerre_n(n,s,t)*exp(-im*t*(x-1)*R*ωc/c_phonon)*t^s*exp(-s)
                return sqrt(2*α*gamma(n+s + 1)/gamma(n+1))*ωc*quadgk(polynomial0, 0, 1)[1]
            else
                polynomial(t) = x==1 ? jacobi(2*t-1,n-1, 0, s)*t^s : jacobi(2*t-1,n-1, 0, s)*exp(-im*t*(Ra + (x-1)*R)*ωc/c_phonon)*t^s
		        return sqrt(2*α*(2*(n-1) + s + 1))*ωc*quadgk(polynomial, 0, 1)[1]
            end
        elseif beta!="inf"
           polynomial1(t) = polybeta(t,n-1,a_chain,b_chain,[t])
           integrand(t) = x==1 ? polynomial1(t)*SDTOhmic(t) : polynomial1(t)*exp(im*t*(Ra + (x-1)*R)*ωc/c_phonon)*SDTOhmic(t)
           N2(t) = polynomial1(t)^2*SDTOhmic(t)
           if Norm[n]==0
               Norm[n] = sqrt(quadgk(N2,-1,1)[1])
           end
           return (ωc/Norm[n])*quadgk(integrand, -1, 1)[1]

        end
    end


    # Construction of the MPO
    W = Any[] # list of the MPO's tensors
    d = 2 #N # Hilbert Space dimension of the sites operators
    u = unitmat(d)

    ### Construction of the system sites MPO ###

    if length(E)<N
        print("The on-site energies of the system parts of the interaction Hamiltonian where not all specified: zero energy will be assumed. \n")
        E = vcat(E,zeros(N-length(E)))
    end

    if length(signature)<N
        print("The signs of the system parts of the interaction Hamiltonian where not all specified: positive signs will be assumed. \n")
        signature = vcat(signature, ones(N-length(signature)))
    end

    # Definition of the single excitation creation, anihilation operators on site x and the projector on site x
    cd = crea(2)
    c = anih(2)
    P = [0 0;0 1]

    print("Sites MPO \n")
    if N==1
        error("The number of sites N is too small.")
    else
        for x = 1:N-1
            print("Site ",x,"\n")

            D = 2*(x+2) # Bond dimension
            M = zeros(ComplexF64,D-2,D,d,d)
            M[1,1,:,:] = M[D-2,D,:,:] = u

            i = 2 # index counter
            M[1,i,:,:] = x==1 ? zeros(ComplexF64, size(c)...) : J*c # The first site (the switch) is NOT coupled to the other sites
            M[1,i+1,:,:] = x==1 ? zeros(ComplexF64, size(c)...) : J*cd
            M[1,i+2*x,:,:] = M[1,i+1+2*x,:,:] = signature[x]*P # system operator coupled to the environment 
            M[1,D,:,:] = E[x]*P # onsite energy

            M[i,D,:,:] = cd
            M[i+1,D,:,:] = c
            i += 2

            while i<D-2
                M[i,i,:,:] = u
                i += 1
            end
            push!(W, M)
        end
        print("Last site before the chain \n")
        ### Last site before the bath chain doesn't have any coupling with the rest of the sites, so is size is Dx(D-2)xNxN

        D = 2*(N+1)
        M = zeros(ComplexF64,D, D, d, d)
        M[1,1,:,:] = M[D,D,:,:] = u
        i = 2 # index counter
        M[1,D-2,:,:] = M[1,D-1,:,:] = signature[N]*P 
        M[1,D,:,:] = E[N]*P

        M[i,D,:,:] = cd
        M[i+1,D,:,:] = c

        while i<D-2
            M[i+2,i,:,:] = u
            i += 1
        end
        push!(W, M)
    end

    print("Chain MPO \n")
    ### Construction of the bosonic bath MPO ###
    e = chainparams[1] # list of the energy of each modes
    t = chainparams[2] # list of the hopping energy between modes
    d = Int64(dhilbert) # Hilbert space dimension of the bosonic bath operators
    u = unitmat(d)
    # modes creation, anihilation and number operators
    bd = crea(d)
    b = anih(d)
    n = numb(d)

    coupling_stored = zeros(ComplexF64,N,Nm) # just need NxNm because the two chains have the same coupling coeff up to complex conjugation
    fnamecc = "chaincouplings_ohmic_N$(N)_Ra$(Ra)_R$(R)_a$(α)_wc$(ωc)_xc$(ωc/c_phonon)_beta$(beta).csv"
    arestored = 1 # are the coupling coefficient stored
    Nstored, Nmstored = N, Nm # number of stored coeff.
    couplinglist = []
    try
        couplinglist = readdlm(fnamecc,',',ComplexF64,'\n')
    catch stored_error
        if isa(stored_error, ArgumentError)
            print("Coupling coefficient not found. They will be computed and stored.\n")
            arestored = 0
        end
    end
    if couplinglist != []
        Nstored, Nmstored = size(couplinglist)
        if Nmstored>=Nm && Nstored>=N
            coupling_stored = couplinglist[1:N,1:Nm]
        else
            coupling_stored[1:min(N,Nstored),1:min(Nm,Nmstored)] = couplinglist[1:min(N,Nstored),1:min(Nm,Nmstored)]
            print("Less coupling coefficient stored than needed. Available ones will be used and missing one will be computed and stored.\n")
        end
    end

    if Nm != 1
        print("First Mode \n")
        ## First chain MPO
        D = 2*(N + 2) #Bond dimension
        M = zeros(ComplexF64,D-2, D, d, d)
        M[1,1,:,:] = M[D-2,D,:,:] = u
        M[1,D,:,:] = e[1]*n
        i = 2
        M[1,i,:,:] = t[1]*bd
        M[1,i+1,:,:] = t[1]*b

        a = 0 #site counter
        while i<D-2
            a+=1
            if arestored==1
                couplingcoeff = coupling_stored[a,1]
            else
                couplingcoeff = γ(a,1,issoft, beta=beta)
                coupling_stored[a,1] = couplingcoeff
            end
            couplingcoeff = a==1 ? sqrt(κ)*couplingcoeff : couplingcoeff
            M[i,D,:,:] = couplingcoeff*b
            M[i,i+2,:,:] = u
            i+=1
            M[i,D,:,:] = conj(couplingcoeff)*bd
            M[i,i+2,:,:] = u
            i+=1
        end
        M = reshape(M,D-2,D,d,d)
        push!(W, M)

        for m = 2:Nm-1
            print("Chain mode ",m,"\n")
            D = 2*(N + 2)
            M = zeros(ComplexF64,D, D, d, d)
            M[1,1,:,:] = M[D,D,:,:] = u
            M[1,D,:,:] = e[m]*n
            i = 2
            M[1,i,:,:] = t[m]*bd
            M[1,i+1,:,:] = t[m]*b
            M[i,D,:,:] = b
            M[i+1,D,:,:] = bd
            i += 2

            a = 0 #site counter
            while i<D
                a+=1
                if arestored==1 && m<=Nmstored && a<=Nstored
                    couplingcoeff = coupling_stored[a,m]
                else
                    couplingcoeff = γ(a, m, issoft,beta=beta)
                    coupling_stored[a,m] = couplingcoeff
                end
                M[i,D,:,:] = couplingcoeff*b
                if i<D-1
                    M[i,i,:,:] = u
                end
                i+=1
                M[i,D,:,:] = conj(couplingcoeff)*bd
                if i<D
                    M[i,i,:,:] = u
                end
                i+=1
            end

            M = reshape(M,D,D,d,d)
            push!(W, M)
        end

	print("Last Mode of the First Chain \n")
        D = 2*(N + 2)
        M = zeros(ComplexF64,D, D, d, d)
        M[1,1,:,:] = M[D,D,:,:] = u
        M[1,D,:,:] = e[Nm]*n
        i = 2
        M[i,D,:,:] = b
        M[i+1,D,:,:] = bd
        i += 2

        a = 0 #site counter
        while i<D
            a+=1
            if arestored==1 && Nm<=Nmstored && a<=Nstored
                couplingcoeff = coupling_stored[a,Nm]
            else
                couplingcoeff = γ(a, Nm, issoft,beta=beta)
                coupling_stored[a,Nm] = couplingcoeff
            end
            M[i,D,:,:] = couplingcoeff*b
            if i<D-1
                M[i,i,:,:] = u
            end
            i+=1
            M[i,D,:,:] = conj(couplingcoeff)*bd
            if i<D
                M[i,i,:,:] = u
            end
            i+=1
        end

            M = reshape(M,D,D,d,d)
            push!(W, M)

	print("Second chain \n")
	for m = 1:Nm-1
            print("Chain mode ",m,"\n")
            D = 2*(N + 2)
            M = zeros(ComplexF64,D, D, d, d)
            M[1,1,:,:] = M[D,D,:,:] = u
            M[1,D,:,:] = e[m]*n
            i = 2
            M[1,i,:,:] = t[m]*bd
            M[1,i+1,:,:] = t[m]*b
            M[i,D,:,:] = b
            M[i+1,D,:,:] = bd
            i += 2

            a = 0 #site counter
            while i<D
                a+=1
                couplingcoeff = a==1 ? sqrt(κ)*coupling_stored[a,m] : coupling_stored[a,m]
                M[i,D,:,:] = couplingcoeff*bd
                if i<D-1
                    M[i,i,:,:] = u
                end
                i+=1
                M[i,D,:,:] = conj(couplingcoeff)*b
                if i<D
                    M[i,i,:,:] = u
                end
                i+=1
            end

            M = reshape(M,D,D,d,d)
            push!(W, M)
        end

        print("Last mode \n")
        WNm = zeros(ComplexF64,D, 1, d, d)
        WNm[1,1,:,:] = e[Nm]*n
        WNm[2,1,:,:] = b
        WNm[3,1,:,:] = bd
        a = 0 #site counter
        i = 4 #row index
        while i<D
            a+=1
            couplingcoeff = coupling_stored[a, Nm]
            WNm[i,1,:,:] = couplingcoeff*bd
            i+=1
            WNm[i,1,:,:] = conj(couplingcoeff)*b
            i+=1
        end
        WNm[D,1,:,:] = u
        WNm = reshape(WNm,D,1,d,d)

    elseif Nm==1

        print("First and Last mode \n")
        D = 2*(N + 1) #Bond dimension
        WNm = zeros(ComplexF64,D, 1, d, d)
        WNm[1,1,:,:] = e[Nm]*n
        a = 0 #site counter
        i = 2 #row index
        while i<D
            a+=1
            couplingcoeff = γ(a,Nm, issoft, beta=beta)
            WNm[i,1,:,:] = couplingcoeff*b
            i+=1
            WNm[i,1,:,:] = conj(couplingcoeff)*bd
            i+=1
        end
        WNm[D,1,:,:] = u
        WNm = reshape(WNm,D,1,d,d)
    end

    W1 = W[1]

    if arestored==0 || Nmstored<Nm || Nstored<N
        writedlm(fnamecc, coupling_stored, ',')
    end

    return Any[W1[1:1,:,:,:], W[2:(N+2*Nm-1)]..., WNm]

end

function displacedchainmps(A::Vector{Any}, N::Int, Nm::Int; γ=nothing, chainparams=[fill(1.0,Nm),fill(1.0,Nm-1), 1.0], s=1, α=0.02, ωc=1, R=1, c_phonon=1, beta ="inf", issoft=false)
"""
For a displacement gamma of the bath modes, compute the corresponding displaced operator on the 2*Nm-long chain and apply it to a given mps A.
"""
    coupling_stored = zeros(ComplexF64,N,Nm) # just need NxNm because the two chains have the same coupling coeff up to complex conjugation
    fnamecc = "chaincouplings_ohmic_N$(N)_Ra$(Ra)_R$(R)_a$(α)_wc$(ωc)_xc$(ωc/c_phonon)_beta$(beta).csv"
    lcouplings = Nm # number of stored coeff.
    chaincouplings = readdlm(fnamecc,',',ComplexF64,'\n')

    v = chaincouplings[2,1:Nm]

    M = Tridiagonal(chainparams[2], chainparams[1], chainparams[2])

    ι = M\v
    ι = vcat(-1*conj(ι), -1*ι)

    B = Any[] # displaced chain MPS

    # The system part of the MPS should be unaffected by the displacement operator
    for i=1:N
        push!(B, A[i])
    end

    # Displacement operator for the chain
    for n=1:2*Nm
        d1, d2, dhilbert = size(A[N+n]) # left bond dim, right bond dim and system physical dim
        χ = d1*d2

        bd = crea(dhilbert) # chain creation operator
        b = anih(dhilbert) # chain anihilation operator

        W = ι[n]*bd - conj(ι[n])*b  # Argument of the displacement operator

        Ap = permutedims(reshape(A[N+n], χ, dhilbert), [2,1]) # reshape the MPS tensor as a (d, χ)-matrix

        As = permutedims(exp(W)*Ap, [2,1]) # matrix multiplication B_n = exp(W)*A_n

        push!(B, reshape(As, d1, d2, dhilbert))
    end

    return B
end
