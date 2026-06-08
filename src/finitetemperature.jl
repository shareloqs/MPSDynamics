using HDF5
include("../ChainOhmT/quadohmT.jl")
include("../ChainOhmT/mcdis2.jl")


"""
    chaincoeffs_finiteT(nummodes, β, ohmic=true; α, s, J, ωc=1, mc=4, mp=0, AB=nothing, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=true)

Generate chain coefficients ``[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]`` for a harmonic bath at the inverse temperature β. It can also save the generated coefficients in ../ChainOhmT/ohmicT/

By default a Ohmic spectral density ``J(ω) = 2αω_c (\\frac{ω}{ω_c})^s θ(ω-ω_c)`` is considered.
Users can provide their own spectral density.

# Arguments
* nummodes: Number of bath modes
* β: inverse temperature
* ohmic: true if the spectral density is Ohmic, false if the user provides its own spectral density
* α: Kondo parameter of the Ohmic spectral density 
* s: ohmicity
* J: user-provided spectral density. Should be a function f(x,i) where x is the frequency and i ∈ {1,...,mc} labels the intervals on which the SD is defined
* ωc: the maximum frequency of the Ohmic spectral density
* mc: the number of component intervals
* mp: the number of points in the discrete part of the measure (mp=0 if there is none)
* iq: a parameter to be set equal to 1, if the user provides his or her own quadrature routine, and different from 1 otherwise
* idelta: a parameter whose default value is 1, but is preferably set equal to 2, if iq=1 and the user provides Gauss-type quadrature routines
* procedure: choice between the Stieltjes and the Lanczos procedure
* AB: component intervals. The defaut intervals are [[-Inf -ωc];[-ωc 0];[0 ωc];[ωc Inf]].
* Mmax: maximum number of integration points
* save: if true the coefficients are saved
* smooth: if true the spectral density is multiplied by a decreasing exponential instead of a Heaviside function.
* jacobi: Jacobi weight function parameters for each interval of AB. Can help with integration of diverging SDs. If the (thermalized) SD behaves as ``\\sim ω^s`` at the left and as ``\\sim ω^q`` to the right of an interval i, the mc-by-2 matrix `jacobi` should have `jacobi[i,:] = [s, q]` (all other entries should be 0). 
"""
function chaincoeffs_finiteT(nummodes, β, ohmic=true; α=1, s=1, J=nothing, ωc=1, mc=4, mp=0, AB=nothing, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=true, smooth=false, jacobi=nothing)

    N = nummodes #Number of bath modes

    # Set the interval for the spectral density
    if AB==nothing 
        if mc==4
            AB = [[-Inf -ωc];[-ωc 0];[0 ωc];[ωc Inf]]
            smooth && @warn "It is advised to specify the AB intervals with smooth=true, current segments imply a hard cutoff at ωc=$ωc. \
            Example use: AB=[[-Inf -5*ωc];[-5*ωc 0];[0 5*ωc];[5*ωc Inf]]"
        else
            throw(ArgumentError("An interval AB with mc = $mc components should have been provided."))
        end
    elseif size(AB)[1] != mc
        throw(ArgumentError("AB has a different number of intervals than mc = $mc."))             
    end
 
    # Express the spectral density according to the intervals
    if ohmic==true
        jacobi = [[0 0];[0 s-1];[s-1 0];[0 0]]
        wf = (x,i) -> ohmicspectraldensity_finiteT(x,i,α,s,ωc,β; smooth=smooth, jacobi=true) 
    elseif isnothing(J)
        throw(ArgumentError("A spectral density should have been provided."))
    else
        jacobi = isnothing(jacobi) ? zeros(mc, 2) : jacobi
        @assert size(jacobi) == (mc, 2) "jacobi should be a mc=$mc by 2 matrix with the 2 Jacobi weight parameters for each interval of AB"
        @assert !any(any(jacobi .!= 0, dims=2) .&& any(isinf, AB, dims=2)) "non-zero jacobi can only be used on finite intervals"
        wf = (x,i) -> J(x,i)./jacweight(x,i,AB,jacobi) 
    end
    
    # Choose the procedure to calculate the chain coefficients
    if procedure==:Lanczos  # choice between the Stieltjes and the Lanczos procedure
        irout = 2 
    elseif procedure==:Stieltjes
        irout = 1 
    else
        throw(ArgumentError("Procedure should be either Lanczos or Stieltjes."))
    end
    
    eps0=1e7*eps(Float64)

    jacerg = zeros(N,2)

    ab = 0.
    ab, Mcap, kount, suc, uv = mcdis(N,eps0,quadfinT,Mmax,idelta,mc,AB,wf,mp,irout,jacobi) # Calculate the chain coefficients
    for m = 1:N-1
        jacerg[m,1] = ab[m,1] #site energy e
        jacerg[m,2] = sqrt(ab[m+1,2]) #hopping parameter t
    end
    jacerg[N,1] = ab[N,1]

    # Calculate the integral of the spectral density to get system-chain coupling c
    eta = 0.
    for i = 1:mc
        xw = quadfinT(Mcap,i,uv[i],mc,AB,wf)
        eta += sum(xw[:,2])
    end
    jacerg[N,2] = sqrt(eta) # system-chain coupling c

    if save==true
        # Write a HDF5 file
        #curdir = @__DIR__
        dir = @__DIR__
        curdir = abspath(joinpath(dir, "../ChainOhmT"))

        if ohmic==true
            Nstr=string(N)
            astr=string(α)
            sstr=string(s)
            bstr=string(β)
            # the "path" to the data inside of the h5 file is beta -> alpha -> s -> data (e, t or c)

            # Write onsite energies
            h5write("$curdir/ohmicT/chaincoeffs.h5", string("/", Nstr,"/",astr,"/",sstr,"/",bstr,"/e"), jacerg[1:N,1])
            # Write hopping energies
            h5write("$curdir/ohmicT/chaincoeffs.h5", string("/", Nstr,"/",astr,"/",sstr,"/",bstr,"/t"), jacerg[1:N-1,2])
            # Write coupling coefficient
            h5write("$curdir/ohmicT/chaincoeffs.h5", string("/", Nstr,"/",astr,"/",sstr,"/",bstr,"/c"), jacerg[N,2])


        else
            Nstr = string(N)
            wstr = string(ωc)
            bstr = string(β)
            # the "path" to the data inside of the h5 file is N -> ωc -> beta -> data (e, t or c)
          
            # Write onsite energies
            h5write("$curdir/ohmicT/chaincoeffs.h5", string("/", Nstr, "/", wstr, "/", bstr, "/e"), jacerg[1:N,1])
            # Write hopping energies
            h5write("$curdir/ohmicT/chaincoeffs.h5", string("/", Nstr, "/", wstr, "/", bstr,"/t"), jacerg[1:N-1,2])
            # Write coupling coefficient
            h5write("$curdir/ohmicT/chaincoeffs.h5", string("/", Nstr, "/", wstr, "/", bstr,"/c"), jacerg[N,2])
 
        end
    end

    return [jacerg[:,1], jacerg[1:N-1,2],jacerg[N,2]]
end

"""
    chaincoeffs_fermionic(nummodes, β, chain; ϵ=x, ωc=1, mc=4, mp=0, AB=nothing, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=true)

Generate chain coefficients ``[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]`` for a fermionic bath at the inverse temperature β. It can also save the generated coefficients in ../ChainOhmT/fermionicT/


# Arguments
* nummodes: Number of bath modes
* β: inverse temperature
* chain: 1 if the chain modes are empty, 2 if the chain modes are filled
* ϵ: user-provided dispersion relation. Should be a function f(x) where x is the wavenumber
* J: user-provided spectral density. Should be a function f(x) where x is the wavenumber
* ωc: the maximum frequency allowed in the spectral density
* mc: the number of component intervals
* mp: the number of points in the discrete part of the measure (mp=0 if there is none)
* iq: a parameter to be set equal to 1, if the user provides his or her own quadrature routine, and different from 1 otherwise
* idelta: a parameter whose default value is 1, but is preferably set equal to 2, if iq=1 and the user provides Gauss-type quadrature routines
* procedure: choice between the Stieltjes and the Lanczos procedure
* AB: component intervals. The defaut intervals are [[-Inf -ωc];[-ωc 0];[0 ωc];[ωc Inf]].
* Mmax: maximum number of integration points
* save: if true the coefficients are saved
"""
function chaincoeffs_fermionic(nummodes, β, chain; ϵ=nothing, J=nothing, ωc=1, mc=4, mp=0, AB=nothing, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=true)

    N = nummodes # Number of bath modes

    # Set the interval for the spectral density
    if AB==nothing 
        if mc==4
            AB = [[-Inf -ωc];[-ωc 0];[0 ωc];[ωc Inf]]
        else
            throw(ArgumentError("An interval AB with mc = $mc components should have been provided."))
        end
    elseif size(AB)[1] != mc
        throw(ArgumentError("AB has a different number of intervals than mc = $mc."))             
    end

    # Express the spectral density according to the intervals
    if ϵ==nothing
        throw(ArgumentError("A dispersion relation should have been provided."))
    elseif J==nothing
        throw(ArgumentError("The spectral density should be provided."))
    else
        wf(x,i) = fermionicspectraldensity_finiteT(x, i , β, chain, ϵ, J)
    end

    # Choose the procedure to calculate the chain coefficients
    if procedure==:Lanczos  # choice between the Stieltjes (irout = 1) and the Lanczos procedure (irout != 1)
        irout = 2 
    elseif procedure==:Stieltjes
        irout = 1 
    else
        throw(ArgumentError("Procedure should be either Lanczos or Stieltjes."))
    end
    
    eps0=1e7*eps(Float64)

    jacerg = zeros(N,2)

    ab = 0.
    ab, Mcap, kount, suc, uv = mcdis(N,eps0,quadfinT,Mmax,idelta,mc,AB,wf,mp,irout) # Calculate the chain coefficients
    for m = 1:N-1
        jacerg[m,1] = ab[m,1] #site energy e
        jacerg[m,2] = sqrt(ab[m+1,2]) #hopping parameter t
    end
    jacerg[N,1] = ab[N,1]

    # Calculate the integral of the spectral density to get system-chain coupling c
    eta = 0.
    for i = 1:mc
        xw = quadfinT(Mcap,i,uv,mc,AB,wf)
        eta += sum(xw[:,2])
    end
    jacerg[N,2] = sqrt(eta) # system-chain coupling c

    if save==true
        # Write a HDF5 file
        #curdir = @__DIR__
        dir = @__DIR__
        curdir = abspath(joinpath(dir, "../ChainOhmT"))

        Nstr=string(N)
        cstr=string(chain)
        bstr=string(β)
        # the "path" to the data inside of the h5 file is beta -> alpha -> s -> data (e, t or c)

        # Write onsite energies
        h5write("$curdir/fermionicT/chaincoeffs.h5", string("/", Nstr, "/", bstr, "/", cstr, "/e"), jacerg[1:N,1])
        # Write hopping energies 
        h5write("$curdir/fermionicT/chaincoeffs.h5", string("/", Nstr, "/", bstr, "/", cstr, "/t"), jacerg[1:N-1,2])
        # Write coupling coefficient
        h5write("$curdir/fermionicT/chaincoeffs.h5", string("/", Nstr, "/", bstr, "/", cstr, "/c"), jacerg[N,2])
    end

    return [jacerg[:,1], jacerg[1:N-1,2],jacerg[N,2]]
end


"""
    chaincoeffs_finiteT_discrete(β, ωdiscrete, Jωdiscrete; procedure=:Lanczos, Mmax=5000, save=true)

Generate chain coefficients ``[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]`` for a discrete harmonic bath at the inverse temperature β.

# Arguments
* β: inverse temperature. For 0 Kelvin, put Inf as β
* ωdiscrete: discrete frequency corresponding to the Jωdiscrete values
* Jωdiscrete: amplitude of the spectral density at the corresponding ωdiscrete
* N: number of chain sites for which coefficients are calculated. Default keeps the same number of modes as the input. When T !=0K, the exact mapping is for N=2*length(ωdiscrete)
* procedure: choice between the Stieltjes and the Lanczos procedure
* Mmax: maximum number of integration points
* save: if true the coefficients are saved

Warning : works only for equal spacing between frequencies.
"""
function chaincoeffs_finiteT_discrete(β, ωdiscrete, Jωdiscrete;N=length(ωdiscrete), procedure=:Lanczos, Mmax=5000, save=true)

    ω = ωdiscrete
    Jω = Jωdiscrete
    length(ω)== length(Jω) || throw(ErrorException("J(ω) has $(length(Jω)) values while there is $(length(ω)) frequecencies"))
    if β != Inf
       ω_pos = ω
       Jω_pos = Jω .* (coth.((β/2).*ω_pos[:]) .+ 1) ./2

       ω_neg = -ω
       Jω_neg = -Jω .* (coth.((β/2).*ω_neg[:]) .+ 1) ./2

       if ω_pos[1]==0.0 #discard freq=zero because leads to Nan
            Jω_neg[1] = Jω_pos[1]= 0.0
       end

       ω = vcat(reverse(ω_neg),ω_pos)
       Jω =  vcat(reverse(Jω_neg),Jω_pos)
    end

    mp=length(ω) # the number of points in the discrete part of the measure

    DM =Array{Float64}(undef,mp,2)
    for i=1:mp
       DM[i,1] = ω[i]
       DM[i,2] = Jω[i]
    end

    jacerg = zeros(N,2)

    ab = Array{Float64}(undef, N, 2)
    ab[:,2] = zeros(N,1) 

    if procedure==:Lanczos  # choice between the Stieltjes and the Lanczos procedure
        ab = lanczos(N,DM)
    elseif procedure==:Stieltjes
        ab = stieltjes(N,DM)
    else
        throw(ArgumentError("Procedure should be either Lanczos or Stieltjes."))
    end

    for m = 1:N-1
       jacerg[m,1] = ab[m,1] #site energy
       jacerg[m,2] = sqrt(ab[m+1,2]) #hopping parameter
    end
    jacerg[N,1] = ab[N,1]

    eta = sum((ωdiscrete[2]-ωdiscrete[1]) .* Jω) # quadrature, could be a source of numerical error, works for frequency equal spacing

    jacerg[N,2] = sqrt(eta) #coupling coeficient

    if save==true
        # Write a HDF5 file
        #curdir = @__DIR__
        dir = @__DIR__
        curdir = abspath(joinpath(dir, "../ChainOhmT"))

        Nstr = string(N)
        bstr = string(β)
        # the "path" to the data inside of the h5 file is N -> ωc -> beta -> data (e, t or c)
         
        # Write onsite energies
        h5write("$curdir/discreteT/chaincoeffs.h5", string("/", Nstr, "/", bstr, "/e"), jacerg[1:N,1])
        # Write opping energies
        h5write("$curdir/discreteT/chaincoeffs.h5", string("/", Nstr, "/", bstr,"/t"), jacerg[1:N-1,2])
        # Write coupling coefficient
        h5write("$curdir/discreteT/chaincoeffs.h5", string("/", Nstr, "/", bstr,"/c"), jacerg[N,2])
 
    end

    return [jacerg[:,1], jacerg[1:N-1,2],jacerg[N,2]]
end

