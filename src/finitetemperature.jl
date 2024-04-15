using HDF5
include("../ChainOhmT/quadohmT.jl")
include("../ChainOhmT/mcdis2.jl")


"""
    chaincoeffs_finiteT(nummodes, β, ohmic=true; α, s, J, ωc=1, mc=4, mp=0, AB=nothing, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=true)

Generate chain coefficients ``[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0]`` for a harmonic bath at the inverse temperature β.

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
* AB: component intervals
* Mmax: maximum number of integration points
* save: if true the coefficients are saved
"""
function chaincoeffs_finiteT(nummodes, β, ohmic=true; α=1, s=1, J=nothing, ωc=1, mc=4, mp=0, AB=nothing, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=true)

    N = nummodes #Number of bath modes

    if AB==nothing 
        if mc==4
            AB = [[-Inf -ωc];[-ωc 0];[0 ωc];[ωc Inf]]
        else
            throw(ArgumentError("An interval AB with mc = $mc components should have been provided."))
        end
    elseif length(AB) != mc
        throw(ArgumentError("AB has a different number of intervals than mc = $mc."))             
    end

    if ohmic==true
        wf(x,i) = ohmicspectraldensity_finiteT(x,i,α,s,ωc,β)
    elseif J==nothing
        throw(ArgumentError("A spectral density should have been provided."))
    else
        wf = J
    end
    
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
    ab, Mcap, kount, suc, uv = mcdis(N,eps0,quadfinT,Mmax,idelta,mc,AB,wf,mp,irout)
    for m = 1:N-1
        jacerg[m,1] = ab[m,1] #site energy e
        jacerg[m,2] = sqrt(ab[m+1,2]) #hopping parameter t
    end
    jacerg[N,1] = ab[N,1]

    eta = 0.
    for i = 1:mc
        xw = quadfinT(Mcap,i,uv,mc,AB,wf)
        eta += sum(xw[:,2])
    end
    jacerg[N,2] = sqrt(eta) # system-chain coupling c

    if save==true
        # Write a HDF5 file
        curdir = @__DIR__

        if ohmic==true
            Nstr=string(N)
            astr=string(α)
            sstr=string(s)
            bstr=string(β)
            # the "path" to the data inside of the h5 file is beta -> alpha -> s -> data (e, t or c)

            # Write onsite energies
            h5write("$curdir/ohmicT/chaincoeffs.h5", string("/N=",Nstr,"/α=",astr,"/s=",sstr,"/β=",bstr,"/e"), jacerg[1:N,1])
            # Write hopping energies
            h5write("$curdir/ohmicT/chaincoeffs.h5", string("/N=",Nstr,"/α=",astr,"/s=",sstr,"/β=",bstr,"/t"), jacerg[1:N-1,2])
            # Write coupling coefficient
            h5write("$curdir/ohmicT/chaincoeffs.h5", string("/N=",Nstr,"/α=",astr,"/s=",sstr,"/β=",bstr,"/c"), jacerg[N,2])

        else
            Nstr = string(N)
            wstr = string(ωc)
            bstr = string(β)
            # the "path" to the data inside of the h5 file is N -> ωc -> beta -> data (e, t or c)

            # Write onsite energies
            h5write("$curdir/ohmicT/chaincoeffs.h5", string("/N=",Nstr,"/ωc=",wstr,"/β=",bstr,"/e"), jacerg[1:N,1])
            # Write hopping energies
            h5write("$curdir/ohmicT/chaincoeffs.h5", string("/N=",Nstr,"/ωc=",wstr,"/β=",bstr,"/t"), jacerg[1:N-1,2])
            # Write coupling coefficient
            h5write("$curdir/ohmicT/chaincoeffs.h5", string("/N=",Nstr,"/ωc=",wstr,"/β=",bstr,"/c"), jacerg[N,2])
        end
    end

    return [jacerg[:,1], jacerg[1:N-1,2],jacerg[N,2]]
end
