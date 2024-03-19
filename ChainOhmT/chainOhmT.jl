module Chaincoeffs
using CSV
using HDF5
using Tables
include("quadohmT.jl")
include("mcdis2.jl")

"""
Code translated in Julia by Thibaut Lacroix in 10/2020 from a previous Matlab code
"""

export mc, mp, iq, idelta, irout, AB, a, wc, beta, DM, uv
## Spectral density parameters
a = 0.01 # Coupling parameter
wc = 0.035 # Cutoff frequency. Often set up as 1 for arbitrary units
s = 1 #Ohmicity parameter
beta = 100 # 1/(kB T)

xc = wc # Limit of definition segments

## discretisation parameters
mc=4 # the number of component intervals
mp=0 # the number of points in the discrete part of the measure (mp=0 if there is none)
iq=1 # a parameter to be set equal to 1, if the user provides his or her own quadrature routine, and different from 1 otherwise
idelta=2 # a parameter whose default value is 1, but is preferably set equal to 2, if iq=1 and the user provides Gauss-type quadrature routines
irout=2 # choice between the Stieltjes (irout = 1) and the Lanczos procedure (irout != 1)
AB =[[-Inf -xc];[-xc 0];[0 xc];[xc Inf]] # component intervals
N=30 #Number of bath modes
Mmax=5000
eps0=1e7*eps(Float64)

jacerg = zeros(N,2)

ab = 0.
ab, Mcap, kount, suc, uv = mcdis(N,eps0,quadfinT,Mmax)
for m = 1:N-1
    jacerg[m,1] = ab[m,1] #site energy
    jacerg[m,2] = sqrt(ab[m+1,2]) #hopping parameter
end
jacerg[N,1] = ab[N,1]

eta = 0.
for i = 1:mc
    xw = quadfinT(Mcap,i,uv)
    global eta += sum(xw[:,2])
end
jacerg[N,2] = sqrt(eta)

# Write a CSV file
curdir = @__DIR__

header = Vector([Symbol("α_n"), Symbol("β_n and η")])
CSV.write("$curdir/ohmicT/chaincoeffs_ohmic_a$(a)wc$(wc)xc$(xc)beta$(beta).csv", Tables.table(jacerg, header=header))

Nstr=string(N)
astr=string(a)
sstr=string(s)
bstr=string(beta)
# the "path" to the data inside of the h5 file is beta -> alpha -> s -> data (e, t or c)

# Write onsite energies
h5write("$curdir/ohmicT/chaincoeffs.h5", string("/",Nstr,"/",astr,"/",sstr,"/",bstr,"/e"), jacerg[1:N,1])
# Write hopping energies
h5write("$curdir/ohmicT/chaincoeffs.h5", string("/",Nstr,"/",astr,"/",sstr,"/",bstr,"/t"), jacerg[1:N-1,2])
# Write coupling coefficient
h5write("$curdir/ohmicT/chaincoeffs.h5", string("/",Nstr,"/",astr,"/",sstr,"/",bstr,"/c"), jacerg[N,2])

end

