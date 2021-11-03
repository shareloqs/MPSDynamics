module Chaincoeffs
using CSV
using Tables
include("quadohmT.jl")
include("mcdis2.jl")

"""
Code translated in Julia by Thibaut Lacroix in 10/2020 from a previous Matlab code
"""

export mc, mp, iq, idelta, irout, AB, a, wc, beta, DM, uv
## Spectral density parameters
a = 0.03
wc = 1
beta = 1
xc = wc

## discretisation parameters
mc=4 # the number of component intervals
mp=0 # the number of points in the discrete part of the measure (mp=0 if there is none)
iq=1 # a parameter to be set equal to 1, if the user provides his or her own quadrature routine, and different from 1 otherwise
idelta=2 # a parameter whose default value is 1, but is preferably set equal to 2, if iq=1 and the user provides Gauss-type quadrature routines
irout=2 # choice between the Stieltjes (irout = 1) and the Lanczos procedure (irout != 1)
AB =[[-Inf -xc];[-xc 0];[0 xc];[xc Inf]] # component intervals
N=1000 #Number of bath modes
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
jacerg[N,2] = sqrt(eta/pi)
header = Vector([Symbol("α_n"), Symbol("β_n and η")])
CSV.write("chaincoeffs_ohmic_a$(a)wc$(wc)xc$(xc)beta$(beta).csv", Tables.table(jacerg, header=header))

end
