include("include/MPSDynamics.jl")

fdir="/media/angus/DATA/methyleneBlue/ohmic_a0100/"

dt=1
tmax=600
N=500

cpars = [
    ["d",[10,15,20]],
    ["Dmax",[5,10,15]]
]

chainDat = readdlm(joinpath(MATDIR, "chaincoeffs_ohmic_a0100wc0001xc0006b20.csv"), ',', Float64, '\n')

cp = [chainDat[1:N,1], chainDat[1:N-1,2], chainDat[end,2]]

H(d, Dmax) = methylblue_S1_mpo(0, N, d, cp)

s1 = unitcol(1, 2)
s0 = unitcol(2, 2)
ψ = (1/sqrt(2))*(s1 + s0)

A(d, Dmax) = productstatemps(physdimsmpo(H(d, Dmax)), Dmax, statelist=[ψ, fill(unitcol(1, d), N)...])

obs = [
    #        ["s0s0", (psi, args...)->real(measure1siteoperator(psi, s0*s0', 1))],
    #        ["s1s1", (psi, args...)->real(measure1siteoperator(psi, s1*s1', 1))],
]
cobs = [
    ["realdcf", (psi, args...)->real(measure1siteoperator(psi, s0*s1', 1))],
    ["imagdcf", (psi, args...)->imag(measure1siteoperator(psi, s0*s1', 1))]
]
pars = [
    ["beta",20],
    ["xc",0.06],
    ["wc",0.01],
    ["a",1.0],
    ["dt",dt],
    ["tmax",tmax],
    ["ψ",ψ]
]

B, convdat, dat, plts, deltas, unid =
    convtdvp(dt, tmax, A, H, fdir;
             params=pars,
             observables=obs,
             convobs=cobs,
             convparams=cpars,
             save=true,
             verbose=false,
             timed=false,
             lightcone=true
             )
