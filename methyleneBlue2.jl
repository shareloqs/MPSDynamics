using MPSDynamics
MPSDynamics.MATDIR="/home/angus/Documents/MATLAB/OPO/"
MPSDynamics.DEFSAVEDIR="/home/angus/Documents/Julia/test/"

fdir="/media/angus/DATA/methyleneBlue/B3LYP_QM/"

MACHINES = [
    Machine("alex-pc","/home/angus/bin/julia", 1)
]

dt=20.
tmax=20000

e1 = 0.08449053845529315
e2 = 0.09095192787848483
δ = -8.442963074920613e-05

dm1 = 3.63548038842 
dm2 = 0.201000934417

T = 300
chainDat1 = readdlm(joinpath(MATDIR, "chaincoeffs_B3LYP_QM_S1S2_T300.csv"), ',', Float64, '\n')
chainDat2 = readdlm(joinpath(MATDIR, "chaincoeffs_B3LYP_QM_S1_T300.csv"), ',', Float64, '\n')
chainDat3 = readdlm(joinpath(MATDIR, "chaincoeffs_B3LYP_QM_S2_T300.csv"), ',', Float64, '\n')

N = 250

cp_coupling = [chainDat1[1:N,1], chainDat1[1:N-1,2], chainDat1[end,2]]
cp_s1 = [chainDat2[1:N,1], chainDat2[1:N-1,2], chainDat2[end,2]]
cp_s2 = [chainDat3[1:N,1], chainDat3[1:N-1,2], chainDat3[end,2]]

d = 20

H = methylbluempo(e1, e2, δ, N, N, N, d, d, d, cp_s1, cp_s2, cp_coupling);

s2 = unitcol(1, 3)
s1 = unitcol(2, 3)
s0 = unitcol(3, 3)
a = sqrt(1/(1+dm1^2+dm2^2))
ψ = a*(s0 + dm1*s1 + dm2*s2)

μ = (1/a^2)*(dm1*s0*s1' + dm2*s0*s2')

A = productstatemps(H.tree, physdims(H); state=[ψ, fill(unitcol(1, d), 3*N)...]);

ob1 = TwoSiteObservable("adaga_s1", crea(d), anih(d), (2, N+1))
ob2 = TwoSiteObservable("adaga_s2", crea(d), anih(d), (N+2, 2N+1))
ob3 = TwoSiteObservable("adaga_s1s2", crea(d), anih(d), (2N+2, 3N+1))

ob4 = OneSiteObservable("occ_s1", numb(d), (2, N+1))
ob5 = OneSiteObservable("occ_s2", numb(d), (N+2, 2N+1))
ob6 = OneSiteObservable("occ_s1s2", numb(d), (2N+2, 3N+1))

ob7 = OneSiteObservable("dcf", μ, 1)
ob8 = OneSiteObservable("s0", s0*s0', 1)
ob9 = OneSiteObservable("s1", s1*s1', 1)
ob10 = OneSiteObservable("s2", s2*s2', 1)
ob11 = OneSiteObservable("s1s2", s1*s2', 1)
ob12 = OneSiteObservable("s0s2", s0*s2', 1)
ob13 = OneSiteObservable("s0s1", s0*s1', 1)


out = runtdvp_fixed!(dt, tmax, A, H;
                     obs = [ob1,ob2,ob3,ob4,ob5,ob6,ob7],
                     convobs = [ob7,ob8,ob9,ob10,ob11,ob12,ob13],
                     Dmax=[15,25],
                     verbose=false,
                     savedir=fdir,
                     lightcone=true,
                     save=true,
                     saveplot=true,
                     params = @LogParams(ψ, e1, e2, δ, dm1, dm2, T, N, d)
                     )

A, dat, convdat = out
