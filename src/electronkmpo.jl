mpo = AutoMPO()

N = 20
e = 1.0
t = 1.0

chainparams = [fill(e, N), fill(t, N-1)]

es = chainparams[1]
ts = chainparams[2]
hmat = diagm(0=>es, 1=>ts, -1=>ts)
F = eigen(hmat)
U = F.vectors

k = [1,3,10,15]

npart = length(k)

for n=1:N
    println("n=$n")
    for m=1:N
        for l=1:N
            for o=1:N
                global mpo
                mpo += (U[n,k[1]]*U[m,k[2]]*U[l,k[3]]*U[o,k[4]], "Cdagup", n, "Cdagup", m, "Cdagup", l, "Cdagup", o)
            end
        end
    end
end

sites = [Index(4,"Electron") for i=1:N]

H = MPOtoVector(MPO(mpo, sites))

A = productstatemps(N, 4)

applympo!(A, H)
