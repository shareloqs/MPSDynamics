spinSX(n) = OneSiteObservable("spinSX", sx, n)
spinSY(n) = OneSiteObservable("spinSY", sy, n)
spinSZ(n) = OneSiteObservable("spinSZ", sz, n)

SX = OneSiteObservable("SX", sx)
SY = OneSiteObservable("SY", sy)
SZ = OneSiteObservable("SZ", sz)

struct FockError <: Observable
    name
    cparams
    sysop
    sites
    eigensys
    FockError(cpars, sysop, sites) = new("FockError", cpars, sysop, sites, eigenchain(cpars))
end
FockError(sysop, N) = FockError([fill(0.0, N), fill(1.0, N-1), 0.5], sysop, (2, N+1))
Base.ndims(::FockError) = 0

"""
    measure(A, obs::FockError; t=0, kwargs...)
    
Return the measure of the observable obs on the MPS A. 

"""
function measure(A, obs::FockError; t=0, kwargs...)
    d = physdims(A)[obs.sites[1]:obs.sites[2]][1]
    all(x->x==d, physdims(A)[obs.sites[1]:obs.sites[2]]) || error("MPS has non-uniform local Hilbert space dimensions")
    x = disp(d)
    p = mome(d)
    h = obs.cparams[3] * obs.sysop
    U = (obs.eigensys.vectors)'
    S = obs.eigensys.values
    N = length(S)
    ca = transpose(U'*diagm(0=>exp.(im*t.*S))*U)
    c = reshape([real.(ca), -imag.(ca), imag.(ca), real.(ca)], 2, 2)
    γ = reshape(
        [
            ((transpose(transpose(c[1,1][1,:])*c[1,r] + transpose(c[1,2][1,:])*c[2,r]) *
            (transpose(c[1,1][1,:])*c[1,s] + transpose(c[1,2][1,:])*c[2,s]))) -
            (r==1 ? unitcol(1,N)*(transpose(c[1,1][1,:])*c[1,s] + transpose(c[1,2][1,:])*c[2,s]) : zero(c[1,1])) -
            (s==1 ? transpose(transpose(c[1,1][1,:])*c[1,r] + transpose(c[1,2][1,:])*c[2,r])*unitrow(1,N) : zero(c[1,1])) +
            (s==1 && r==1 ? unitcol(1,N) * unitrow(1,N) : zero(c[1,1]))
         for s=1:2 for r=1:2
        ], 2, 2)

    ρ = leftcontractmps(A, [h^2])

    u = unitmat(d+1)
    ud = diagm(0=>[fill(1,d)...,0])
    o = ud*anih(d+1)*(u-ud)*crea(d+1)*ud
    o = o[1:d,1:d]
        
    aad = measure(A, OneSiteObservable("", o, obs.sites), ρ=ρ)
    oxx = aad
    oxp = im .* aad
    opx = -im .* aad
    opp = aad
    ors = reshape([oxx, opx, oxp, opp], 2, 2)
    
    h2xx = measure(A, TwoSiteObservable("", x, x, obs.sites), ρ=ρ)
    h2pp = measure(A, TwoSiteObservable("", p, p, obs.sites), ρ=ρ)
    h2xp = measure(A, TwoSiteObservable("", x, p, obs.sites), ρ=ρ)
    h2px = h2xp'
    h2rs = reshape([h2xx, h2px, h2xp, h2pp], 2, 2)
    ϵ1 = sum([γ[r,s][l,q] * h2rs[r,s][l,q] for r=1:2 for s=1:2 for l=1:N for q=1:N])
    ϵ2 = sum([c[1,r][1,k] * c[1,s][1,k] * ors[r,s][k] for r=1:2 for s=1:2 for k=1:N])
    return [γ[r,s][l,q] * h2rs[r,s][l,q] for r=1:2 for s=1:2 for l=1:N for q=1:N]
    #    return ϵ1+ϵ2
#    return h2pp
end

function errorbar(op, e, t)
    Nt = length(t)
    eave = [(t[i+1]-t[i]) * sqrt((e[i]+e[i+1])/2) for i=1:Nt-1]
    int = [sum(eave[1:i]) for i=0:Nt-1]
    return sqrt(8) * opnorm(op) .* sqrt.(int) 
end
