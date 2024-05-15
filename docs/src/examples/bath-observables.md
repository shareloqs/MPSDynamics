# Inspecting the bath by undoing the chain mapping

Here we give some context on the example script provided in `MPSDynamics/example/bath-observables.jl`. This example demonstrates the setup and execution of a simulation for a two-level system coupled to an Ohmic bath at finite temperature, where we exploit the access to the chain observables to:
- undo the chain mapping[^chin_exact_2010], thus obtaining their representation in the extended bath of T-TEDOPA [^tamascelli_efficient_2019], characterized by $J(\omega, \beta)$
- inverting the thermofield transformation[^devega_thermo_2015], thus obtaining the representation of the physical frequencies in the original environment, characterized by $J(\omega)$

T-TEDOPA allows to substitute a thermally occupied bath by an extended one, in the pure state of the vacuum, extending the bath of frequencies to negative values: the _creation_ from the system of a mode of _negative frequency_ in the extended bath of frequencies corresponds to the _absorption_ of energy for the system from the thermally occupied modes of the environment. However, by exploiting the thermofield transformation, it is possible to recover the picture of the occupations of the steady states in the physical bath of frequencies[^riva_thermal_2023]. 

## The code

We start by defining the parameters of the simulation:

```julia
d = 10      # number of Fock states of the chain modes
N = 60      # length of the chain
α = 0.01    # coupling strength
ω0 = 0.2  # TLS gap
s = 1       # ohmicity
ωc = 1.  # Cut-off of the spectral density J(ω)
β = 20    # Thermalized environment
```
We set the specifics of the simulation:
```julia
method = :TDVP1         # time-evolution method
conv = 3                # bond dimension for the TDVP1
dt = 0.5                # time step
tfinal = 60.0           # simulation time
```
And then compute the chain coefficients, i.e. on-site energies $\epsilon_i$, hopping energies $t_i$, and system-chain coupling $c_0$, that define the chain representation of the pure-dephasing model:
```julia
cpars = chaincoeffs_finiteT(N, β; α=α, s=s, J=nothing, ωc=ωc, mc=4, mp=0, AB=nothing, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=false)
```
with this, we can compute the corresponding MPO, and the initial state in MPS form ($1/\sqrt{2}(|0\rangle + |1\rangle)$):
```julia
H = puredephasingmpo(ω0, d, N, cpars)

ψ = zeros(2)
ψ[1] = 1/sqrt(2)
ψ[2] = 1/sqrt(2)

A = productstatemps(physdims(H), state=[ψ, fill(unitcol(1,d), N)...])
```
We can now define the observables we are interested in computing. Importantly, we not only compute observables related to the system, but also to the chain modes, so that we can inspect the environment.
```julia
ob1 = OneSiteObservable("sz", sz, 1)
ob2 = OneSiteObservable("sx", sx, 1)
ob3 = OneSiteObservable("chain_mode_occupation", numb(d), (2,N+1))
ob4 = OneSiteObservable("c", crea(d), collect(2:N+1))
ob5 = OneSiteObservable("cdag", crea(d), collect(2:N+1))
ob6 = TwoSiteObservable("cdagc", crea(d), anih(d), collect(2:N+1), collect(2:N+1))
ob7 = TwoSiteObservable("cdagcdag", crea(d), crea(d), collect(2:N+1), collect(2:N+1))
ob8 = TwoSiteObservable("cc", anih(d), anih(d), collect(2:N+1), collect(2:N+1))
```
with this we run the simulation:
```julia
A, dat = runsim(dt, tfinal, A, H, prec=1E-4;
                name = "Bath observables in the pure dephasing model",
                method = method,
                obs = [ob1, ob2, ob3, ob4, ob5, ob6, ob7, ob8],
                convobs = [ob1],
                params = @LogParams(ω0, N, d, α, s, ψ),
                convparams = conv,
                reduceddensity = true,
                verbose = false,
                save = false,
                plot = true,
                );
```

After the simulation, we need to post-process the data a bit, in order to map the correlations from the chain representation to the original environment. We start by defining the matrices associated to the mean values of $c$ and $c^\dagger$:
```julia
T = length(dat["data/times"])

constr = Array{ComplexF64}(undef, N, N, T)
destr = Array{ComplexF64}(undef, N, N, T)
for t in 1:T
    constr[:,:,t] = diagm(0 => dat["data/cdag"][:,t])
    destr[:,:,t] = diagm(0 => dat["data/c"][:,t])
end
```

During a numerical simulation we work with chain of finite length $N$. This truncation on chain modes, introduces a sampling on the modes in the original star-like environment. To recover the frequency modes that are implicitly sampled, one has to diagonalize the tri-diagonal $N\times N$ matrix $H^\text{chain}$, where the diagonal is formed by the $e_n$ coefficients, that is the chain's frequencies, and the upper and lower diagonals by the $N-1$ hopping coefficients $t_n$. The unitary matrix that defines the change of basis from the star-like to the chain-like environment is $U_n$, constituted by the eigenvectors of $H^\text{chain}$. In the code, we use the `eigenchain` function:

```julia
omeg = eigenchain(cpars, nummodes=N).values
```

At each time step of the simulation, a number of one-site and two-sites observables where evaluated on the chain. To obtain their value in the extended bath of T-tedopa, characterized by $J(\omega,\beta)$, the unitary transformation that maps the extended bath Hamiltonian into the chain representation has to be reversed. For instance, when measuring the single site $\hat n^c_i=\hat c_i^\dagger \hat c_i$ occupation number, we are not measuring the occupation number of the bosonic mode associated to the $\omega_i$ frequency, but the occupation number of the $i-$th chain mode.
Therefore to calculate the number of modes of the environment associated to a specific frequency $\omega_i$, the mapping must be reversed, to obtain the diagonal representation of the bosonic number operator:

$$
\hat n^b_{i} = \hat b_i^\dagger \hat b_i = \sum_{k,l} U_{ik}^* \hat c_k^\dagger \hat c_l U_{li}.  
$$

This is done in the code using the `measuremodes(X, cpars[1], cpars[2])` function, which outputs the vector of the diagonal elements of the operators, in the following way:

```julia
bath_occup = mapslices(X -> measuremodes(X, cpars[1], cpars[2]), dat["data/cdagc"], dims = [1,2])
cdag_average = mapslices(X -> measuremodes(X, cpars[1], cpars[2]), constr, dims = [1,2])
c_average = mapslices(X -> measuremodes(X, cpars[1], cpars[2]), destr, dims = [1,2])
```

To compute the correlators, we need the full matrix in the original basis. We therefore use the `measurecorrs` function:

```julia
cc_average = mapslices(X -> measurecorrs(X, cpars[1], cpars[2]), dat["data/cc"], dims = [1,2])
cdagcdag_average = mapslices(X -> measurecorrs(X, cpars[1], cpars[2]), dat["data/cdagcdag"], dims = [1,2])

correlations_c = [
    cc_average[i, j, t] - c_average[i, 1, t] .* c_average[j, 1, t]
    for i in 1:size(cc_average, 1), j in 1:size(cc_average, 2), t in 1:size(cc_average, 3)
]
correlations_cdag = [
    cdagcdag_average[i, j, t] - cdag_average[i, 1, t] .* cdag_average[j, 1, t]
    for i in 1:size(cdagcdag_average, 1), j in 1:size(cdagcdag_average, 2), t in 1:size(cdagcdag_average,3)
]
```

It is possible to invert the thermofield transformation (details in [^riva_thermal_2023]). The expression of the mean value of the number operator for the physical modes can be expressed as a function of mean values in the extended bath, which we denote $\langle \hat a_{2k}^\dagger \hat a_{2k} \rangle$:

$$
    \langle \hat b_k^\dagger \hat b_k \rangle = \cosh{\theta_k}\sinh{\theta_k} (\langle \hat a_{2k}\hat a_{1k}\rangle + \langle \hat a_{1k}^\dagger\hat a_{2k}^\dagger\rangle ) + \sinh^2{\theta_k} (1+ \langle \hat a_{2k}^\dagger \hat a_{2k} \rangle ) ++ \cosh^2{\theta_k} \langle \hat a_{1k}^\dagger \hat a_{1k} \rangle
$$

We remark that in the thermofield case, a negative frequency $\omega_{2k}$ is associated to each positive frequency $\omega_{1k}$. The sampling is therefore symmetric around zero. This marks a difference with T-TEDOPA, where the sampling of frequencies was obtained through the thermalized measure $d\mu(\beta) = \sqrt{J(\omega, \beta)}d\omega$, and was not symmetric. To recover the results for the physical bath of frequencies starting from the results of our simulations, that were conducted using the T-TEDOPA chain mapping, we need to do an extrapolation for all of the mean values appearing in Eq. \ref{eq:physical_occupations}, in order to have their values for each $\omega$ at $-\omega$ as well. This is done in the code with the `physical_occup` function:

```julia
bath_occup_phys = physical_occup(correlations_cdag[:,:,T], correlations_c[:,:,T], omeg, bath_occup[:,:,T], β, N)
```

Finally, in the pure dephasing case, it is also possible to obtain the analytical prediction of the time evolution of the occupations of the bath's modes, so that we can compare our numerical results with the analytical ones, exploiting the Heisenberg time evolution relation:
$$
\frac{d \langle \hat b_\omega \rangle}{dt} = -i \langle[ \hat b_\omega, \hat H] \rangle = - i \omega \langle\hat b_\omega \rangle - i \frac{\langle \hat \sigma_x \rangle}{2} \sqrt{J(\omega, \beta)},     \\
\frac{d \langle \hat n_\omega \rangle}{dt} = -i \langle[\hat b_\omega^\dagger \hat b_\omega, \hat H] \rangle= 2 \frac{|J(\omega,\beta)|}{\omega} \sin(\omega t).
$$
To this end, it is convenient to choose one of the eigenstates of $\hat \sigma_z$ as the initial state, so that $\langle \hat \sigma_x \rangle = \pm 1$. By solving these differential equations, one obtains the time evolved theoretical behavior of the bath. We define the function for the comparison with analytical predictions:
```julia
Johmic(ω,s) = (2*α*ω^s)/(ωc^(s-1))

time_analytical = LinRange(0.0, tfinal, Int(tfinal))

Γohmic(t) = - quadgk(x -> Johmic(x,s)*(1 - cos(x*t))*coth(β*x/2)/x^2, 0, ωc)[1]

Decoherence_ohmic(t) = 0.5 * exp(Γohmic(t))


α_theo = 0.25 * α
function Jtherm(x)
    if 1 >= x >= 0
        return +α_theo * abs(x)^s * (1 + coth(β*0.5*x))
    elseif -1 <= x <= 0
        return -α_theo * abs(x)^s * (1 + coth(β*0.5*x))
    else
        return 0
    end
end

bath_occup_analytical(ω, t) = abs(Jtherm(ω))/(ω^2)*2*(1-cos(ω*t)) 
```
We conclude the example by plotting.
```julia
ρ12 = abs.(dat["data/Reduced ρ"][1,2,:])

p1 = plot(time_analytical, t->Decoherence_ohmic(t), label="Analytics", title=L"Pure Dephasing, Ohmic $s=%$s$, $\beta = %$β ~\mathrm{K}$", linecolor=:black, xlabel="Time (arb. units)", ylabel=L"Coherence $|\rho_{12}(t)|$", linewidth=4, titlefontsize=16, legend=:best, legendfontsize=16, xguidefontsize=16, yguidefontsize=16, tickfontsize=10)
p1 = plot!(dat["data/times"], ρ12, lw=4, ls=:dash, label="Numerics")

cumul = [bath_occup_analytical(omeg[i], tfinal)*(omeg[i+1]-omeg[i]) for i in 1:(length(omeg)-1)]

p2 = plot(omeg[1:length(omeg)-1], cumul, lw = 4, linecolor=:black,
             xlabel=L"\omega", ylabel=L"\langle n^b_\omega \rangle", label="Analytics",
             title="Mode occupation in the extended bath")
p2 = plot!(omeg, bath_occup[:, :, T], lw=4, ls=:dash, label="Numerics")

p3 = heatmap(omeg, omeg, abs.(real.(correlations_cdag[:,:,T]) .+ im*imag.(correlations_cdag[:,:,T])), 
            xlabel=L"\omega",
            ylabel=L"\omega", title="Environmental correlations")


Mhalf = Int(length(omeg)*0.5)+1
M = length(omeg)

p4 = plot(omeg[Mhalf:M], bath_occup_phys, lw=4,
            xlabel=L"\omega", ylabel=L"\langle n^b_\omega \rangle",
            title="Mode occupation in the physical bath")

plot(p1, p2, p3, p4, layout = (2, 2), size = (1400, 1200))
```

___________________
# References

[^chin_exact_2010]:
    > Chin, A. W.; Rivas, Á.; Huelga, S. F.; Plenio, M. B. Exact Mapping between System-Reservoir Quantum Models and Semi-Infinite Discrete Chains Using Orthogonal Polynomials. Journal of Mathematical Physics 2010, 51 (9), 092109. https://doi.org/10.1063/1.3490188.

[^tamascelli_efficient_2019]:
    > Tamascelli, D.; Smirne, A.; Lim, J.; Huelga, S. F.; Plenio, M. B. Efficient Simulation of Finite-Temperature Open Quantum Systems. Phys. Rev. Lett. 2019, 123 (9), 090402. https://doi.org/10.1103/PhysRevLett.123.090402.

[^devega_thermo_2015]:
    > de Vega, I.; Banuls, M-.C. Thermofield-based chain-mapping approach for open quantum systems. Phys. Rev. A 2015, 92 (5), 052116. https://doi.org/10.1103/PhysRevA.92.052116.

[^riva_thermal_2023]:
    > Riva, A.; Tamascelli, D.; Dunnett, A. J.; Chin, A. W. Thermal cycle and polaron formation in structured bosonic environments. Phys. Rev. B 2023, 108, 195138, https://doi.org/10.1103/PhysRevB.108.195138.
