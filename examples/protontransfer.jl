#=
    Example of a Proton Transfer model at zero temperature for an isolated system or with an environment characterised by an hard cut-off Ohmic spectral density J(ω) = 2αω when ω < ωc and 0 otherwise#

    The dynamics is simulated using the T-TEDOPA method that maps the normal modes environment into a non-uniform tight-binding chain.

    The RC tensor is initially displaced at the value γ corresponding to a space coordinate.

    The initial electronic system is initialized in the adiabatic LOW surface at the RC displacement. 

    H_S + H_RC + H_int^{S-RC} = ω^0_{e} |e⟩ ⟨e| + ω^0_{k} |k⟩ ⟨k| + \\Delta (|e⟩ ⟨k| + |k⟩ ⟨ e|) + ω_{RC} (d^{\dagger}d + \\frac{1}{2}) + g_{e} |e⟩ ⟨e|( d + d^{\dagger})+ g_{k} |k⟩ ⟨k|( d + d^{\dagger})
``
    H_B + H_int^{RC-B} = ∫_{-∞}^{+∞} dk ω_k b_k^\dagger b_k - (d + d^{\dagger})∫_0^∞ dω\\sqrt{J(ω)}(b_ω^\\dagger+b_ω) + λ_{reorg}(d + d^{\\dagger})^2
``.
    λ_{reorg} = ∫ \frac{J(ω)}{ω}dω

=#

using MPSDynamics, Plots, LaTeXStrings, ColorSchemes, PolyChaos, LinearAlgebra

#----------------------------
# Physical parameters
#----------------------------
# Enol / Keto

ω0e= 0.8  # Enol energy at x=0
ω0k= 0.8  # Keto energy at x=0

x0e = -0.25 # Enol Equilibrium displacement 
x0k = 0.25 # Keto Equilibrium displacement 

Δ = 0.05 # Direct coupling between enol and keto

m=1.83E3 # Mass of the reaction coordinate particle. 1.83E3 is the proton mass in atomic units. 
ħ=1 # Atomic Units convention

dFockRC = 25 # Fock space of the RC tensor

ωRC = 0.0347 # Frequency of the RC tensor

γ = sqrt(m*ωRC/2)*x0e # γ = displacement RC ; γ = \sqrt{m\omega/2} x_disp

cparsRC = [ωRC,ωRC*sqrt(m*ωRC/2)] # Energy and g RC coupling parameter

isolated = true # isolated = true : no environment ; isolated = false : system coupled to a bosonic environment
# Creates the chain depending on the isolated condition. The parameters for isolated = false can be modified as desired.
if isolated
    d=1; N=2; α = 0.0  # number of Fock states of the chain modes ; length of the chain ; coupling strength
else
    d=4; N=40; α = 0.008 # number of Fock states of the chain modes ; length of the chain ; coupling strength
end

s = 1 # ohmicity

ωc = 2*ωRC # Cut-off of the spectral density J(ω)

λreorg = (2*α*ωc)/s # Reorganisation Energy taken into account in the Hamiltonian

cpars = chaincoeffs_ohmic(N, α, s; ωc=ωc) # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0

#-----------------------
# Simulation parameters
#-----------------------

dt = 10.0 # time step

tfinal = 2000.0 # simulation time

numsteps = length(collect(0:dt:tfinal))-1

method = :TDVP1 # time-evolution method

D = 5 # MPS bond dimension


#---------------------------
# MPO and initial state MPS
#---------------------------

# Initial electronic system in the left well, at mean displacement . It results to a superposition of 1 and 2
# This state is obtained with the diagonalization of the diabatic Hamiltonian at x=γ
ψ = zeros(2)
X2 = zeros(2)

Ae = ω0e  - 0.5*m*ωRC^2*(x0e)^2
Ak = ω0k  - 0.5*m*ωRC^2*(x0k)^2
xini = γ/(sqrt(m*ωRC/2))

Adumb = (ω0e-ω0k+0.5*m*ωRC^2*((xini-x0e)^2-(xini-x0k)^2))
Bdumb = sqrt((Adumb)^2 + 4*Δ^2)
X2[1] = (Adumb-Bdumb)/(2*Δ)
X2[2] = 1

mod_X2 = sqrt(sum(X2[i]^2 for i=1:2))

ψ[1] = X2[1] / mod_X2
ψ[2] = X2[2] / mod_X2

#### Coherent state . Initialise the RC oscillator at displacement γ thanks to the Taylor expression of the displacement operator

coherent_RC = [fill(unitcol(1,dFockRC), 1)...]
chainpop_RCini =zeros(dFockRC,2)

chainpop_RCini[1,2]=1
chainpop_RCini[1,1]=exp(-(γ^2)/2)

for j=2:dFockRC            
    chainpop_RCini[j,1]=chainpop_RCini[j-1,1]*(γ/sqrt(j-1))
end

coherent_RC[1] = hcat(chainpop_RCini[:,1]...) #Manipulation to get good format for MPS, does not change value

H = protontransfermpo(ω0e, ω0k, x0e, x0k, Δ, dFockRC, d, N, cpars, cparsRC, λreorg)

A = productstatemps(physdims(H), state=[ψ, coherent_RC..., fill(unitcol(1,d), N)...]) # MPS representation of |ψ>|Displaced RC>|Vacuum>

Eini = measurempo(A,H)
print("\n Energy =",Eini)

#---------------------------
# Definition of observables
#---------------------------

ob1 = OneSiteObservable("sz", sz, 1)
#-------------
ob2 = OneSiteObservable("RC displacement", MPSDynamics.disp(dFockRC,ωRC,m), 2)

#-------------
# Simulation
#------------

A, dat = runsim(dt, tfinal, A, H;
                name = "proton transfer and tunneling at zero temperature",
                method = method,
                obs = [ob1, ob2],
                convobs = [ob1, ob2],
                params = @LogParams(N, d, α, s),
                convparams = D,
                Nrho = 2, # Need to precise that the reduced density matrix is
                # calculated for the 2nd tensor. The value by default is 1. 
                reduceddensity= true,
                verbose = false,
                save = true,
                plot = false,
                );

#-----------------------
# Plot parameters
#-----------------------

xlist_rho =collect(-1.2:0.05:1.2) # x definition for the reduced density matrix expressed in space
# other values of xlist_rho can be chosen to gain numerical time

palette = ColorSchemes.okabe_ito # color palette for plots

#-------------
# Plots
#------------

#### Plot the doublewell in the adiabatic basis. It is obtained with the diagonalization of the diabatic space Hamiltonian #####

function fp(ɸ)
    H11 = Ae + 0.5*m*ωRC^2*(ɸ-x0e)^2
    H22 = Ak + 0.5*m*ωRC^2*(ɸ-x0k)^2 
    H12 = Δ
    λp = ((H11+H22)/2)+0.5*sqrt((H11-H22)^2+4*H12^2)
    return λp
end

function fm(ɸ)
    H11 = Ae + 0.5*m*ωRC^2*(ɸ-x0e)^2
    H22 = Ak + 0.5*m*ωRC^2*(ɸ-x0k)^2
    H12 = Δ
    λm = ((H11+H22)/2)-0.5*sqrt((H11-H22)^2+4*H12^2)
    return λm
end

function gaussian(x)
    return Eini .+((1/0.06*sqrt(2*pi)).*exp.(-(x.+(0.25)).^2 ./(2*0.06^2))) .*0.001
end

ɸ = [i for i in -0.75:0.001:0.75]
ɸwp = [i for i in -0.45:0.001:-0.05]
xmin = -0.75 ; xmax = 0.75 ; ymin = 0.7 ; ymax = 0.9 

plt=plot(ɸ,fm, label="LOW",xlabel="Reaction Coordinate (arb. units)",ylabel=L"$\omega$ (arb. units)",linewidth=4,thickness_scaling = 1, bg_legend=:transparent, legendfont=16, legendfontsize=16, xguidefontsize=16, yguidefontsize=16, tickfontsize=(12), xlims=(xmin,xmax),ylims=(ymin,ymax),xticks=(xmin:0.25:xmax),yticks=(ymin:0.05:ymax),legend=(0.6,0.6), grid=false,linecolor =palette[8])

plt=plot!(ɸ,fp,label="UP",lc =palette[3], linewidth=3,thickness_scaling = 1)

plt = hline!([Eini], lw=2,  lc=:black, ls=:dash, label=L"\omega_{system}")

plt=plot!(ɸwp,gaussian,linecolor=:black,linewidth=3,  label="", fillrange = (Eini , gaussian(ɸwp)),c = :black)

display(plt)

##### Results #####

display(plot(dat["data/times"], dat["data/RC displacement"],label="<X> (arb. units)", linecolor =:black, xlabel="Time (arb. units)",ylabel="<X>", title="", linewidth=4, legend=:none, tickfontsize=10,xlims=(0,2000)))

#### Reduced Density Matrix in space#####

# Loads the Hermite coefficients to calculate the eigenstates of the harmonic oscillator. The coefficient 2 comes from the definition of the physics field. 
αherm, βherm = 2.0 .*rm_hermite(dFockRC)

#Fct building the nth harmonic eigenstate at coordinate x
# sqrt(m*ωRC/ħ)*sqrt(2)*x in the Hermite function to recover the form of the physical Hermite polynomial in the atomic units
function Ψ(n,x)
return Acoef(n)*(m*ωRC/(ħ*π))^(1/4)*exp(-m*ωRC*(x)^2/(2*ħ))*PolyChaos.evaluate(n,(sqrt(m*ωRC/ħ)*sqrt(2)*x),αherm,βherm)
end


#Function to compute the pre factor of Ψ(n,x)
  function Acoef(n)
    if n == 0
        return 1.0
    else
        return (1/sqrt(n*2))*2^(1/2)*Acoef(n-1)
    end
end

#### Calculation of functions to express the reduced densty matrix in space. These functions build the matrix of oscillator eigenstates
# Ψ(i,x) Ψ(j,x) ∀i,j ∈ dFockRC^2 and ∀ x ∈ xlist_rho
Pmatrix = zeros(length(xlist_rho),dFockRC)
Pmatrixnormed = zeros(length(xlist_rho),dFockRC)
normalization = zeros(dFockRC)

for n=1:dFockRC
    for x=1:length(xlist_rho)
        Pmatrix[x,n] = Ψ(n-1,xlist_rho[x])
    end
end

for n=1:dFockRC
    normalization[n] = sqrt(sum(Pmatrix[k,n]^2 for k=1:length(xlist_rho)))
end


for n=1:dFockRC
    Pmatrixnormed[:,n] = Pmatrix[:,n]/normalization[n]
end

### Converts the reduced density matrix expressed in Fock states in a reduced density matrix expressed in space.
ρx=zeros(ComplexF64,length(xlist_rho),length(xlist_rho),length(dat["data/times"]))
for t=1:length(dat["data/times"])
    print(" \n t = ", dat["data/times"][t])
    for x=1:length(xlist_rho)
        for y=1:length(xlist_rho)
            for n=1:dFockRC
                for m=1:dFockRC
                    ρx[x,y,t] +=  Pmatrixnormed[x,n]*Pmatrixnormed[y,m]*dat["data/Reduced ρ"][n,m,t]
                end
            end
        end
    end
end

#Creates a GIF of the diagonal elements of the reduced density over time, ie the population dynamics. Here, the reduced density matrix is expressed in space.

anim = @animate for t=1:length(dat["data/times"])
    k=(dat["data/times"][t])
    print("\n k = ", k)
    plot(xlist_rho, abs.(diag(ρx[:,:,t])), title="Time = $k (arb. units)", xlabel=L"Diag( $||\rho_{RC}(x,x)||$ )", ylabel="Amplitude", xlims=(-0.5, 0.5), ylims=(0,0.25),legend=:none)
end
display(gif(anim, "gif_population.gif", fps = 2.5))

#Creates a GIF of the reduced density expressed over time. Here, the reduced density matrix is expressed in space. Diagonal elements are populations 
#whereas anti-diagonal elements  represent coherences

anim = @animate for t=1:length(dat["data/times"])
    k=(dat["data/times"][t])
    print("\n k = ", k)
    (plot(heatmap(xlist_rho, xlist_rho, abs.(ρx[:,:,t]), c=:Blues ), title="Time = $k (arb. units)", xlabel=L"$x$ (arb. units)",ylabel=L"$x$ (arb. units)", tickfontsize=(12),colorbar_title = L"||\rho_{RC}(x,x)||", legend=:none, xlims=(-0.5, 0.5), ylims=(-0.5,0.5), clims=(0,0.18),aspect_ratio=:equal))
end
display(gif(anim, "gif_reducedrho.gif", fps = 2.5))
