# The Anderson Impurity Model 

In these two examples, we use the fermionic chain mapping proposed in [^khon_efficient_2021] to perform tensor network simulations of the Single Impurity Anderson Model (SIAM). The SIAM Hamiltonian is defined as:

```math
    \hat H^\text{SIAM}  = \hat H_\text{loc} + \hat H_\text{hyb} + \hat H_\text{cond} = \overbrace{\epsilon_d \hat d^\dagger \hat d}^{\hat H_\text{loc}} + \underbrace{\sum_{k} V_k \Big( \hat d^\dagger \hat c_k + \hat c_k^\dagger \hat d \Big)}_{H_\text{hyb}} + \underbrace{\sum_k \epsilon_k \hat c_k^\dagger \hat c_k}_{H_I^\text{chain}}.
```
All of the operators obey to the usual fermionic anti-commutation relations: $\{\hat c_i, \hat c_j^\dagger \} = \delta_{ij}$, $\{\hat c_i, \hat c_j \} =\{\hat c_i^\dagger, \hat c_j^\dagger \} =0$ $\forall i,j$. The chain mapping is based on a thermofield-like transformation 
[^devega_thermo_2015], performed with fermions: ancillary fermionic operators $\hat c_{2k}$ are defined, one for each of the original fermionic modes $\hat c_{1k}$. A Bogoliubov transformation is then applied, so that two new fermionic modes $\hat f_{1k}$ and $\hat f_{2k}$ are defined as a linear combination of $\hat c_{1k}$ and $\hat c_{2k}$. Two chains are defined: the chain labelled $1$ for the empty modes, the chain labelled $2$ for the filled modes.
The following relations are used to define the functions equivalent to the spectral density of the bosonic case, one for each chain:
```math
\begin{split}
   &V_{1k} = V_{k} \sin \theta_k = \sqrt{\frac{1}{e^{\beta \epsilon_k}+1}} \\
   &V_{2k} = V_{k} \cos \theta_k = \sqrt{\frac{1}{e^{-\beta \epsilon_k}+1}}, 
\end{split}
```    
where we choose the spectral function that characterizes the fermionic bath to be: $V_k= \sqrt{1-k^2}$, and we define the dispersion relation as: $e_k = k$, that is, a linear dispersion relation with propagation speed equal to $1$. This latter choice corresponds to a model of metals (gapless energy spectrum). We select a filled state as the initial state of the defect.
Using the mapping proposed in [1], the chain Hamiltonian becomes:
```math
    \begin{split}
    \hat H^\text{chain}  = \hat H_\text{loc} &+ \sum_{i = \{1,2\}}\bigg[ J_{i,0} \Big(\hat d^\dagger \hat a_{i,0} + \hat d \hat a_{i,0}^\dagger \Big) + \\ &+ \sum_{n=1}^\infty  \Big( J_{i,n} \hat a_{i,n}^\dagger \hat a_{i,n-1} +  J_{i,n} \hat  a_{i,n-1}^\dagger \hat a_{i,n} \Big) + \sum_{n=0}^\infty E_{i,n} \hat  a_{i,n}^\dagger \hat a_{i,n} \bigg],
    \end{split}
```
where the $J_{i,n}$ coefficients are the couplings between the chain sites and the $E_{i,n}$ coefficients are the energies associated to each chain site. Clearly, the interactions are between nearest neighbors. This, combined with the fact that the fermions in our model are spinless, enables a straightforward mapping into fermionic operators of the bosonic creation and annihilation operators, that on their part obey to the bosonic commutation relations: $[\hat b_i, \hat b_j^\dagger] = \delta_{ij}$, $[\hat b_i, \hat b_j] =[\hat b_i^\dagger, \hat b_j^\dagger] =0$ $\forall i,j$. The mapping derived from Jordan-Wigner transformations for spinless fermions is:
```math
    \hat a_{i}^\dagger \hat a_{i+1} + \hat a_{i+1}^\dagger \hat a_{i} = \hat b_{i}^\dagger \hat b_{i+1} + \hat b_{i+1}^\dagger \hat b_{i}.  
```


## Double chain mapping


In the double chain geometry of Fig. \ref{subfig:double_ferm}, the MPOs bond dimension is: $\chi = 4$. The MPO has the structure defined in Eq. \ref{eq:mpo} and therefore can be seen as the product of the following matrices:
```math
    H = W_{1N} \cdot...\cdot W_{1 0} \cdot W_d \cdot W_{20} \cdot ... \cdot W_{2 N},  
```
where the matrices are defined as:
```math
\begin{split}
& W_{1N} = 
\begin{bmatrix}
\hat{\mathbb I} & J_{2,N} \hat b_{2,N}^\dagger & J_{2,N} \hat b_{2,N} & E_{2,N} \hat b_{2,N}^\dagger \hat b_{2,N} 
\end{bmatrix}, \quad   W_{1 0} = 
\begin{bmatrix}
 \hat{ \mathbb I} & J_{2,0} \hat b_{2,0}^\dagger & J_{2,0} \hat b_{2,0} & E_{2,0} \hat b_{2,0}^\dagger \hat b_{2,0}\\
0 &0 & 0 & \hat b_{2,0} \\
0 &0 & 0 & \hat b_{2,0}^\dagger \\
0 &0 & 0 & \hat{\mathbb I}
\end{bmatrix}, \\
& W_d = 
\begin{bmatrix}
 \hat{ \mathbb I} & \hat d^\dagger & \hat d & \epsilon_d \hat d^\dagger \hat d\\
0 &0 & 0 & \hat d \\
0 &0 & 0 & \hat d^\dagger \\
0 &0 & 0 & \hat{\mathbb I}
\end{bmatrix}\\ ,
& W_{2 0} =
\begin{bmatrix}
 \hat{ \mathbb I} & \hat b_{1,0}^\dagger & \hat b_{1,0} & E_{1,0} \hat b_{1,0}^\dagger \hat b_{1,0}\\
0 &0 & 0 & \hat J_{1,0}b_{1,0} \\
0 &0 & 0 & \hat J_{1,0}b_{1,0}^\dagger \\
0 &0 & 0 & \hat{\mathbb I}
\end{bmatrix}
, \quad W_{2 N} =
\begin{bmatrix}
 E_{2,N} \hat b_{2,N}^\dagger \hat b_{2,N} \\ J_{2,N} \hat b_{2,N} \\ J_{2,N} \hat b_{2,N}^\dagger \\ \hat{\mathbb I}
\end{bmatrix}.
\end{split}
```

The system starts from a filled state, the chain starts as in the Fermi sea.

## Interleaved chain mapping

The drawback of the double chain representation is that the particle-hole pairs are spatially separated in the MPS, creating correlations and therefore leading to a dramatic increase in the bond dimensions. This is why Kohn and Santoro propose an interleaved geometry, the advantages of which are thoroughly explained in [^khon_eff_2022]. Exploiting the interleaved representation, the interaction comes to be between next-nearest neighbors: a string operator appears in the Jordan-Wigner transformation from bosons to fermions:
```math
    \hat a_{i}^\dagger \hat a_{i+2} + \hat a_{i+2}^\dagger \hat a_{i} = \hat b_{i}^\dagger \hat F_{i+1} \hat b_{i+2} + \hat b_{i} \hat F_{i+1} \hat b_{i+2}^\dagger,
```
where the string operator $\hat F_i$ is defined as: $\hat F_i = (-1)^{\hat n_i} = \hat{\mathbb I} -2 \hat n_i = \hat{\mathbb I}-2 \hat b_i^\dagger \hat b_i$.  It is possible to find the analytical form also for MPOs with long range interaction \cite{mpo}. In the case of next-nearest neighbors interactions between spinless fermions, in the interleaved geometry of Fig. \ref{subfig:folded}, the MPO representation will require a bond dimension $\chi=6$. We explicitly write it as:
```math
    H = W_{d} \cdot W_{2 0} \cdot W_{1 0} \cdot...\cdot W_{2N} \cdot W_{1 N},  
```
where the matrices are defined as: 

```math
\begin{split}
& W_d = 
\begin{bmatrix}
 \hat{\mathbb I} & \hat d & \hat d^\dagger & 0 & 0 & E_{d} \hat d^\dagger \hat d 
\end{bmatrix}, \quad W_{2 0} = 
\begin{bmatrix}
 \hat{ \mathbb I} & \hat b_{2,0} & \hat b_{2,0}^\dagger & 0 & 0 & E_{2,0} \hat b_{2,0}^\dagger \hat b_{2,0}\\
0 &0 & 0 & \hat{F}_{2,0} & 0 & J_{2,0} \hat b_{2,0}^\dagger \\
0 &0 & 0 & 0 & \hat{F}_{2,0} & J_{2,0} \hat b_{2,0} \\
0 &0 & 0 & 0 & 0 &  0\\
0 &0 & 0 & 0 & 0 & 0 \\
0 &0 & 0 & 0 & 0 & \hat{\mathbb I}
\end{bmatrix}, \\
& W_{1 0} = 
\begin{bmatrix}
 \hat{ \mathbb I} & \hat b_{1,0} & \hat b_{1,0}^\dagger & 0 & 0 & E_{1,0} \hat b_{1,0}^\dagger \hat b_{1,0}\\
0 &0 & 0 & \hat{ F}_{1,0} & 0 & 0 \\
0 &0 & 0 & 0 & \hat{F}_{1,0} & 0 \\
0 &0 & 0 & 0 & 0 & J_{1,0} \hat b_{1,0}^\dagger \\
0 &0 & 0 & 0 & 0 & J_{1,0} \hat b_{1,0} \\
0 &0 & 0 & 0 & 0 & \hat{\mathbb I}
\end{bmatrix}, \\
& W_{2,N} = 
\begin{bmatrix}
 \hat{ \mathbb I} & \hat b_{2,N} & \hat b_{2,N}^\dagger & 0 & 0 & E_{2,N} \hat b_{2,N}^\dagger \hat b_{2,N}\\
0 &0 & 0 & \hat{F}_{2,N} & 0 & 0 \\
0 &0 & 0 & 0 & \hat{F}_{2,N} & 0 \\
0 &0 & 0 & 0 & 0 & J_{2,N} \hat b_{2,N}^\dagger \\
0 &0 & 0 & 0 & 0 & J_{2,N} \hat b_{2,N} \\
0 &0 & 0 & 0 & 0 & \hat{\mathbb I}
\end{bmatrix}, 
\quad W_{1 N}
\begin{bmatrix}
 E_{1,N} \hat b_{1,N}^\dagger \hat b_{1,N} \\ 0 \\0 \\ J_{1,N} \hat b_{1,N}^\dagger \\ J_{1,N} \hat b_{1,N} \\ \hat{\mathbb I}
\end{bmatrix} .
\end{split}
```
________________
### References

[^khon_efficient_2021]: 
    > Kohn, L.; Santoro, G. E. Efficient mapping for anderson impurity problems with matrix product states. Phys. Rev. B 2021, 104 (1), 014303. https://doi.org/10.1103/PhysRevB.104.014303.

[^devega_thermo_2015]:
    > de Vega, I.; Banuls, M-.C. Thermofield-based chain-mapping approach for open quantum systems. Phys. Rev. A 2015, 92 (5), 052116. https://doi.org/10.1103/PhysRevA.92.052116.

[^khon_eff_2022]:
    > Kohn L.; Santoro G. E. J. Stat. Mech. (2022) 063102 DOI 10.1088/1742-5468/ac729b.
    