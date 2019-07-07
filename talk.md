name: inverse
layout: true
class: center, middle, inverse
---
## Molecular Properties in Solution 

### The Good, the bad and the ugly

.author[Roberto Di Remigio]

.institution[UiT - The Arctic University of Norway]

.institution[Department of Chemistry - Virginia Tech]

.date[8 July 2019, Lille]

.footnote[[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) licensed.
Browse slides at [http://tinyurl.com/talk-lille](http://tinyurl.com/talk-lille)]

???

A whirlwind tour of solvation models for the calculation of molecular properties
in solution, showing how variational formulations can be hugely advantageous.
I'll conclude by showing recent advances in stochastic methods for CC theory.

---
layout: false

## The Problem of Solvation

Aspects computational models need to capture
- Dynamical features of solvation
- Average nature of observables
- Variety of interactions

<p style="text-align:center;"><img src="images/solvation-dynamics.gif" style="width: 70%"></p>
<p style="clear: both;">

???

- Chemistry is a wet science: experiments mostly happen in solution.
- We can classify solvent effects as:
  * Direct: these effects stem straightforwardly from the modification
    underwent by the solute electronic density when interacting with the
    environment.
  * Indirect: it is common for solutes to exhibit different minimum-energy
    conformations in different environments. These effects are commonly
    labelled as indirect.
  * Local field: light-matter interactions are also affected by the
    environment. Local modifications of externally applied fields subtly
    influence molecular responses.
  * Dynamic: the presence of the environment radically influences excited
    states, since relaxation processes in the medium become important.
  * Specific:  This catch-all category includes all effects stemming from the
    peculiar solute-solvent pair interactions that cannot be fully described
    under any of the previous labels. In general, modelling such effects
    demands an atomistic level of detail.

<font color="red">
Problems
</font>
<p style="clear: both;">

- Size \\(\Leftrightarrow\\) accurate _ab initio_ not possible
- Conformations \\(\Leftrightarrow\\) statistical sampling

- Size is a huge problem. It is simply not possible to model large systems very
  accurately.
- Moreover, the conformational space is very large and needs to be sampled
  extensively to have meaningful results.

---
class: split-50-50

## Quantum/Classical Multiscale Models.red[<sup>1</sup>]

- Use **quantum** and **classical** physics together
  * quantum for chemically relevant subsystem
  * classical for the environment

.column[
<p style="text-align:center;"><img src="images/gfp_barrel.png" style="width: 100%"></p>
]
.column[
<p style="text-align:center;"><img src="images/pyridine+water_12AA_QMMM.png" style="width: 50%"></p>
]
<p style="clear: both;">
.footnote-cite[.red[<sup>1</sup>] Senn, H. M.; Thiel, W. _Angew. Chem. Int. Ed._ (2009), __48__, 1198]

???

- The idea of multiscale models is to _focus_ on the chemically relevant part
  of the system and treat the environment approximately.

---
class: split-60-40

## QM/Continuum.red[<sup>2</sup>]

- Use **quantum** and **classical** physics together
  * quantum for chemically relevant subsystem
  * _continuum dielectric_ for the environment

.column[
.red[Pros]
- **Exact** electrostatics
- Self-consistent polarization
- Statistically averaged _by construction_

.red[Cons]
- **No** chemical detail in the environment
- Dispersion and repulsion approximate
]
.column[
<p style="text-align:right;"><img src="images/pyridine_Continuum.png" style="width: 30%"></p>
]
<p style="clear: both;">
.footnote-cite[.red[<sup>2</sup>] Tomasi, J.; Mennucci, B.; Cammi, R. _Chem. Rev._ (2005), __105__, 2999]

???

---
class: split-60-40

## QM/Continuum: The Polarizable Continuum Model.red[<sup>2</sup>]

.column[
### Transmission problem
`$$
  \begin{align}
  L_\mathrm{i} u &= \nabla^2 u = -4\pi\rho \,\, \text{in}\,\, \Omega_\mathrm{i} \label{eq:internal} \\
  L_\mathrm{e} u &= 0 \,\, \text{in}\,\, \Omega_\mathrm{e} \label{eq:external} \\
  [u](\mathbf{s}) &= u_\mathrm{e} - u_\mathrm{i} = 0 \,\, \text{on}\,\, \Gamma
  \label{eq:trace-jump} \\
[\partial_L u](\mathbf{s}) &= \partial_{L_\mathrm{e}} u - \partial_{L_\mathrm{i}} u = 0 \,\, \text{on}\,\, \Gamma \label{eq:conormal-jump} \\
|u(\mathbf{r})| &\leq C \|\mathbf{r} \|^{-1} \,\,\text{for}\,\,\| \mathbf{r} \|\rightarrow\infty
\label{eq:radiation}
\end{align}
$$`
]
.column[
<p style="text-align:right;"><img src="images/alanine.svg" style="width: 30%"></p>
]

<p style="clear: both;">
.footnote-cite[.red[<sup>2</sup>] Tomasi, J.; Mennucci, B.; Cammi, R. _Chem. Rev._ (2005), __105__, 2999]

???

- Replace environment with continuum \\(\varepsilon\\)
- Create cavity in continuum \\(\Omega_i\\)
- Vacuum inside cavity \\(\varepsilon=1\\)
- Solute charge density _entirely_ in \\(\Omega_i\\)

- Model the solvent as a polarizable dielectric continuum
- Parameters for the definition of the boundary, i.e. the cavity
- Parameters describing the solvent: permittivity (static and optical)
- Notice that the solvent parameters are, _by definition_, averaged!
- \\(L_\star\\) are elliptic differential operators
- Trace operators are the extension of the concept of restriction of a function
  over a boundary to generalized functions in Sobolev space
- Conormal derivative extends the notion of a normal derivative to functions in
  Sobolev spaces
- \\(\sigma(\mathbf{s})\\) is called the apparent surface charge (ASC)

* \\(L_\star\\) are elliptic differential operators
* Dirichlet condition: \\([u] (\mathbf{s})\\)
* Neumann condition: \\([\partial_L u] (\mathbf{s})\\)

---
## Mathematics of PCM.red[<sup>3</sup>]

- Define the .red[reaction potential]
`$$
 u(\mathbf{r}) = \color{Blue}{\varphi(\mathbf{r})} + \color{Red}{\xi(\mathbf{r})}
 = \int_C \mathop{}\!\mathrm{d}\mathbf{r}^\prime \frac{\color{Blue}{\rho(\mathbf{r}^\prime)}}{|\mathbf{r} - \mathbf{r}^\prime|} +
  \int_{\partial C} \mathop{}\!\mathrm{d}\mathbf{s} \frac{\color{Red}{\sigma(\mathbf{s})}}{|\mathbf{r} - \mathbf{s}|}
$$`

- .red[Apparent surface charge (ASC)]
`$$
\color{Green}{\mathcal{T}}\color{Red}{\sigma(\mathbf{s})} = -\color{Green}{\mathcal{R}}\color{Blue}{\varphi(\mathbf{s})}
$$`

- Green's functions for \\(L_\star\\) define integral operators

`$$
 \begin{align}
  (\color{Green}{\mathcal{S}_\star} f)(\mathbf{s}) &=
\int_{\partial C} \mathop{}\!\mathrm{d}\mathbf{s}^\prime \color{Green}{G_\star(\mathbf{s}, \mathbf{s}^\prime)}f(\mathbf{s}^\prime) \\
  (\color{Green}{\mathcal{D}_\star} f)(\mathbf{s}) &=
\int_{\partial C} \mathop{}\!\mathrm{d}\mathbf{s}^\prime [\partial_{L_\star}^\prime\color{Green}{G_\star(\mathbf{s}, \mathbf{s}^\prime)}]f(\mathbf{s}^\prime) \\
  (\color{Green}{\mathcal{D}^\dagger_\star} f)(\mathbf{s}) &=
\int_{\partial C} \mathop{}\!\mathrm{d}\mathbf{s}^\prime [\partial_{L_\star}\color{Green}{G_\star(\mathbf{s}, \mathbf{s}^\prime)}]f(\mathbf{s}^\prime)
 \end{align}
$$`

.footnote-cite[.red[<sup>3</sup>] Cancès, E.; Mennucci, B. _J. Math. Chem._ (1998), __23__, 309]

???

- We've transformed a boundary value problem (BVP) into a boundary integral equation (BIE)
- Integral operators are defined in terms of traces and conormal derivatives
- The integral operators have well-defined mapping properties between Sobolev
  spaces of fractional order
- Knowledge of the Green's functions inside and outside the cavity is key to
  the method

---

## A PCM for All Seasons.red[<sup>3</sup>]

- Environments with anisotropies and nonuniformities
`$$
\begin{align}
\left[\left(2\pi\mathcal{I} - \color{Green}{\mathcal{D}_\mathrm{e}}\right)\color{Green}{\mathcal{S}_\mathrm{i}} +
      \color{Green}{\mathcal{S}_\mathrm{e}}\left(2\pi\mathcal{I} + \color{Green}{\mathcal{D}^\dagger_\mathrm{i}}\right) \right]&\color{Red}{\sigma(\mathbf{s})} = \\
&-\left[\left(2\pi\mathcal{I} - \color{Green}{\mathcal{D}_\mathrm{e}}\right) -
\color{Green}{\mathcal{S}_\mathrm{e}}\color{Green}{\mathcal{S}_\mathcal{i}^{-1}}\left(2\pi\mathcal{I} - \color{Green}{\mathcal{D}_\mathrm{i}}\right)\right]\color{Blue}{\varphi(\mathbf{s})}
\end{align}
$$`
- Isotropic environments
`$$
\left[ 2\pi \left(\frac{\varepsilon+1}{\varepsilon-1}\right)\mathcal{I} - \color{Green}{\mathcal{D}}\right]
\color{Green}{\mathcal{S}}\color{Red}{\sigma(\mathbf{s})} = -\left( 2\pi\mathcal{I} - \color{Green}{\mathcal{D}} \right)
\color{Blue}{\varphi(\mathbf{s})}
$$`
- Conductor-like boundary conditions
`$$
\color{Green}{\mathcal{S}}\color{Red}{\sigma(\mathbf{s})} = -f(\varepsilon)\color{Blue}{\varphi(\mathbf{s})}$$`

<font color="red">
<center>How to solve these?</center>
</font>
<p style="clear: both;">

.footnote-cite[.red[<sup>3</sup>] Cancès, E.; Mennucci, B. _J. Math. Chem._ (1998), __23__, 309]

???

---
class: split-70-30

## Boundary Element Method and PCM.red[<sup>4</sup>]

.column[
Solution by a boundary element method (BEM)

* Cavity of interlocking, atom-centered spheres
* \\(N_\mathrm{ts}\\) finite elements on the cavity boundary
* Form boundary integral operators
`$$
 \color{Green}{\mathbf{T}}\color{Red}{\mathbf{q}} = - \color{Green}{\mathbf{R}}\color{Blue}{\mathbf{v}}
$$`
* Solve linear system
]
.column[
<p style="text-align:right;"><img src="images/benzene_GePol.png" style="width: 30%"></p>
]
<p style="clear: both;">

- Similar equation for IEF-PCM, isotropic PCM and COSMO
- _Independent_ of QM method!

.footnote-cite[.red[<sup>4</sup>] Ern, A; Guermond, J.-L. _Theory and Practice of Finite Elements_, Springer, 2004]

???

- _Galerkin_ or _collocation_ method
- Direct inversion or iterative solver

---
## More Mathematics of PCM.red[<sup>5</sup>]

`$$
\left\{
\begin{aligned}
  &\text{Seek $\zeta \in V$ such that:}\\
  &\forall \zeta \in V \,\,\,\,\,
  \color{Red}{a(\zeta, \sigma)} = \color{Blue}{b(\zeta)}
\end{aligned}
\right.
$$`

**Lax-Milgram Lemma**

If the bilinear form \\(a\\) is _continous_ and _coercive_
in \\(V\\), then, for any linear form \\(b\\), the Problem is
_well-posed_.

**Variational Corollary**

If the bilinear form is _symmetric_ and _positive_
the unique solution to the Problem is
the _unique minimum_ on \\(V\\) of the functional:
`$$
  \mathcal{F}(\zeta) = \frac{1}{2}\color{Red}{a(\zeta, \zeta)} - \color{Blue}{b(\zeta)}
$$`

.footnote-cite[.red[<sup>5</sup>] Lipparini, F.; Scalmani, G.; Mennucci, B.; Cancès, E.; Caricato, M.; Frisch, M. J. _J. Chem. Phys._ (2010), __133__, 014106]

???

- Introduce a suitable vector space, equipped with a norm
  derived from a scalar product.
  The solution belongs to such a vector space, which is
  bigger that in the strong formulation (less constraints on
  the regularity of the solution)
- Hadamard well-posedness definition:
  The Problemis well-posed if it admits one and only
  one solution and the solution is bounded by the initial data, by
  an a priori estimate (see Ern&Guermond, p.82, Def. 2.1)
- Def. Continuity: A bilinear form on a normed vector space is
  bounded, or continuous, if there is a constant \\(C\\) such that for
  all \\(u, v \in V\\):
  `$$
    a(u, v) \leq C \lVert u\rVert \lVert v \rVert
  $$`
- Def. Coercivity: A bilinear form on a vector space is coercive,
  or elliptic, if there is a constant \\(c > 0\\) such that \\(\forall u
  \in V\\):
  `$$
    a(u, u) \geq c \lVert u\rVert^2
  $$`
  This implies that no eigenvalue of the linear operator associated
  to the bilinear form can be zero, hence that the operator is
  invertible!!

---

## Variational Formulation of the PCM.red[<sup>5</sup>]

Find \\(\color{Red}{\sigma(\mathbf{s})}\\) minimizing
`$$
U_\mathrm{PCM}[\color{Red}{\sigma(\mathbf{s})}]
=\frac{1}{2}\int_{\partial C} \mathop{}\!\mathrm{d}\mathbf{s} \color{Red}{\sigma(\mathbf{s})} \color{Green}{[\mathcal{R}^{-1}\mathcal{T}]} \color{Red}{\sigma(\mathbf{s})}
+\int_{\partial C} \mathop{}\!\mathrm{d}\mathbf{s} \color{Blue}{\varphi(\mathbf{s})}\color{Red}{\sigma(\mathbf{s})}
$$`

1. Physically meaningful
   - Unfavourable self-interaction of \\(\color{Red}{\sigma(\mathbf{s})}\\)
   - Favourable interaction of \\(\color{Red}{\sigma(\mathbf{s})}\\) and  \\(\color{Blue}{\varphi(\mathbf{s})}\\)

2. Equilibrium values are energies
3. Convexity

.footnote-cite[.red[<sup>5</sup>] Lipparini, F.; Scalmani, G.; Mennucci, B.; Cancès, E.; Caricato, M.; Frisch, M. J. _J. Chem. Phys._ (2010), __133__, 014106]

???

---

## The Importance of Being Variational.red[<sup>6</sup>]

<p style="text-align:center;"><img src="images/variational-brain.png" style="width: 85%"></p>
<p style="clear: both;">

1. Avoids nonlinear couplings in the QM Hamiltonian
2. _Classical_ Hellmann-Feynman theorem
`$$
\frac{\mathop{}\!\mathrm{d}U_\mathrm{pol}}{\mathop{}\!\mathrm{d} F} =
\frac{\partial U_\mathrm{pol}}{\partial F}
+\cancelto{0}{\frac{\partial U_\mathrm{pol}}{\partial \color{Red}{\mathrm{p}}}\frac{\partial \color{Red}{\mathrm{p}}}{\partial F}}
$$`
3. **Valid for polarizable MM models**
4. Extended Lagrangian dynamics


.footnote-cite[.red[<sup>6</sup>] Lipparini, F.; Lagardère, L.; Raynaud, C.; Stamm, B.; Cancès, E.; Mennucci, B.; Schnieders, M.; Ren, P.; Maday, Y.; Piquemal, J.-P.
_J. Chem. Theory Comput._ (2015), __11__, 623]

???

Classical polarization energy of mixed discrete-continuum system.
`$$
U_\mathrm{pol}[\mathrm{p}] = \frac{1}{2}\mathrm{p}\color{Green}{\mathbb{V}}\mathrm{p} + \mathrm{p}\mathrm{s}
$$`
- p is the polarization degree of freedom
- s is the source term
- V is the classical interaction operator
- \\(\mathrm{p}\\) is an _independent_ degree of freedom
- _linear_ QM/classical coupling through

---

## Self-Consistent Field with PCM...

* .blue[Molecular electrostatic potential (MEP)]
  `$$
  \varphi(\mathbf{s}) = \sum_{\kappa\lambda}^{N_\mathrm{AO}} \color{Blue}{\varphi_{\kappa\lambda}(\mathbf{s})} D_{\lambda\kappa},
  \quad
  \color{Blue}{\varphi_{\kappa\lambda}(\mathbf{s}) }= \int\mathop{}\!\mathrm{d}\mathbf{r} \frac{-\chi^*_\kappa(\mathbf{r})\chi_\lambda(\mathbf{r}) }{|\mathbf{r} - \mathbf{s}|}
  $$`
* .red[Apparent surface charge (ASC)]
  `$$
   \color{Green}{\mathcal{T}}\color{Red}{\sigma(\mathbf{s})} = - \color{Green}{\mathcal{R}}\color{Blue}{\varphi(\mathbf{s})}
  $$`
* Fock matrix
  `$$
   f_{\kappa\lambda} =
   f^\mathrm{vac}_{\kappa\lambda} +
    (\color{Red}{\sigma(\mathbf{s})}, \color{Blue}{\varphi_{\kappa\lambda}(\mathbf{s})})_{
  \partial C}
  $$`

<p style="text-align:center;"><img src="images/scf.png" style="width: 100%"></p>
<p style="clear: both;">

???

---

## ...or indeed _any_ classical polarizable model

* .blue[Source term]
  `$$
  \mathrm{s} = \sum_{\kappa\lambda}^{N_\mathrm{AO}} \color{Blue}{\mathrm{s}_{\kappa\lambda}} D_{\lambda\kappa}
  $$`
* .red[Classical polarization]
  `$$
   \color{Green}{\mathbb{V}}\color{Red}{\mathrm{p}} = - \color{Blue}{\mathrm{s}}
  $$`
* Fock matrix
  `$$
   f_{\kappa\lambda} =
   f^\mathrm{vac}_{\kappa\lambda} +
    \color{Red}{\mathrm{p}} \color{Blue}{\mathrm{s}_{\kappa\lambda}}
  $$`

???

---

## Solvation in the Relativistic Regime

<p style="text-align:center;"><img src="images/paper-i.png" style="width: 90%"></p>
<p style="clear: both;">

<p style="text-align:center;"><img src="images/rela-mep.png" style="width: 80%"></p>
<p style="clear: both;">

???

Attacking the Hamiltonian axis

---

## Relativity and Solvation: EPR and pNMR parameters

<p style="text-align:left;"><img src="images/paper-iv.png" style="width: 110%"></p>
<p style="clear: both;">

`$$
\begin{align}
g_{uv}   = \frac{2c} {\langle \tilde{S}_v \rangle} \, \mathrm{ Tr} \left\{\mathbf{\Lambda}_{B_u} \mathbf{D}^{ (J_v)} \right\},
\quad&\quad
  \left(\mathbf{\Lambda}_{B_u}\right)_{\lambda\tau} = \frac{1}{2}
  \left\langle
       \mathbf{X}_\lambda |
       \left( \mathbf{r}_{G} \times \mathbf{\alpha} \right)_u
       |
       \mathbf{X}_\tau
  \right\rangle.
\\
A^M_{uv} = \frac{1} {\langle \tilde{S}_v \rangle} \,
           \mathrm{Tr} \left\{\mathbf{\Lambda}_{I^M_u} \mathbf{D}^{(J_v)}\right\},
\quad&\quad
\left(\mathbf{\Lambda}^{\mathrm{PN}}_{I^M_u}\right)_{\lambda\tau}
    = \gamma^M
    \left\langle\mathbf{X}_\lambda
       \left|
          \left(
          \frac{\mathbf{r}_{M} \times \mathbf{\alpha}}{r^{3}_{M}}
          \right)_u
       \right|
       \mathbf{X}_\tau
    \right\rangle
\\
\delta_M^\mathrm{para} &= \frac{\mu_e}{
  12\gamma^MkT} \,\text{Tr} \left( \mathbf{g}\mathbf{A}^\mathrm{T}_M \right)
\end{align}
$$`

???

- We want properties
- Expectation values
- Particularly sensitive to the environment

---

## Relativity and Solvation: EPR and pNMR parameters

<p style="text-align:center;"><img src="images/pNMR-tab3.png" style="width: 100%"></p>
<p style="clear: both;">

???

---

## Arbitrary Order Molecular Response Properties in Solution

<div class="imageWrapper">
  <img class="overlayImage" src="images/pcm_openrsp.png" style="width: 90%" align="middle">
--
  <img class="overlayImage" src="images/paper-v.png" style="width: 110%" align="left">
</div>

---

## Perturbed Quantum/Classical Polarizable SCF Energy Functional.red[<sup>8</sup>]

- Energy functional in \\(\mathbf{C}\\)-representation:
`$$
\tilde{G}(\tilde{\mathbf{C}}, \color{Red}{\tilde{\mathrm{p}}}, t) =
\tilde{E}(\tilde{\mathbf{C}}, t) +
\frac{1}{2} \color{Red}{\tilde{\mathrm{p}}} \color{Green}{\mathbb{V}} \color{Red}{\tilde{\mathrm{p}}}+
\color{Blue}{\tilde{\mathrm{s}}}(\tilde{\mathbf{C}}) \color{Red}{\tilde{\mathrm{p}}}
$$`
- Quasienergy in \\(\mathbf{C}\\)-representation:
`$$
\tilde{Q}(\tilde{\mathbf{C}}, \color{Red}{\tilde{\mathrm{p}}}, t) =
\tilde{G}(\tilde{\mathbf{C}}, \color{Red}{\tilde{\mathrm{p}}}, t) - \mathrm{i}\sum_I \langle \tilde{\phi}_I| \dot{\tilde{\phi}}_I \rangle
$$`
- Quasienergy Lagrangian in \\(\mathbf{C}\\)-representation:
`$$
\tilde{L}(\tilde{\mathbf{C}}, \color{Red}{\tilde{\mathrm{p}}}, \tilde{\mathbf{\lambda}}, t) =
\tilde{Q}(\tilde{\mathbf{C}}, \color{Red}{\tilde{\mathrm{p}}}, t) - \sum_{IJ}\tilde{\lambda}_{JI}(\langle \tilde{\phi}_I|\tilde{\phi}_J \rangle - \delta_{IJ})
$$`

.footnote-cite[.red[<sup>8</sup>] Thorvaldsen, A. J.; Ruud, K. _et al._, _J. Chem. Phys._ (2008), __129__, 214108]

???

---

**Variational parameters \\(\tilde{\mathbf{C}}\\), \\(\tilde{\mathrm{p}}\\), \\(\tilde{\mathbf{\lambda}}\\)**

`$$
\begin{align}
  \tilde{\mathbf{F}}\tilde{\mathbf{C}}\mathbf{\rho}
  -\mathrm{i}(\tilde{\mathbf{R}}\tilde{\mathbf{C}}+\tilde{\mathbf{S}}\dot{\tilde{\mathbf{C}}})\mathbf{\rho} -
  \tilde{\mathbf{S}}\tilde{\mathbf{C}}\mathbf{\rho}\tilde{\mathbf{\lambda}}\mathbf{\rho} & = 0
  \quad\quad\text{TDSCF} \\
  \mathbf{\rho}\tilde{\mathbf{C}}^\dagger \tilde{\mathbf{S}}\tilde{\mathbf{C}}\mathbf{\rho}
  & = \mathbf{\rho} \quad\quad\text{Idempotency} \\
  \color{Green}{\mathbb{V}}\color{Red}{\mathrm{p}} + \color{Blue}{\mathrm{s}} & = 0 \quad\quad\text{Classical polarization}
\end{align}
$$`

`$$
\tilde{\mathbf{F}}
= \tilde{\mathbf{F}}^\mathrm{vac} + \color{Red}{\tilde{\mathrm{p}}}\color{Blue}{\tilde{\mathrm{s}}}
$$`

**Generalized KS energy, \\(\mathbf{D}\\)-representation**

`$$
\tilde{\mathcal{G}}
\overset{\mathrm{tr}}{=}
(\tilde{\mathbf{h}} + \tilde{\mathbf{V}}^{t} + \frac{1}{2}\tilde{\mathbf{G}}^{\gamma}(\tilde{\mathbf{D}}) - \frac{\mathrm{i}}{2}\tilde{\mathbf{T}} )\tilde{\mathbf{D}} +
 \tilde{E}_\mathrm{xc}[\tilde{\rho}(\tilde{\mathbf{D}})] + h_\mathrm{nuc} +
 \frac{1}{2} \color{Red}{\tilde{\mathrm{p}}} \color{Green}{\mathbb{V}} \color{Red}{\tilde{\mathrm{p}}} + \color{Blue}{\tilde{\mathrm{s}}} \color{Red}{\tilde{\mathrm{p}}}\tilde{\mathbf{D}}
$$`

**Time-averaged quasienergy derivative**

`
$$
 \tilde{L}^a(\tilde{\mathbf{D}}, \color{Red}{\tilde{\mathrm{p}}}, t)
 \overset{\lbrace\text{Tr}\rbrace_T}{=}
 \tilde{\mathcal{G}}^{00,a} - \tilde{\mathbf{S}}^a\tilde{\mathbf{W}}
$$
`
???

---

## Exploiting Our Variational Advantages

`
$$
 \tilde{L}^a(\tilde{\mathbf{D}}, \color{Red}{\tilde{\mathrm{p}}}, t)
 \overset{\lbrace\text{Tr}\rbrace_T}{=}
 \tilde{\mathcal{G}}^{00,a} - \tilde{\mathbf{S}}^a\tilde{\mathbf{W}}
$$
`

**Fixed cavity approximation**
`
$$
\tilde{\mathcal{G}}^{00,a}(\tilde{\mathbf{D}}, \color{Red}{\tilde{\mathrm{p}}}, t)
\overset{\text{Tr}}{=}
(\tilde{\mathbf{h}}^a + \tilde{\mathbf{V}}^{t,a} +
\frac{1}{2}\tilde{\mathbf{G}}^{\gamma, a}(\tilde{\mathbf{D}}) +
\tilde{\mathbf{F}}^{\Omega_a}_\mathrm{xc}
-\frac{\mathrm{i}}{2}\tilde{\mathbf{T}}^a)\tilde{\mathbf{D}} + h_\mathrm{nuc}^a
+\color{Red}{\tilde{\mathrm{p}}}\color{Blue}{\tilde{\mathbf{\mathrm{s}}}^a}\tilde{\mathbf{D}}
$$
`

**Response functions**
`
$$
\langle\langle A; B \rangle\rangle_{\omega_b} =
    L^{ab}
 \overset{\lbrace\text{Tr}\rbrace_T}{=}
\mathcal{G}^{00,ab}
+\mathcal{G}^{10,a}\mathbf{D}^{b}
+\mathcal{G}^{01,a}\color{Red}{\mathrm{p}^{b}}
-\mathbf{S}^{ab}\mathbf{W} - \mathbf{S}^{a}\mathbf{W}^{b}
$$
`

**PCM response equation**
`
$$
\begin{alignat*}{2}
  \color{Green}{\mathbb{V}}\color{SkyBlue}{\mathrm{p}_\mathrm{H}^{b_N}}
  +\text{Tr}\,\color{Blue}{\mathrm{s}}\color{SkyBlue}{\mathbf{D}_\mathrm{H}^{b_N}} = 0
  \quad&\quad
  \color{Green}{\mathbb{V}}\color{Orange}{\mathrm{p}_\mathrm{P}^{b_N}}
  +\text{Tr}\,\color{Blue}{\mathrm{s}}\color{Orange}{\mathbf{D}_\mathrm{P}^{b_N}} =
  \color{Orange}{\mathbb{S}^{(n-1)}_\omega}
\end{alignat*}
$$
`

**\\(N\\)-the order response equations**
`
$$
\bbox[SkyBlue, 5pt]{\left[\mathbf{E}^{[2]} - \omega_{b_N}\mathbf{S}^{[2]} \right]
\mathbf{X}^{b_N}}
=\bbox[Orange, 5pt]{\mathbf{M}_\mathrm{RHS}^{b_N}}
$$
`

---

## Fluctuating Charge Model: One- and Two-Photon Absorption in Water

`
$$
\begin{aligned}
  F[\mathbf{q},\boldsymbol{\lambda}]
  &=
\sum_{\alpha,i} q_{\alpha i}\chi_{\alpha i} + \frac{1}{2}\sum_{\alpha,i} \sum_{\beta, j}q_{\alpha i} J_{\alpha i,\beta j}q_{\beta j} + \sum_{\alpha}\lambda_{\alpha} \left(\sum_{i} q_{\alpha i} - Q_{\alpha}\right) \\
  &= \mathbf{q}^\dagger\boldsymbol{\chi} + \frac{1}{2}\mathbf{q}^\dagger\mathbf{J}\mathbf{q} + \mathbf{q}^{\dagger}\boldsymbol{\lambda}
\end{aligned}
$$
`

- Classical sites with electronegativities (\\(\chi\\)) and hardnesses (\\(\eta\\))
- Multipliers (\\(\lambda\\)) to ensure charge conservation
- Polarization kernel:

`
$$
J_{\alpha i, \beta j} = \eta_{\alpha i}\delta_{\alpha i, \beta j} + \frac{\eta_{\alpha i, \beta j}}{\sqrt{1 + \eta_{\alpha i, \beta j}^{2}r^{2}_{\alpha i, \beta j}}}(1 - \delta_{\alpha i, \beta j})
$$
`

--

<p style="text-align:center;"><img src="images/paper-fq.png" style="width: 100%"></p>
<p style="clear: both;">

???

- Electronegativity equalization principle.
- Only parametrized for water.
- Coupling with QM same as continuum, as it only requires the MEP.

---
class: split-50-50

.column[
<p style="text-align:center;"><img src="images/OPA-methods-multiplot.png" style="width: 100%"></p>
]
.column[
<p style="text-align:center;"><img src="images/TPA-methods-multiplot-gaus-range.png" style="width: 50%"></p>
]
<p style="clear: both;">

<p style="text-align:center;"><img src="images/R6G_angle_label.png" style="width: 30%"></p>
<p style="clear: both;">

???

- OPA to the left, TPA to the right.

- Rhodamine 6G is a promising 2-photon dye. The `$S_{0} \rightarrow S_{2}$` transition is dark in OPA, but allowed in TPA.
- \\(25~\AA\\) sphere defines the MM region.
- Average over 100 snapshots.
- We can reproduce deltas of excitation energies for OPA. Absolute it's more difficult as S3 has significant CT character.
- For TPA comparison with experiment is challenging, as experiments are all over the place.

---
class: split-50-50

## PCM: Multiphoton Absorption in Solution

.column[
<p style="text-align:center;"><img src="images/PDNB_energy.png" style="width: 100%"></p>
<p style="clear: both;">

- Nonequilibrium solvation
- Odd-order MPA
- Even-order MPA
]
.column[
<p style="text-align:center;"><img src="images/PDNB_MPA.png" style="width: 50%"></p>
<p style="clear: both;">
]


---
class: split-50-50

## Coupled cluster theory

.column[
<p style="text-align:center;"><img src="images/excitations.svg" style="width: 100%"></p>
<p style="clear: both;">
]
.column[
`
$$
|\mathrm{CC}\rangle = \mathrm{e}^T |D_{\mathbf{0}}\rangle,
\quad
T = \sum_{u=1} T_u
$$
`
]

`
$$
T_u =
  \sum_{\substack{a_1\\ i_1}} t^{a_1}_{i_1}\tau_{i_1}^{a_1} +
  \frac{1}{4} \sum_{\substack{a_1, a_2 \\ i_1, i_2}} t^{a_1a_2}_{i_1i_2}\tau_{i_1i_2}^{a_1a_2} +
  \ldots +
  \frac{1}{(u!)^2} \sum_{\substack{a_1, a_2, \ldots, a_u \\ i_1, i_2, \ldots, i_u}}
t^{a_1a_2\ldots a_u}_{i_1i_2\ldots i_u}
\tau_{i_1i_2\ldots i_u}^{a_1a_2\ldots a_u}
$$
`

- Size-consistent and size-extensive
- Systematically improvable
- Polynomial scaling

???

- Construct the many-electron wavefunction as a multideterminantal expansion.
- Generate excited determinants from a single reference determinant.
- Excitations
- Wave operator is **not** unitary: `$\mathrm{e}^{-T} \neq \mathrm{e}^{T^\dagger}$`

---

## Coupled cluster equations: _unlinked_

`
$$
H_{\mathrm{N}}|\mathrm{CC} \rangle = \Delta E_{\mathrm{CC}}| \mathrm{CC} \rangle
$$
`

Assume _intermediate normalization_ `$ \langle D_{\mathbf{0}}| \mathrm{CC} \rangle = 1$`:

- Energy equation:

`
$$
\langle D_{\mathbf{0}} | H_{\mathrm{N}} |\mathrm{CC} \rangle =
\langle D_{\mathbf{0}} | H_{\mathrm{N}} (T_2 + T_1 + \frac{1}{2}T_1^{2})|D_\mathbf{0} \rangle =
\Delta E_{\mathrm{CC}}
$$
`

- Amplitude equations

`
$$
\langle D_{\mathbf{n}} | H_{\mathrm{N}} - \Delta E_{\mathrm{CC}} |\mathrm{CC} \rangle = 0
$$
`

- Size-extensive order-by-order, _not_ term-by-term

???

- Start from the normal-ordered Hamiltonian
- Project on reference determinant -> Get the energy (only singles and doubles needed)
- Project on excitation manifold -> Get amplitudes

---

## Coupled cluster equations: _linked_

`
$$
H_{\mathrm{N}}|\mathrm{CC} \rangle = \Delta E_{\mathrm{CC}}| \mathrm{CC} \rangle
$$
`

Assume _intermediate normalization_ `$ \langle D_{\mathbf{0}}| \mathrm{CC} \rangle = 1$`:

- Multiply by `$\mathrm{e}^{-T}$`:

`
$$
 \bar{H}_{\mathrm{N}} |D_\mathbf{0} \rangle =
\mathrm{e}^{-T} H_{\mathrm{N}} \mathrm{e}^{T} |D_\mathbf{0} \rangle = \Delta E_{\mathrm{CC}}| D_{\mathbf{0}} \rangle
$$
`

- Energy equation:

`
$$
\langle D_{\mathbf{0}} | \bar{H}_{\mathrm{N}} |\mathrm{CC} \rangle = \Delta E_{\mathrm{CC}}
$$
`

- Amplitude equations

`
$$
\langle D_{\mathbf{n}} | \bar{H}_{\mathrm{N}} | D_\mathbf{0} \rangle = 0
$$
`

- Size-extensive order-by-order _and_ term-by-term

???

- Start from the normal-ordered Hamiltonian
- Project on reference determinant -> Get the energy (only singles and doubles needed)
- Project on excitation manifold -> Get amplitudes

---

## From `$\bar{H}_\mathrm{N}$` to diagrams

- Baker-Campbell-Hausdorff expansion truncates

`
$$
\begin{aligned}
 \bar{H}_{\mathrm{N}} &=
       H_{\mathrm{N}} +
       [H_{\mathrm{N}}, T] +
       \frac{1}{2!}[[H_{\mathrm{N}}, T], T] \\ &+
       \frac{1}{3!}[[[H_{\mathrm{N}}, T], T], T] +
       \frac{1}{4!}[[[[H_{\mathrm{N}}, T], T], T], T]
\end{aligned}
$$
`

- Only connected terms contribute

<p style="text-align:center;"><img src="images/wick.svg" style="width: 60%"></p>
<p style="clear: both;">

???

- BCH truncates because the Hamiltonian is a two-body operator
- normal ordering, generalised Wick's theorem

---
count: false

## From `$\bar{H}_\mathrm{N}$` to diagrams

<p style="text-align:center;"><img src="images/bartlett-shavitt.jpg" style="width: 50%"></p>
<p style="clear: both;">

---
class: split-50-50

## Diagrammatic representation

.column[
### Excitors
<p style="text-align:left;"><img src="images/excitors.png" style="width: 60%"></p>
]
.column[
### Hamiltonian
<p style="text-align:left;"><img src="images/operators.png" style="width: 50%"></p>
]
<p style="clear: both;">

???

- Hole and particle lines
- Hamiltonian pieces

---
class: split-50-50

## Diagrammatic representation, contd.

.column[
### Excitors
<p style="text-align:left;"><img src="images/colored_excitors.png" style="width: 80%"></p>
]

.column[
### Hamiltonian
<p style="text-align:left;"><img src="images/operators.png" style="width: 50%"></p>
]
<p style="clear: both;">

---

<p style="text-align:left;"><img src="images/FT1.png" style="width: 90%"></p>
<p style="clear: both;">

<p style="text-align:left;"><img src="images/colored_FT1.png" style="width: 90%"></p>
<p style="clear: both;">

<p style="text-align:left;"><img src="images/FT1D0.png" style="width: 100%"></p>
<p style="clear: both;">

???

- We represent contractions with connecting lines
- Internal _vs_ external lines. The former are summed over.
- Loops, equivalent lines, signs, permutations

---

## Solving the coupled cluster equations

- Derive equations from diagrams
- Factorise nonlinear terms

<p style="text-align:left;"><img src="images/diagrams_and_equations.png" style="width: 100%"></p>
<p style="clear: both;">

???

- An explosion of terms, can be automated, but factorisation is nontrivial
- Dense algebra, nontrivial to parallelise

---
layout: false

## From deterministic to stochastic.red[<sup>7</sup>]

- Finite-different imaginary-time propagation

`
$$
| \Psi(\tau+\delta\tau) \rangle = [1 - \delta\tau(H - S)] | \Psi(\tau) \rangle
$$
`

- Modified CC Ansatz

`
$$
| \mathrm{CCMC} \rangle = N_{0}\mathrm{e}^{\frac{T^\prime}{N_0}} |D_{\mathbf{0}}\rangle
$$
`

where `$ \langle D_{\mathbf{0}}| \mathrm{CCMC}(\tau) \rangle = N_0(\tau)$`

- _Unlinked_ dynamic equation

`
$$
t_{\mathbf{n}}(\tau+\delta\tau) = t_{\mathbf{n}}(\tau) - \delta\tau
  \langle D_{\mathbf{n}} | H - S | \mathrm{CCMC}(\tau) \rangle
$$
`

.footnote-cite[.red[<sup>7</sup>] Thom, A. J. W. _Phys. Rev. Lett._ (2010), **105**, 263004]

???

- Sample the wavefunction
- Sample action of Hamiltonian
- Whereas results will be size-extensive and consistent, terms to sample may
not scale linearly with system size, even for perfectly noninteracting systems.

---
## Size extensivity and memory cost

<p style="text-align:center;"><img src="images/Be_replica_nstates_ccmc_fciqmc_fraction_hilbert_spaces.png" style="width: 100%"></p>
<p style="clear: both;">

???

- These are noninteracting replicas of Be atoms. Easy system.
- Cost is emphatically growing faster than linear.
- Why so? We are sampling terms that cancel out exactly in the end, but we're
still doing that part of the work and using memory for it.

---
## Learning new tricks from an old dog

- Imaginary-time Schrödinger equation.red[<sup>8</sup>]
`
$$
\frac{\mathrm{d}}{\mathrm{d}\tau} | \mathrm{CC} \rangle = -H | \mathrm{CC} \rangle
$$
`

- _Constant_ intermediate normalization
`
$$
\frac{\mathrm{d}t_{\mathbf{n}}}{\mathrm{d}\tau} = - \langle D_{\mathbf{n}} | \bar{H}_{\mathrm{N}}(\tau) | D_\mathbf{0} \rangle
$$
`

- Finite difference
`
$$
t_{\mathbf{n}}(\tau+\delta\tau) = t_{\mathbf{n}}(\tau) - \delta\tau\langle D_{\mathbf{n}} | \bar{H}_{\mathrm{N}}(\tau) | D_\mathbf{0} \rangle
$$
`

.footnote-cite[.red[<sup>2</sup>] Pigg, D. A.; Hagen, G.; Nam, H.; Papenbrock, T. _Phys. Rev. C Nucl. Phys._ (2012), **86**, 014308]

???

- Use the linked formulation, because terms are size-extensive by construction
- Rather than sample the commutator expansion, as done by Franklin _et al._,
enforce connectedness by using diagrammatic techniques.

---
## Diagrammatic Coupled Cluster Monte Carlo

`
$$
t_{\mathbf{n}}(\tau+\delta\tau) = t_{\mathbf{n}}(\tau) - \delta\tau\textcolor{red}{\langle D_{\mathbf{n}} | \bar{H}_{\mathrm{N}}(\tau) | D_\mathbf{0} \rangle}
$$
`

<div class="imageWrapper">
### Idea and Plan

- The _residual_ integral `$\textcolor{red}{\langle D_{\mathbf{n}} | \bar{H}_{\mathrm{N}}(\tau) | D_\mathbf{0} \rangle} = \sum \mathrm{diagrams}$`
- Build diagrams on the fly, _stochastically_
- _Stochastic_ rounding of amplitudes.

<p style="text-align:center;"><img src="images/propagate.png" style="width: 70%"></p>
<p style="clear: both;">
--
  <img class="overlayImage" src="images/paper-diagccmc.png" style="width: 110%" align="left">
</div>

???

- The red terms is a sum of diagrams. We can devise an algorithm to sample this
integral by building diagrams on the fly
- Stochastic rounding ensures sparse representation.
- Additionally this can be interpreted on par with deterministic approaches!

---

<p style="text-align:center;"><img src="images/fast-algorithm.gif" style="width: 100%"></p>
<p style="clear: both;">

.footnote-cite[Animation courtesy Katia Di Antonio]

---
## Size extensivity and memory cost, contd.

<p style="text-align:center"><img src="images/Be_replica_nstates_ccmc_fciqmc_diagccmc_fraction_hilbert_spaces.png" style="width: 100%"></p>
<p style="clear: both;">

---
## CPU cost

<p style="text-align:center;"><img src="images/Be_replica_nattempts.png" style="width: 100%"></p>
<p style="clear: both;">

---
## Locality

<p style="vertical-align:text-top"><img src="images/He_pentamer_nstates_dissoc.png" style="width: 100%"></p>
<p style="clear: both;">

---
class: split-50-50

## Acknowledgements

.column[
__Luca Frediani__ .cite[Hylleraas Centre, University of Tromsø]

__T. Daniel Crawford__ .cite[Virginia Tech]

__Charles J. C. Scott__ .cite[Univesity of Cambridge]

__Michal Repisky__ .cite[Hylleraas Centre, University of Tromsø]
]
.column[
__Filippo Lipparini__ .cite[Università di Pisa]

__Radovan Bast__ .cite[University of Tromsø]

__Alex J. W. Thom__ .cite[Univesity of Cambridge]

__Lori A. Burns__ .cite[Georgia Tech]
]

<p style="text-align:center;"><img src="images/NFR_merke_horisontal_oransje_eng.png" style="width: 60%"></p>
<p style="clear: both;">

---
name: last-page
template: inverse

# Thanks for your attention!

Slideshow created using [remark] and served using [cicero]

Slides available on [GitHub](https://github.com/robertodr/talk-lille)

Browse slides at [http://tinyurl.com/talk-lille](http://tinyurl.com/talk-lille)

[remark]: https://github.com/gnab/remark
[cicero]: https://github.com/bast/cicero
