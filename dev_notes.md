# CondensateDynamics.jl dev notes

### CUFFT

CUFFT: in-place + planned fft is fastest

#### RAM saturation in plotting 
Due to 

### Soliton velocity and barrier energy
Tiling parameters are set by soliton units, described as follows
In SolitonBEC.jl we have the barrier set by the energy
$$
V = E_B \exp[-\frac{x^2}{2w^2}]
$$
with 
$$
\max{E_B} \approx 7500 \ \hbar
$$
corresponding to $1.2246 \ g_{1D}^2$ in harmonic units with $g_{1D} = -0.587$.

In [Helm-Gardiner], they use the potential 
$$
V = E_B \exp[-2\frac{x^2}{w'^2}]
$$
with a barrier parameter set by $q$ in soliton units
$$
 q = \sqrt{\frac{\pi}{2}} w'\frac{E_B}{gN}
$$
so that the integral of the potential is 
$$
qgN
$$
translated in our potential the $q$ param is inserted as
$$
q = \sqrt{2\pi} w\frac{E_B}{gN}
$$
$q$ is adimensional.
This implies 
$$
E_B = q \frac{gN}{\sqrt{2\pi} w}
$$
in harmonic units the relation is exactly the same.


## Units
We work in normalized harmonic units.

### Healing length. Scattering energy

Energy scale of interaction
$$
g_{1D}N|\psi_0|^2 = g_{1D}n_0
$$
this relation holds for the normalization 
$$
    \int dx  \ |\psi_0|^2 = 1
$$
since we want to represent the total number of particles. In the GPE, the coefficient g includes the total number of particle already, i.e.
$$
 g = g_{GPE-1D} = g_{1D}N 
$$
The length associated to interaction energy scale is the healing length
$$
\xi = \frac{\hbar}{\sqrt{2m g_{1D} n_0}}
$$

From this input we define the for a finite (localized) system, with a finite number of particles. In [Helm-Gardiner] the nonlinearity is quantified as 
$$
\tilde{g} = 2\hbar \omega_\perp |a_s| = g_{1D}
$$
So the corresponding units are: *soliton length unit*:
$$
\mathrm{L}_S =  \frac{\hbar^2}{m\tilde{g}N}
$$
And the *soliton time unit*:
$$
\mathrm{T}_S = \frac{\hbar^3}{m \tilde{g}^2N^2}
$$
and finally the *soliton energy unit*
$$
\mathrm{T}_S = \frac{m\tilde{g}N^2}{\hbar^2}
$$
Simulation include a total maximum velocity corresponding to one soliton velocity.

We need to express soliton units in unit of trap units to get the conversion factors and obtain the right velocity and barrier.
Thus
$$
\mathrm{L}_S =  \frac{1}{g_{1D}N} \mathrm{L}_H
$$
$$
\mathrm{T}_S =\frac{1}{g_{1D}^2N^2} \mathrm{T}_H  
$$
so we have
$$
\mathrm{V}_S= g_{GPE}\mathrm{V}_H
$$
Equivalently, energy is written in harmonic units as
$$
\mathrm{E}_S= g_{GPE}^2 \mathrm{E}_H 
$$


## Split-step solution of the 3D problem exploiting cylindrical symmetry
We can expand the GPE equation in cylindrical coordinates.
Assuming axial symmetry of the wavefunction we can neglect the azhimutal equation, and focus on the solution of the radial and axial equations. This leads to a 2D problem with great simplifications.
The normalized 3D GPE reads
$$
i\partial_t \psi = \left[-\dfrac{1}{2}\nabla^2 + V + g|\psi|^2\right]\psi
$$
Assume the decomposition
$$
\psi(r, x, \theta) = \phi(r, x) f(\theta)
$$
with 
$$
f(\theta) = \mathrm{const.}
$$
and
$$
V(r, x, \theta) = V(r, x)
$$
we have the radial-axial equation
$$
i\partial_t \phi = \left[-\frac{1}{2}\left(\frac{1}{r}\partial_r \left(r \partial_r \right) + \partial^2_z\right) + V + gf^2|\phi|^2\right]\phi
$$

> Remember the Laplacian in cylindrical coordinates
$$
    \nabla^2 = \frac{1}{r}\partial_r \left(r \partial_r \right) + \frac{1}{r^2} \partial^2_\theta + \partial^2_z
$$

How to implement a split-step approximation of this equation?
The most difficult computation to seems to be the radial (linear) component.

NB: in [Helm-Gardiner], they do not use $\mathrm{E}_S$ to set the barrier, they use the parameter $q$.

# June 2023: Major refactoring
We changed radically the workflow, using scripts as standalone units in which we import CondensateDynamics

##### Ground state issues
- 1D ground state convergence fails
 
1D-GPE with $N=1024, \, a_{tol}=10^{-4}, \, dt = 0.01 $ 
---> no advantage in decreasing tolerance or $dt$.
 Error on $\mu$ of around $1\%$. The only sensible improvement is in increasing $N$.

- 3D ground state convergence is bound to the initial state, too much

- NPSE_plus smash the RAM

And we don't know why. Using GC.gc() helps but cost too much time.


##### Tiling 
Almost ok, but needs to be tested before the armageddon-like simulations.

#### TODO
- [ ] 1D ground state convergence (rinuncio)
- [ ] Fix automatic ground state (wish)
- [ ] Fix CrankNicholson ground state (wish)

#### COLLAPSES
- The collapse gamma found for NPSE is found to be around 0.751, strange enough
This behavious depends strongly on the spatial discretization.
N=1024 -> 0.73qualcosa
N=2048 -> 0.751
N=4096 -> 0.74533

- NPSE_plus seems to have a higher threshold for collapse.
N=1024 -> gamma_c>0.75
> (are we sure we are even collapsing at all in NPSE_plus?) 


#### DO WE STILL HAVE CORRECTNESS PROBLEMS IN THE TRANSFORMS? 
 SEEMS SO, SINCE INCREMENTING THE STEPS IN THE SOLITON SIMULATION LEADS TO A DIFFERENT RESULT. _yes, but how do they SCALE?_
 - Parseval is well tested now.
- Maybe the transform domain is somehow corrupted and not taking into account the right nonlinearity?

---> working in kspace, we are "inverted"

### 27/06/2023
### Debug of strange behaviour of ground states
- [x] TEST: CrankNicholson -> seems to be without kinetic term: check initialization of tridiagonal matrices
- [x] - TEST: CrankNicholson initialized with soliton solution -> broken (seems without kinetic term)
- [x] TEST: SplitStep initialized with soliton solution -> After relatively few iterations, it converges to the wrong function (slightly below)  
- [x] BackwardEuler is shifting the ground state (not correct)

- [ ] Quick FIX CrankNicholson (nota available)

- [ ] Check correctness of the tranforms


**NB: we still do not have 3D CrankNicholson**
##### Curiosity
- [ ] obtain the Thomas-Fermi function simulating without kinetics

#### Why SplitStep is broken? Transforms??
We should check the kspace, the xspace should be ok (the transforms are inverse of each other)
- [ ] check using $dk = 2\pi/L$. It is used like that!
- [ ] justify theoretically check the correct version for $dk$ using sampling arguments

**NB: in our code the multiplicative constant to have the correct sampling of the transform are written in a counterintuitive way (DX, DK)**.

- in measures() the $dk$ is not corresponding to $2/pi$
- [ ] check also chempot functions

#### we are using a naive kvec!! 
- [x] try to replace the correct version, the tests survives
- [x] try to see if something changes in the soliton simulation: NOTHING CHANGES

There are two different $dk$, one used for norm in kspace - **modified** ($dk=\frac{2\pi}{dx \,N }$), and one which is the official one of the discretization **official** ($dk = \frac{2\pi}{L}$). 
WHY?

- the one used in kvecs (official) affects also the phyisics: it is contained in k^2. 
- In the FFT the only effect of the choice of $dk$ is the prefactor chosen for the normalization of the transform. In such a factor the "$dk$" is the modified



##### find the scaling of the error
Error scales strangely like $\sim 1/N^2$. This is expected for the SplitStep ??
- [ ] check in the literature: is this scaling expected?