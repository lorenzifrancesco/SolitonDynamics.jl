# CondensateDynamics.jl dev notes

### CUFFT

CUFFT: in-place + planned fft is fastest

#### RAM saturation in plotting 
Due to 

### Soliton velocity and barrier energy
Tiling parameters are set by soliton units, described as follows


### Units
We work in normalized 
##### Healing length. Scattering energy

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