using Pkg
Pkg.activate(".")

using CondensateDynamics, OrdinaryDiffEq, LSODA
import CondensateDynamics.V
using CUDA
using LaTeXStrings, Plots
import GR
using CUDA.CUFFT

gr()
GR.usecolorscheme(1)

## Solve the 1D harmonic oscillator
# problem with 1D-GPE 
L = (40.0,40.0,40.0)
N = (128,128,128)
sim = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)

@unpack_Sim sim
g = 0.0
equation = GPE_1D
iswitch = 1
x = X[1]
y = X[2]
z = X[3]
dV= volume_element(L, N)
reltol = 1e-4
tf = 10.0
maxiters = 1000
@. psi_0 = exp(-(x^2+y^2+z^2)/2) * exp(-im*x*10)

psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
initial_state = psi_0
display(typeof(psi_0)) 
fft!(psi_0)
alg = Tsit5()
@. V0 = 1/2 * (x^2+y^2+z^2)/10
V(x,y,z,t) = x^2
@pack_Sim! sim

# =============== method testing 
@info typeof(sim.T.Tkx!)
Tl = makeT(X,K, CuArray{Complex{Float64}},flags=nothing)
typeof(CUDA.CUFFT.plan_fft(CuArray(rand(45, 45, 45))))
typeof(Tl.Tkx)
display(sim.T.Tkx!*psi_0)



# ===================== simulation
sol = runsim(sim; info=false)
final = CuArray(sum(abs2.(sol[end]), dims=(2, 3)))
#display(fft!(final))
# NO IN PLACE OPERATIONS
final = Array((final))
psi_0 = Array(ifft(psi_0))
@info "final distribution norm squared: " ns(final, sim)

# p = plot(real.(x), abs2.(psi_0), label="initial")
display(abs2.(final)[:, 1, 1]|> collect)
p = plot(real.(Array(x)), abs2.(final[:, 1, 1]), label="final")
display(p)
ifft!(initial_state)

ht = heatmap(real.(x), sol.t, final)
display(ht)