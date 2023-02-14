using Makie, LaTeXStrings
using CondensateDynamics
gr(colorbar=false,size=(600,150),legend=false,grid=false)

L = (4.0)
N = (256)
sim = InitSim(L, N)
@unpack_Sim sim

@pack_Sim
x = X[1]

runsim(sim; info=true)