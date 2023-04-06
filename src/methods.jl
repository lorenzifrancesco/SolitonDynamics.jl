## ==== Arrays


V(x, t) = 0.0
V(x, y, t) = 0.0
V(x, y, z, t) = 0.0

# array methods
# generate arrays 
"""
    x = xvec(λ,N)

Create `x` values with correct periodicity for box specified by length `λ`, using `N` points.
"""
xvec(L,N) = LinRange(-L/2,L/2,N) |> collect

"""
    k = kvec(λ,N)

Create `k` values with correct periodicity for box specified by length `L` for number of points `N`.
"""
function kvec(L,N)
    # @assert iseven(N)
    # nk = 0:Int(N/2)
    # k = [nk[1:end-1];-reverse(nk[2:end])]*2*π/λ
    k = fftfreq(N)*N*2*π/L
    return k
end

"""
    X = xvecs(L,N)

Create a tuple containing the spatial coordinate array for each spatial dimension.
"""
function xvecs(L,N)
    X = []
    for (λ,ν) in zip(L,N)
        x = xvec(λ,ν)
        push!(X,x)
    end
    return X |> Tuple
end

"""
    K = kvecs(L,N)

Create a tuple containing the spatial coordinate array for each spatial dimension.
"""
function kvecs(L,N)
    K = []
    for (λ,ν) in zip(L,N)
        k = kvec(λ,ν)
        push!(K,k)
    end
    return K |> Tuple
end

"""
    k² = k2(K)

Create the kinetic energy array `k²` on the `k`-grid defined by the tuple `K`.
"""
function k2(K, A)
    kind = Iterators.product(K...)
    tmp::A = map(k-> sum(abs2.(k)),kind)
    return tmp
end

"""
return spatial volume element
"""
function volume_element(L, N)
    dV=1
    for i in eachindex(L)
        dV *= L[i]/N[i]
    end
    return dV
end


"""
compute normalization
"""
function ns(psi, sim)
    return sum(abs2.(psi)) * sim.dV
end

"""
compute normalization in k-space
"""
function nsk(psi, sim)
    return sum(abs2.(psi)) / sim.Vol
end

"""
compute norm squared of a region
"""
function ns(psi, sim, mask)
    return sum(abs2.(psi).* mask) * sim.dV
end

"""
return σ^2(ψ) of the NPSE
"""
function init_sigma2(g::Float64)
    function sigma2(psi::ComplexF64)
        result = psi
        try
            result = sqrt(maximum([1 + g * abs2(psi), 0.05]))
        catch  err
            if isa(err, DomainError)
                result = NaN
                throw(NpseCollapse(g * maximum(abs2.(psi))))
            else
                throw(err)
            end
        end
        return result
    end
    return sigma2
end

# """
# check notes, Sigma is sigma^2
# """
# function deSigmaLambda!(dsigma, sigma, psi, t)
#     return 
# end

# """
# return the sigma2 function to be called in the ODE loop
# """
# function init_sigma2_ode(g::Float64, x::Array{Float64})
#     function sigma2(psi::Array{ComplexF64}, sigma2_0::Float64, lambda::Float64)
#         problem = ODEProblem(deSigmaLambda!, Array([sigma2_0, lambda]), (x[1], x[end]), psi)
#         sol = solve(problem,
#                     alg=BS3(),
#                     dense=false,
#                     maxiters=2000,
#                     progress=true, 
#                     )
#         result = sol.u
#         try
#             tmp = sqrt.(result)
#         catch  err
#             if isa(err, DomainError)
#                 result = NaN
#                 throw(NpseCollapse(NaN))
#             else
#                 throw(err)
#             end
#         end
#         return result
#     end
#     return sigma2
# end

"""
chemical potential in a given configuration
"""
function chempotk(psi, sim)
    @unpack ksquared,dV,V0,Vol,g = sim 
    mu = 1/Vol * sum(1/2 * ksquared .* abs2.(psi))
    tmp = xspace(psi, sim)
    mu += dV * sum((V0 + g*abs2.(tmp)) .* abs2.(tmp))
    mu *= 1/ns(tmp, sim) 
    #mu += 1 # add one transverse energy unit (1D-GPE case)
    return mu
end

"""
chemical potential in a given configuration
"""
function chempot(psi, sim)
    @unpack ksquared,dV,V0,Vol,g = sim 
    mu = dV * sum((V0 + g*abs2.(psi)) .* abs2.(psi))
    tmp = kspace(psi, sim)
    mu += 1/Vol * sum(1/2 * ksquared .* abs2.(tmp))
    mu *= 1/nsk(tmp, sim) 
    #mu += 1 # add one transverse energy unit (1D-GPE case)
    return mu
end


"""
    X,K,dX,dK = makearrays(L,N)

Create all `x` and `k` arrays for box specified by tuples `L=(Lx,...)` and `N=(Nx,...)`.
For convenience, differentials `dX`, `dK` are also reaturned. `L` and `N` must be tuples of equal length.
"""
function makearrays(L,N)
    @assert length(L) == length(N)
    X = xvecs(L,N)
    K = kvecs(L,N)
    dX = Float64[]; dK = Float64[]
    for j ∈ eachindex(X)
        x = X[j]; k = K[j]
        dx = x[2]-x[1]; dk = k[2]-k[1]
        push!(dX,dx)
        push!(dK,dk)
    end
    dX = dX |> Tuple
    dK = dK |> Tuple
    return X,K,dX,dK
end

"""
    A = crandn_array(M)

Make placeholder `2x2x...` complex `randn()` array of `M` dimensions."""
function crandn_array(N,T) 
    a::T = rand(N...) |> complex
    return a
end

"""
    A = crandnpartition(D,M)

Make placeholder ArrayPartition vector of length `M`, containing `2x2x...` rank D complex matrices.
"""
function crandnpartition(N,A)
    a = crandn_array(N,A)
    args = []
    for j = 1:N
        push!(args,a)
    end
    return ArrayPartition(args...)
end

## ==== Transforms

"""
    Dx,Dk = dfft(x,k)

Measures that make `fft`, `ifft` 2-norm preserving.
Correct measures for mapping between `x`- and `k`-space.
"""
function dfft(x,k)
    dx = x[2]-x[1]
    Dx = dx
    Dk = 1/Dx
    return Dx, Dk
end

"""
    DX,DK = dfftall(X,K)

Evalutes tuple of measures that make `fft`, `ifft` 2-norm preserving for each
`x` or `k` dimension.
"""
function dfftall(X,K)
    M = length(X)
    DX = zeros(M); DK = zeros(M)
    for i ∈ eachindex(X)
        DX[i],DK[i] = dfft(X[i],K[i])
    end
    return DX,DK
end

function xspace(ϕ,sim)
    @unpack T = sim
    return T.Tkx*ϕ
end

function xspace!(ψ,sim)
    @unpack T = sim
    T.Tkx!*ψ
    return nothing
end

function kspace(ψ,sim)
    @unpack T = sim
    return T.Txk*ψ
end

function kspace!(ψ,sim)
    @unpack T = sim
    T.Txk!*ψ
    return nothing
end

# """
#     CuArrays versions
# """
# function xspace(ϕ::CuArray,sim)
#     @unpack T = sim
#     return CUDA.CUFFT.:(*)(T.Tkx::AbstractFFTs.Plan{T}, ϕ)
# end

# function xspace!(ψ::CuArray,sim)
#     @unpack T = sim
#     T.Tkx!*ψ
#     return nothing
# end

# function kspace(ψ::CuArray,sim)
#     @unpack T = sim
#     return T.Txk*ψ
# end

# function kspace!(ψ::CuArray,sim)
#     @unpack T = sim
#     T.Txk!*ψ
#     return nothing
# end

"""
    definetransforms(funcs,args,meas,kwargs)
"""
function definetransforms(funcs,args,meas;kwargs=nothing)
    trans = []
    if kwargs === nothing
        for (fun,arg) in zip(funcs,args)
            push!(trans, fun(arg...))
        end
    else
        for (fun,arg) in zip(funcs,args)
            push!(trans, fun(arg...,flags=kwargs))
        end
    end
    
    return meas.*trans
end

"""
    T = makeT(X,K,j)
"""
function makeT(X,K,T::Type{Array{ComplexF64}};flags=FFTW.MEASURE)
    FFTW.set_num_threads(Sys.CPU_THREADS)
    D = length(X)
    N = length.(X)
    DX,DK = dfftall(X,K)
    dμx = prod(DX)
    dμk = prod(DK)
    psi_test = ones(N...) |> complex

    trans = (plan_fft,plan_fft!,plan_ifft,plan_ifft!)
    meas = (dμx,dμx,dμk,dμk)
    args = ((psi_test,),(psi_test,),(psi_test,),(psi_test,))
    Txk,Txk!,Tkx,Tkx! = definetransforms(trans,args,meas;kwargs=flags)
    return Transforms{D,N,T}(Txk,Txk!,Tkx,Tkx!)
end

"""
    T = makeT(X,K,j) in CUDA
"""
function makeT(X,K,T::Type{CuArray{ComplexF64}};flags=FFTW.MEASURE)
    D = length(X)
    N = length.(X)
    DX,DK = dfftall(X,K)
    dμx = prod(DX)
    dμk = prod(DK)
    psi_test = CuArray(ones(N...)) |> complex
    trans = (CUDA.CUFFT.plan_fft,CUDA.CUFFT.plan_fft!,CUDA.CUFFT.plan_ifft,CUDA.CUFFT.plan_ifft!)
    meas = (dμx,dμx,dμk,dμk)
    args = ((psi_test,),(psi_test,),(psi_test,),(psi_test,))
    Txk,Txk!,Tkx,Tkx! = definetransforms(trans,args,meas)
    return GPUTransforms{D,N,T}(Txk,Txk!,Tkx,Tkx!)
end