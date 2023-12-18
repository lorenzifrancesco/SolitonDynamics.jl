## ==== Arrays


V(x, t) = 0.0
V(x, y, t) = 0.0
V(x, y, z, t) = 0.0

# array methods
# generate arrays 
xvec(L,N) = LinRange(-L/2,L/2,N+1)[1:end-1]

function kvec(L,N)
    # @assert iseven(N)
    # nk = 0:Int(N/2)
    # k = [nk[1:end-1];-reverse(nk[2:end])]*2*π/λ
    k = fftfreq(N)*N*2*π/L 
    return k
end

function xvecs(L,N)
    X = []
    for (λ,ν) in zip(L,N)
        x = xvec(λ,ν)
        push!(X,x)
    end
    return X |> Tuple
end

function kvecs(L,N)
    K = []
    for (λ,ν) in zip(L,N)
        k = kvec(λ,ν)
        push!(K,k)
    end
    return K |> Tuple
end

function k2(K, A)
    kind = Iterators.product(K...)
    tmp::A = map(k-> sum(abs2.(k)),kind)
    return tmp
end

function volume_element(L, N)
    dV=1
    for i in eachindex(L)
        dV *= L[i]/(N[i])
    end
    return dV
end

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
# TODO 1
# --> into transforms
function dfft(x,k)
    dx = x[2]-x[1]
    Dx = dx 
    Dk = 1/Dx
    return Dx, Dk
end

# TODO 2
# --> into nsk, ns
function measures(L, N)
    dX = L ./ (N)
    dK = 1 ./ (dX .* N)
    return  prod(dX), prod(dK)
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
    FFTW.set_num_threads(1)
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