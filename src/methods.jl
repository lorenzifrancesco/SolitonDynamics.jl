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
    return sum(abs2.(psi) * sim.dV)
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
    Dx = dx/sqrt(2*pi)
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

"""
    ψ = xspace(ϕ,sim)

Transform from `k`- to `x`-space using transforms packed into `sim`.
"""
function xspace(ϕ,sim)
    @unpack T = sim
    return T.Tkx*ϕ
end

"""
    xspace!(ϕ,sim)

Mutating transform from `k`- to `x`-space using transforms packed into `sim`.
"""
function xspace!(ψ,sim)
    @unpack T = sim
    T.Tkx!*ψ
    return nothing
end

"""
    kspace(ψ,sim)
"""
function kspace(ψ,sim)
    @unpack T = sim
    return T.Txk*ψ
end

"""
    kspace!(ψ,sim)
"""
function kspace!(ψ,sim)
    @unpack T = sim
    T.Txk!*ψ
    return nothing
end

"""
    definetransforms(funcs,args,meas,kwargs)
"""
function definetransforms(funcs,args,meas,kwargs)
    trans = []
    for (fun,arg) in zip(funcs,args)
        push!(trans, fun(arg...,flags=kwargs))
    end
    return meas.*trans
end

"""
    T = makeT(X,K,j)
"""
function makeT(X,K,T;flags=FFTW.MEASURE)
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
    Txk,Txk!,Tkx,Tkx! = definetransforms(trans,args,meas,flags)
    return Transforms{D,N,T}(Txk,Txk!,Tkx,Tkx!)
end