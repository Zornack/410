#include("euler_methods.jl")
using LinearAlgebra
using SparseArrays
using Plots
# uₜ = kuₓₓ + F(x,t), 0 ≤ x ≤ 1, 0 ≤ t ≤ Tf
# u(x,0) = f(x)
#u(0,t) = u(1,t) = 0
#uₓ`=k[(uⱼ₊₁-2uⱼ+uⱼ₊₁) / Δx²)] + F(xⱼ,t) 
#u₁ = 0
#uₙ₊₁ = 0
#uⱼ(0) = f(xⱼ)

function heat(N, Δx)
    A = makeSparse(N,Δx)
    return A
end

function makeSparse(N, Δx)
    coefficient = 2/(Δx^2)
    A = spzeros(N-1,N-1)
    for i in 1:N-1
        A[i,i] = -2*coefficient
        if i == 1
            A[i,i+1] = 1*coefficient
        elseif i == N-1
            A[i,i-1] = 1*coefficient
        else
            A[i,i+1] = 1*coefficient
            A[i,i-1] = 1*coefficient
        end
    end
    return A
end

function b(t)
    return 2(π^2-1)*exp(-2t)*sin.(pi*x_int)
end

function f(x)
    return sin.(pi*x)
end

function matrix_forward_euler(t0, Tf, delT, y0, f, lambda)
    N = Integer(Tf/delT)  

    t = Vector{Float64}(undef, N+1)
    y = Vector{Float64}(undef, N+1)

    for n = 1:N  # take N time steps

        b = [0; F(t[n])]
    
        y[:, n+1] = y[:, n] + Δt*(A*y[:, n] + b)  # Forward Euler
        t[n+1] = t[n] + Δt
    end
    t[1] = t0
    y[1] = y0

    for n = 1:N 
        y[n+1] = y[n] + delT*f(t[n],y[n], lambda)
        t[n+1] = t[n] + delT
    end
    
    return (t, y)
end


Δx = 0.1
Tf = 1
Δt = 0.01
N = Integer(1/Δx)
M = Integer(Tf/Δt)
x = collect(range(0,1,step=Δx))
x_int = x[2:end-1]
yₒ=f(x_int)
y = zeros(N+1,M+1)
y[:,1] = f.(x)
y[1,1] = 0
y[N+1,1] = 0
t = 0.0
A = makeSparse(N,Δx)
for n = 1:M
    y[2:N,n+1] = y[2:N, n] + Δt*(A*y[2:N,n] + b(t))
    global t = t + Δt
end

# for n = 1:M+1
#     p = plot(x,y[:,n])
#     display(p)
#     sleep(.1)
# end