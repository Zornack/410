using LinearAlgebra
using SparseArrays
using Plots
# U_t = k*U_xx + F(t,x) on o ≤ x ≤, 0 ≤ t ≤ Tf


k = 1

dx = 0.1


λ = (1/2)*(1/2k)
dt = round(λ*dx^2, sigdigits = 10)
dt = 0.01
N = Integer(1/dx)
x = collect(range(0,1,step=dx))
x_int = x[2:end-1]
Tf = 1
M = Integer(Tf/dt)

A = zeros(N-1,N-1)

for i=1:N-1
    A[i,i] = -2
    if i != N-1
        A[i+1,i] = 1
        A[i,i+1] = 1
    end
end

A = (1/(dx^2))*A

function F(t,x)
    return exp(-t) * sin(x)
end

function f(x)
    return sin(pi*x)
end

function b(t)
    return F.(t,x_int)
end

y0 = f.(x)
t = zeros(M+1)
y = zeros(N+1,M+1)
y[:,1] = f.(x)
y[1,1] = 0
y[N+1,1] = 0
t[1] = 0+dt
for n = 1:M
    y[2:N,n+1] = y[2:N,n] + dt*(A*y[2:N,n]+b(t[n]))
    t[n+1] = t[n] + dt
end

# for n = 1:M
#     y[:,n+1] = y[:,n] + dt*(A*F.(t[n],x_int) + b(t[n]))
#     t[n+1] = t[n] + dt
# end

# for i in 1:M+1
#     p = plot(x,y[:,i],ylim=(0.0,1.0))
#     display(p)
#     sleep(0.01)
# end
