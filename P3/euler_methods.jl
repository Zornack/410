using Plots

"""
    my_forward_euler(t0, Tf, delT, y0, f, lambda) 
Computes forward eurler for a general IVP. Takes the initial time, the final time, the time step,
the initial y value, the ODE and the scalar. Returns an array of the time steps and an array of the y values. 
"""
function my_forward_euler(t0, Tf, delT, y0, f, lambda)
    N = Integer(Tf/delT)  

    t = Vector{Float64}(undef, N+1)
    y = Vector{Float64}(undef, N+1)


    t[1] = t0
    y[1] = y0

    for n = 1:N 
        y[n+1] = y[n] + delT*f(t[n],y[n], lambda)
        t[n+1] = t[n] + delT
    end
    
    return (t, y)
end
"""
    my_backward_euler(t0, Tf, delT, y0, f, lambda) 
Computes backwards eurler for a general IVP. Takes the initial time, the final time, the time step,
the initial y value, the ODE and the scalar. Returns an array of the time steps and an array of the y values. 
"""
function my_backward_euler(t0, Tf, delT, y0, f, lambda)
    N = Integer(Tf/delT)  

    t = Vector{Float64}(undef, N+1)
    y = Vector{Float64}(undef, N+1)


    t[1] = t0
    y[1] = y0

    for n = 1:N 
        y[n+1] = y[n] / (1-delT*lambda)
        t[n+1] = t[n] + delT
    end
    
    return (t, y)
end

"""
    ODE(y,lambda)
Returns the value of the ODE with scalar lambda at y. 
"""
function ODE(y, lambda)
    return lambda*y
end

lambda = -3
Tf = 4
delT = (1/7)*2/abs(lambda)
y0 = 17


(T, Y) = my_backward_euler(0, Tf, delT, y0, ODE, lambda)

scatter(T, Y, label = "Forward Euler Approximation", title="Forward Euler", color=:green)


tfine = 0:delT/20:Tf
yexact = 17*exp.(-3.0*tfine)


plot!(tfine, yexact, label = "y(t) = 17e¯³ᵗ")

xlabel!("t")
ylabel!("y")
