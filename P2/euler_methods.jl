using Plots


function my_forward_euler(t0, Tf, Δt, y0, f, λ)
    N = Integer(Tf/Δt)  # N+1 total temporal nodes

    t = Vector{Float64}(undef, N+1)
    y = Vector{Float64}(undef, N+1)

    # fill in the initial condition:
    t[1] = t0
    y[1] = y0

    for n = 1:N # take N time steps
        y[n+1] = y[n] + Δt*f(t[n],y[n])
        t[n+1] = t[n] + Δt
    end
    
    return (t, y)
end

function my_backward_euler(t0, Tf, Δt, y0, f, λ)
    N = Integer(Tf/Δt)  # N+1 total temporal nodes

    t = Vector{Float64}(undef, N+1)
    y = Vector{Float64}(undef, N+1)

    # fill in the initial condition:
    t[1] = t0
    y[1] = y0

    for n = 1:N # take N time steps
        t[n+1] = t[n] + Δt
        y[n+1] = y[n] / (1-Δt*λ)
        # y[n+1] = y[n] + Δt*f(t[n+1],y[n])
    end
    
    return (t, y)
end

function test_equationRHS(t, y)
    return λ*y
end

function my_new_RHS(t, y)
    return -3y
end

λ = -3
Tf = 4
Δt = 2 #2/abs(λ) is limit of stability (λ = -3)
y0 = 17


(T, Y) = my_backward_euler(0, Tf, Δt, y0, my_new_RHS, λ)

scatter(T, Y, label = "approx", shape = :circle, color = :green)

tfine = 0:Δt/20:Tf
yexact = y0*exp.(-3.0*tfine)
# yexact = []
# for i in 0:Δt/20:Tf
#     push!(yexact,17*ℯ^(-3.0*i))
# end
#yexact = y0*exp.(λ*tfine)


plot!(tfine, yexact, label = "exact", color = :red)

