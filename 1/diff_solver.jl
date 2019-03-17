using DifferentialEquations

m = 1
ω = 1

dt = 0.1
t_end = 10
t_vec = [0:dt:t_end;]

function func(dy, y, p, t)
    dy[1] = y[2] / m
    dy[2] = -m * ω^2 * y[1]
end

prob = ODEproblem(func, y0, t0, t_end)
sol = solve(prob, Tsit5(), saveat=t_vec)
y = sol.u
