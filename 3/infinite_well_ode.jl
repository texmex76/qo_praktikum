using LinearAlgebra
using PyPlot
using PyCall
using Printf
using DifferentialEquations
pygui(false)

L = π

Npoints = 32
basis = range(-L, stop = L, length = Npoints)
pos_basis = collect(basis)
Δx = 2 * L / Npoints

pts = size(pos_basis)[1]
V = zeros(pts, pts)
p_sqr = zeros(pts, pts)

V0 = 10000
m = .2

# Function for evaluating potential
function pot(x)
    if abs(x) > L/2
        return V0
    else
        return 0
    end
end

# Calculate potential
for i = 1:pts
    V[i,i] = pot(pos_basis[i])
end


# Momentum operator more Δxs
for i = 1:pts
    if i < pts-2
        p_sqr[i,i] = -49/18 / Δx^2
        p_sqr[i+1,i] = 3/2 / Δx^2
        p_sqr[i,i+1] = 3/2 / Δx^2
        p_sqr[i+2,i] = -3/20 / Δx^2
        p_sqr[i,i+2] = -3/20 / Δx^2
        p_sqr[i+3,i] = 1/90 / Δx^2
        p_sqr[i,i+3] = 1/90 / Δx^2
    elseif i < pts-1
        p_sqr[i,i] = -49/18 / Δx^2
        p_sqr[i+1,i] = 3/2 / Δx^2
        p_sqr[i,i+1] = 3/2 / Δx^2
        p_sqr[i+2,i] = -3/20 / Δx^2
        p_sqr[i,i+2] = -3/20 / Δx^2
    elseif i < pts
        p_sqr[i,i] = -49/18 / Δx^2
        p_sqr[i+1,i] = 3/2 / Δx^2
        p_sqr[i,i+1] = 3/2 / Δx^2
    else
        p_sqr[i,i] = -49/18 / Δx^2
    end
end

# Hamilton operators
H_kin = -p_sqr/2m / (π^2 / L^2)
H_kin = H_kin
H = H_kin + V

eigenvalues, eigenvectors = eigen(H)

ψ0_1 = eigenvectors[:,1]
ψ0_2 = eigenvectors[:,2]
ψ0 = (ψ0_1 + ψ0_2) / sqrt(2)

ψt = []
t_start = 0
t_end = 1
t_step = 0.05
t_points = [t_start:t_step:t_end;]

# for i in t_points
#     push!(ψt, exp(-im*H*i) * ψ0)
# end
#
# fig = figure(figsize=(6, 4))
# clf()
# for i in ψt
#     plot(pos_basis, i)
# end
# suptitle("Exact Time evolution from $(t_start) to $(t_end)")
#
# gcf()
# # savefig("infinite_well_exact.svg")

f(u,p,t) = -im*H*u*ψ0
u0= ψ0
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
using Plots
plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
plot!(sol.t, t->0.5*exp(1.01t),lw=3,ls=:dash,label="True Solution!")
