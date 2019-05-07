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

f(u,p,t) = -im*H*u
tspan = (0.0,3.0)

u0 = complex(ψ0)
prob = ODEProblem(f,u0,tspan)
ψt = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8,saveat=0.1)

u0_1 = complex(ψ0_1)
prob = ODEProblem(f,u0_1,tspan)
ψ_1t = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8,saveat=0.1)

u0_1 = complex(ψ0_2)
prob = ODEProblem(f,u0_1,tspan)
ψ_2t = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8,saveat=0.1)

t = 30
figure(figsize=(6, 5))
clf()
suptitle("Approx Time evolution at " * L"t=" * "$(t)/$(size(ψt)[2])")
subplots_adjust(top=0.85)

subplot(221)
plot(pos_basis, normalize(abs2.(ψt[t])))
xlabel(L"x")
ylabel(L"|\psi(x,t)|^2")
title(L"\psi_t")

subplot(222)
plot(pos_basis, normalize(abs2.((ψ_1t[t] + ψ_2t[t]) / 2)))
xlabel(L"x")
ylabel(L"|\psi(x,t)|^2")
title(L"(\psi_1+\psi_2)/2")

subplot(223)
plot(pos_basis, normalize(abs2.(ψ_1t[t])))
xlabel(L"x")
ylabel(L"|\psi(x,t)|^2")
title(L"\psi_1")

subplot(224)
plot(pos_basis, normalize(abs2.(ψ_2t[t])))
xlabel(L"x")
ylabel(L"|\psi(x,t)|^2")
title(L"\psi_2")

tight_layout(rect=[0, 0, 1, .95])
gcf()

# savefig("infinite_well_approx.svg")
