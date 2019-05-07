using LinearAlgebra
using PyPlot
using PyCall
using Printf
using FFTW
anim = pyimport("matplotlib.animation")
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

# Momentum operator
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
H = H_kin + V

eigenvalues, eigenvectors = eigen(H)
ψ0_1 = eigenvectors[:,1]
ψ0_2 = eigenvectors[:,2]
ψ0 = (ψ0_1 + ψ0_2) / sqrt(2)

p_arr = [-π/Δx:(2*π/(2*L)):π/Δx-1;]

Δt = 0.001
exp_p_left = @. exp(-im * p_arr^2 * Δt / (2 * m * 2))
exp_p_right = @. exp(im * p_arr^2 * Δt / (2 * m * 2))
exp_pot = exp.(-im * diag(V) * Δt)

# ψ_p = fftshift(fft(ψ_x))
# ψ_x = ifft(ifftshift(ψ_p))

function propagate(t, ψ)
    temp = exp_p_right .* fftshift(fft(ψ)) # now in ms
    M = ceil(Int, t / Δt)
    for i in 1:M
        temp = exp_p_left .* temp # now in ms
        temp = exp_pot .* ifft(ifftshift(temp)) # now in ps
        temp = fftshift(fft(temp)) # now in ms
    end
    temp = exp_p_left .* temp
    return ifft(ifftshift(temp)) # now in ps
end

t_start = 0
t_end = 3
t_step = 0.05
t_points = [t_start:t_step:t_end;]

ψt = []
push!(ψt, ψ0)
for t in t_points
    push!(ψt, propagate(t, ψt[end]))
end

ψ_1t = []
push!(ψ_1t, ψ0_1)
for t in t_points
    push!(ψ_1t, propagate(t, ψ_1t[end]))
end

ψ_2t = []
push!(ψ_2t, ψ0_2)
for t in t_points
    push!(ψ_2t, propagate(t, ψ_2t[end]))
end

t = 60
figure(figsize=(6, 5))
clf()
suptitle("Approx Time evolution at " * L"t=" * "$(t)/$(size(ψt)[1])")
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
