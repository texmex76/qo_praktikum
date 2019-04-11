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

ψt = []
t_start = 0
t_end = 1
t_step = 0.05
t_points = [t_start:t_step:t_end;]

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

push!(ψt, ψ0)
for t in t_points
    push!(ψt, propagate(t, ψt[end]))
end


fig = figure(figsize=(6, 4))
clf()
for i in ψt
    plot(pos_basis, i)
end
suptitle("Approx Time evolution from $(t_start) to $(t_end)")

gcf()
# savefig("infinite_well_approx.svg")
