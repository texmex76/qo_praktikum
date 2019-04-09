using LinearAlgebra
using PyPlot
using PyCall
using Printf
pygui(false)

L = π

Npoints = 32
basis = range(-L, stop = L, length = Npoints)
pos_basis = collect(basis)
step = 2 * L / Npoints

pts = size(pos_basis)[1]
V = zeros(pts, pts)
p_sqr = zeros(pts, pts)

V0 = 10000
m = .1999 / sqrt(step)

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


# Momentum operator more steps
for i = 1:pts
    if i < pts-2
        p_sqr[i,i] = -49/18 / step^2
        p_sqr[i+1,i] = 3/2 / step^2
        p_sqr[i,i+1] = 3/2 / step^2
        p_sqr[i+2,i] = -3/20 / step^2
        p_sqr[i,i+2] = -3/20 / step^2
        p_sqr[i+3,i] = 1/90 / step^2
        p_sqr[i,i+3] = 1/90 / step^2
    elseif i < pts-1
        p_sqr[i,i] = -49/18 / step^2
        p_sqr[i+1,i] = 3/2 / step^2
        p_sqr[i,i+1] = 3/2 / step^2
        p_sqr[i+2,i] = -3/20 / step^2
        p_sqr[i,i+2] = -3/20 / step^2
    elseif i < pts
        p_sqr[i,i] = -49/18 / step^2
        p_sqr[i+1,i] = 3/2 / step^2
        p_sqr[i,i+1] = 3/2 / step^2
    else
        p_sqr[i,i] = -49/18 / step^2
    end
end

# Hamilton operators
H_kin = -p_sqr/2m / (π^2 / L^2)
H_kin = H_kin
H = H_kin + V

eigenvalues, eigenvectors = eigen(H)


ψ0 = eigenvectors[:,1]

ψt = []
t_start = 0
t_end = 3
tpoints = [t_start:0.1:t_end;]
for i in tpoints
    push!(ψt, exp(-im*H*i) * ψ0)
end

fig = figure(figsize=(6, 4))
clf()
for i in ψt
    plot(pos_basis, i)
end
suptitle("Time evolution from $(t_start) to $(t_end)")

gcf()
