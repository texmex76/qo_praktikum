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

# # Momentum operator; second derivative
# for i = 1:pts
#     if i < pts
#         p_sqr[i,i] = -2 / step^2
#         p_sqr[i+1,i] = 1 / step^2
#         p_sqr[i,i+1] = 1 / step^2
#     else
#         p_sqr[i,i] = -2 / step^2
#     end
# end

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

# clf()
# for n in 1:3
#     lbl = "n=$(n)"
#     psi = eigenvalues[n] .+ abs2.(eigenvectors[:,n]./sqrt(step))
#     plot(pos_basis, psi, label=lbl)
# end
#
# xlabel(L"x")
# ylabel(L"|\psi_n(x)|^2")
# title("Infinite well wave function densities hand-engineered")
# legend(loc="upper right", bbox_to_anchor=(.95, .8))
#
# gcf()

# println("Points in basis: $(size(pos_basis)[1])")
# println("Energies hand-engineered:")
# for n in 1:3
#     println(@sprintf("%.5f", eigenvalues[n]))
# end
# println("")

ψ0 = eigenvectors[:,1]

ψt = []
tpoints = [0:0.1:2.5;]
for i in tpoints
    push!(ψt, exp(-im*H*i) * ψ0)
end

clf()
for i in ψt
    plot(pos_basis, i)
end

gcf()
