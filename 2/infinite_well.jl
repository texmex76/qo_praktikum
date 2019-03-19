using LinearAlgebra
using PyPlot
using PyCall
using Printf

L = π
V0 = 1000
step = .05
m = .5

pos_basis = [-L:step:L;]
pts = size(pos_basis)[1]
V = zeros(pts, pts)
p_sqr = zeros(pts, pts)

# Function for evaluating potential
function pot(x)
    if abs(x) > L/2
        return V0
    else
        return 0
    end
end

# Potential
for i = 1:pts
    V[i,i] = pot(pos_basis[i])
end

# # Momentum operator
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
        p_sqr[i+2,i] = -3/2 / step^2
        p_sqr[i,i+2] = -3/2 / step^2
        p_sqr[i+3,i] = 1/90 / step^2
        p_sqr[i,i+3] = 1/90 / step^2
    elseif i < pts-1
        p_sqr[i,i] = -49/18 / step^2
        p_sqr[i+1,i] = 3/2 / step^2
        p_sqr[i,i+1] = 3/2 / step^2
        p_sqr[i+2,i] = -3/2 / step^2
        p_sqr[i,i+2] = -3/2 / step^2
    elseif i < pts
        p_sqr[i,i] = -49/18 / step^2
        p_sqr[i+1,i] = 3/2 / step^2
        p_sqr[i,i+1] = 3/2 / step^2
    else
        p_sqr[i,i] = -49/18 / step^2
    end
end

H_kin = -p_sqr/2m
H_kin = H_kin

H = H_kin + V


eigenvalues, eigenvectors = eigen(H)
println(eigenvalues)
println("")

clf()
n = 1
plot(pos_basis, eigenvectors[:,n])
suptitle("Wavefunction for n = " *@sprintf("%.0f", n))
gcf()
