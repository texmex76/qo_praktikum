using QuantumOptics
using PyPlot
using LinearAlgebra

# System Parameters
m = .2
V0 = 10000
L = π

# Position Basis
xmin = -2.5
xmax = 2.5
Npoints = 32
b_position = PositionBasis(xmin, xmax, Npoints)
b_momentum = MomentumBasis(b_position)

# Function for evaluating potential
function pot(x)
    if abs(x) > L/2
        return V0
    else
        return 0
    end
end

T_px = particle.transform(b_momentum, b_position)
T_xp = dagger(T_px)

x = position(b_position)
p = momentum(b_momentum)

H_kin = LazyProduct(T_xp, p^2/2m, T_px)
V = potentialoperator(b_position, pot)
H = LazySum(H_kin, V)
H_dense = dense(LazySum(H_kin, V))

E, ψ_states = eigenstates((H_dense + dagger(H_dense)) / 2, 2)
ψ0 = (ψ_states[1] + ψ_states[2]) / sqrt(2)

# Time evolution
T = [0:0.1:1;]
tout, ψt = timeevolution.schroedinger(T, ψ0, H)

# Plot dynamics of particle density
x_points = particle.samplepoints(b_position)

# Plot eigenstates
# clf()
# plot(x_points, abs2.(ψ_states[1].data))
# plot(x_points, abs2.(ψ_states[2].data))
# gcf()

pygui(false)
clf()
figure(figsize=(6,3))
title("QO Toolbox Time evolution from $(T[1]) to $(T[end])")
xlabel(L"x")
ylabel(L"| \Psi(t) |^2")

for i=1:length(T)
    ψ = ψt[i]
    n = abs2.(real.(ψ.data))
    plot(x_points, n)
end
gcf()
