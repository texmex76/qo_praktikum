using QuantumOptics
using PyPlot

# System Parameters
m = 1.
pot_strength = 1000
L = 2

# Position Basis
xmin = -2
xmax = 2
Npoints = 8
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

# Initial state
x0 = 0
p0 = 0
sigma0 = 0.1
Ψ0 = gaussianstate(b_position, x0, p0, sigma0);

# Time evolution
T = [0:0.1:1;]
tout, Ψt = timeevolution.schroedinger(T, Ψ0, H);

# Plot dynamics of particle density
x_points = particle.samplepoints(b_position)

clf()
figure(figsize=(6,3))
xlabel(L"x")
ylabel(L"| \Psi(t) |^2")

for i=1:length(T)
    Ψ = Ψt[i]
    n = abs.(Ψ.data).^2
    plot(x_points, n)
end
gcf()
