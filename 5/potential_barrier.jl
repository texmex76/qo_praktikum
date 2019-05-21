using QuantumOptics
using PyPlot
pygui(true)

xmin = -100
xmax = 100
Npoints = 256
Δx = (xmax - xmin) / Npoints

b_position = PositionBasis(xmin, xmax, Npoints)
b_momentum = MomentumBasis(b_position)

V0 = 1. # Height of Barrier
d = 5 # Width of Barrier
L = 2 * d
function V_barrier(x)
    if x < -d/2 || x > d/2
        return 0.
    else
        return V0
    end
end
V = potentialoperator(b_position, V_barrier)

Txp = transform(b_position, b_momentum)
Tpx = transform(b_momentum, b_position)
Hkin = LazyProduct(Txp, momentum(b_momentum)^2/2, Tpx)

H = LazySum(Hkin, V)

xpoints = samplepoints(b_position)
pot = @. V_barrier(xpoints)

x0 = -35
sigma0 = 10
E = 6
m = 0.5
k = sqrt(2*m*E)
p0 = k

ψ_0 = gaussianstate(b_position, x0, p0, sigma0)
tmax = 30
timesteps = 25
T = collect(range(0.0, stop=tmax, length=timesteps))
tout, ψt = timeevolution.schroedinger(T, ψ_0, H)

fig = figure(figsize=(8,4))
ax1 = fig.add_subplot(121)
ax2 = ax1.twinx()

axp = fig.add_subplot(122)


# for i in 1:size(ψt)[1]
#     ax1.plot(xpoints, abs2.(ψt[i].data), "C0")
#     ax1.set_ylim(0, 0.2)
#     ax1.set_xlabel(L"x")
#     ax1.set_ylabel(L"|\psi(x)|^2")
#
#     ax2.plot(xpoints, pot, "C1")
#     ax2.set_ylim(0, 2)
#     ax2.set_ylabel(L"V(x)")
#     title("Gaussian state hitting potential step")
#     tight_layout(rect=[0, 0, 1, 1])
#
#     ax1.yaxis.label.set_color("C0")
#     ax2.yaxis.label.set_color("C1")
# end
# gcf()

idx = 25
ax1.plot(xpoints, normalize(abs2.(ψt[idx].data)), "C0")
ax1.set_ylim(0, 0.5)
ax1.set_xlabel(L"x")
ax1.set_ylabel(L"|\psi(x)|^2")

ax2.plot(xpoints, pot, "C1")
ax2.set_ylim(0, 2)
ax2.set_ylabel(L"V(x)")

ax1.yaxis.label.set_color("C0")
ax2.yaxis.label.set_color("C1")

mom = Tpx * ψt[idx]
ppoints = samplepoints(b_momentum)
axp.plot(ppoints, abs2.(mom.data))
axp.set_xlabel(L"p")
axp.set_ylabel(L"|\psi(p)|^2")
axp.yaxis.label.set_color("C0")

suptitle("Gaussian state hitting potential barrier")
tight_layout(rect=[0, 0, 1, .95])

gcf()

ψ_dens = abs2.(ψt[idx].data) ./ Δx
reflect = sum(ψ_dens[1:trunc(Int, size(ψ_dens)[1]/2)] .* Δx)
pass = sum(ψ_dens[trunc(Int, size(ψ_dens)[1]/2):end] .* Δx)
println("\nProbab pass: $(pass)")
println("Probab reflect: $(reflect)")
println("Probab total: $(reflect + pass)")

q = sqrt(2 * m * (E - V0))

trans_coeff = (2*k*q)^2 / ((2*k*q)^2 + (q^2-k^2)*sin(q*L)^2)
println("Analytic pass: $(trans_coeff)")
