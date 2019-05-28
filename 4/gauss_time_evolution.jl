using QuantumOptics
using LinearAlgebra
using PyPlot
using PyCall
pygui(true)

xmin = -100
xmax = 100
Npoints = 256
Δx = (xmax - xmin) / Npoints
b_position = PositionBasis(xmin, xmax, Npoints)
b_momentum = MomentumBasis(b_position)
xpoints = samplepoints(b_position)

x0 = -5
p0 = 10
sigma = 5
psi0 = gaussianstate(b_position, x0, p0, sigma)
m = 0.5

x = position(b_position)
p = momentum(b_position)

H = p^2 / 2m

steps = 0.06
time_end = 4
T = [0:steps:time_end;]
tout, psi_t = timeevolution.schroedinger(T, psi0, H)

function psi_analytic(x, t)
    @. sqrt(sigma) / π^(1/4) * exp(im*p0*(x-x0)-im*p0^2*t/(2m)) / sqrt(sigma^2+im*t/m) * exp(-(x-x0-p0*t/m)^2/(2*(sigma^2+im*t/m)))
end

T = [0:steps:time_end;]
psi_ana = []
for t in T
    push!(psi_ana, psi_analytic(xpoints, t))
end

clf()
# fig = figure(figsize=(5,4))
fig = figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

for i in 1:size(T)[1]
    ax1.plot(xpoints, abs2.(psi_t[i].data), "C0")
    ax1.set_xlabel(L"x")
    ax1.set_ylabel(L"|\psi(x)|^2|")
    ax1.set_title("Simulation")
    # ax1.set_ylim(0, 0.15)

    ax2.plot(xpoints, abs2.(psi_ana[i]) * Δx, "C1")
    ax2.set_xlabel(L"x")
    ax2.set_ylabel(L"|\psi(x)|^2|")
    ax2.set_title("Analytic")
    # ax2.set_ylim(0, 0.15)
end

suptitle("Gaussian state")
# subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
tight_layout(rect=[0, 0, 1, .95])
gcf()


expect_x = real(expect(x, psi_t[1]))
println("\nx expected: $(x0 + p0 * t / m)")
println("x obtained: $(expect_x)")

varr = real(sum(abs2.(psi_t[1].data) .* xpoints.^2) - expect_x^2)
println("\nVariance x expected: $(sigma^2/2*(1+(t/(m*sigma^2))^2))")
println("Variance x obtained: $(varr)")

expect_p = real(expect(p, psi_t[1]))
println("\np expected: $(p0)")
println("p obtained: $(expect_p)")

b_momentum = MomentumBasis(b_position)
Tpx = transform(b_momentum, b_position)
mom = Tpx * psi_t[i]
ppoints = samplepoints(b_momentum)
varr_p = real(sum(abs2.(mom.data) .* ppoints.^2) - expect_p^2)
println("\nVariance p expected: $(1/(2*sigma^2))")
println("Variance p obtained: $(varr_p)")
