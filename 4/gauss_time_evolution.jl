using QuantumOptics
using LinearAlgebra
using PyPlot
using PyCall
pygui(false)

xmin = -10
xmax = 10
Npoints = 256
Δx = (xmax - xmin) / Npoints
b_position = PositionBasis(xmin, xmax, Npoints)
b_momentum = MomentumBasis(b_position)
xpoints = samplepoints(b_position)

x0 = .5
p0 = 0.4
sigma = 0.3
psi0 = gaussianstate(b_position, x0, p0, sigma)
m = 0.5

x = position(b_position)
p = momentum(b_position)
m = 0.5

H = p^2 / 2m

T = [0:0.02:.4;]
tout, psi_t = timeevolution.schroedinger(T, psi0, H)

function psi_analytic(x, t)
    @. sqrt(sigma) / π^(1/4) * exp(im*p0*(x-x0)-im*p0^2*t/(2m)) / sqrt(sigma^2+im*t/m) * exp(-(x-x0-p0*t/m)^2/(2*(sigma^2+im*t/m)))
end

steps = 0.02
T = [0:steps:0.4;]
psi_ana = []
for t in T
    push!(psi_ana, psi_analytic(xpoints, t))
end

t = 0.1
i = trunc(Int, t / steps + 1)

clf()
plot(xpoints, abs2.(psi_t[i].data), "C0")
plot(xpoints, abs2.(psi_ana[i]) * Δx, "C1")
xlabel(L"x")
ylabel(L"|\psi(x)|^2|")
ylim(0, 0.15)
gcf()

expect_x = real(expect(x, psi_t[i]))
println("\nx expected: $(x0 + p0 * t / m)")
println("x obtained: $(expect_x)")

varr = real(sum(abs2.(psi_t[i].data) .* xpoints.^2) - expect_x^2)
println("\nVariance x expected: $(sigma^2/2*(1+(t/(m*sigma^2))^2))")
println("Variance x obtained: $(varr)")

expect_p = real(expect(p, psi_t[i]))
println("\np expected: $(p0)")
println("p obtained: $(expect_p)")

b_momentum = MomentumBasis(b_position)
Tpx = transform(b_momentum, b_position)
mom = Tpx * psi_t[i]
ppoints = samplepoints(b_momentum)
varr_p = real(sum(abs2.(mom.data) .* ppoints.^2) - expect_p^2)
println("\nVariance p expected: $(1/(2*sigma^2))")
println("Variance p obtained: $(varr_p)")
