using QuantumOptics
using LinearAlgebra
using PyPlot
pygui(false)

xmin = -3
xmax = 3
Npoints = 64
Δx = (xmax - xmin) / Npoints
b_position = PositionBasis(xmin, xmax, Npoints)
b_momentum = MomentumBasis(b_position)

x0 = .5
p0 = 0.4
sigma = 0.3
psi = gaussianstate(b_position, x0, p0, sigma)

xpoints = samplepoints(b_position)

x = position(b_position)
p = momentum(b_position)

expect_x = real(expect(x, psi))
println("\nx expected: $(x0)")
println("x obtained: $(expect_x)")

varr = real(sum(abs2.(psi.data) .* xpoints.^2) - expect_x^2)
println("\nVariance x expected: $(sigma^2 / 2)")
println("Variance x obtained: $(varr)")

expect_p = real(expect(p, psi))
println("\np expected: $(p0)")
println("p obtained: $(expect_p)")

b_momentum = MomentumBasis(b_position)
Tpx = transform(b_momentum, b_position)
mom = Tpx * psi
ppoints = samplepoints(b_momentum)
varr_p = real(sum(abs2.(mom.data) .* ppoints.^2) - expect_p^2)
println("\nVariance p expected: $(1/(2*sigma^2))")
println("Variance p obtained: $(varr_p)")

mom_calculated = @. sqrt(sigma) / π^(1/4) * exp(-(ppoints - p0)^2 * sigma^2 / 2) * exp(-im * x0 * ppoints)

figure(figsize=(6, 5))
clf()
# suptitle("Approx Time evolution at " * L"t=" * "$(t)/$(size(ψt)[2])")
subplots_adjust(top=0.85)

subplot(221)
plot(xpoints, abs2.(psi.data))
xlabel(L"x")
ylabel(L"|\psi(x)|^2")
title("Position basis")

subplot(222)
plot(ppoints, abs.(mom.data))
xlabel(L"p")
ylabel(L"|\psi(p)|")
title("Momentum basis")

subplot(223)
plot(xpoints, abs2.(psi.data))
xlabel(L"x")
ylabel(L"|\psi(x,t)|^2")
title(L"\psi_1")

subplot(224)
plot(ppoints, abs.(mom_calculated))
xlabel(L"p")
ylabel(L"|\psi(p)|")
title("Mom calculated by hand")

tight_layout(rect=[0, 0, 1, .95])
gcf()

# savefig("infinite_well_approx.svg")
