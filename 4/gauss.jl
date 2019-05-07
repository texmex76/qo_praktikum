using QuantumOptics
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
clf()
plot(xpoints, normalize(abs2.(psi.data)))
gcf()

x = position(b_position)
p = momentum(b_position)

expect_x = real(expect(x, psi))
println("\nx expected: $(x0)")
println("x obtained: $(expect_x)")

varr = real(sum(abs2.(psi.data) .* xpoints.^2) - expect_x^2)
println("\nVariance expected: $(sigma^2 / 2)")
println("Variance obtained: $(varr)")

expect_p = real(expect(p, psi))
println("\np expected: $(p0)")
println("p obtained: $(expect_p)")
