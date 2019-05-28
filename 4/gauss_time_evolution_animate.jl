using QuantumOptics
using LinearAlgebra
using PyPlot
using PyCall
pygui(false)
anim = pyimport("matplotlib.animation")

xmin = -10
xmax = 10
Npoints = 256
Δx = (xmax - xmin) / Npoints
b_position = PositionBasis(xmin, xmax, Npoints)
b_momentum = MomentumBasis(b_position)
xpoints = samplepoints(b_position)

x0 = .5
p0 = 2
sigma = 0.3
psi0 = gaussianstate(b_position, x0, p0, sigma)

p = momentum(b_position)
m = 0.5

H = p^2 / 2m

b_momentum = MomentumBasis(b_position)
Tpx = transform(b_momentum, b_position)
ppoints = samplepoints(b_momentum)

T = [0:0.02:.8;]
tout, psi_t = timeevolution.schroedinger(T, psi0, H)

function psi_analytic(x, t)
    @. sqrt(sigma) / π^(1/4) * exp(im*p0*(x-x0)-im*p0^2*t/(2m)) / sqrt(sigma^2+im*t/m) * exp(-(x-x0-p0*t/m)^2/(2*(sigma^2+im*t/m)))
end

psi_ana = []
for t in T
    push!(psi_ana, psi_analytic(xpoints, t))
end

fig = figure(figsize=(6,3))
# fig = figure()
clf()

function animate(i)
    clf()
    subplot(121)
    plot(xpoints, abs2.(psi_t[i].data), "C0")
    plot(xpoints, abs2.(psi_ana[i]) * Δx, "C1")
    xlabel(L"x")
    ylabel(L"|\psi(x)|^2")
    ylim(0, 0.15)
    title("Position basis")

    subplot(122)
    mom = Tpx * psi_t[i]
    plot(ppoints, abs.(mom.data))
    xlabel(L"p")
    ylabel(L"|\psi(p)|")
    title("Momentum basis")

    suptitle("Gaussian wave packet dispersion")
    tight_layout(rect=[0, 0, 1, .9])
end

frames = [1:1:size(T)[1];]
movie = anim.FuncAnimation(fig, animate, frames=frames, repeat=false, interval=200)
movie.save("gaussian_thingy.mp4")
