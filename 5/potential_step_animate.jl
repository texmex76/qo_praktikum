using QuantumOptics
using PyPlot
using PyCall
pygui(false)
anim = pyimport("matplotlib.animation")

xmin = -50
xmax = 50
Npoints = 200

b_position = PositionBasis(xmin, xmax, Npoints)
b_momentum = MomentumBasis(b_position)

V0 = 1. # Height of Barrier
function V_barrier(x)
    if x < 0
        return 0.0
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

x0 = -15
sigma0 = 4
p0 = 1.5

ψ_0 = gaussianstate(b_position, x0, p0, sigma0)
tmax = 40
timesteps = 20
T = collect(range(0.0, stop=tmax, length=timesteps))
tout, ψt = timeevolution.schroedinger(T, ψ_0, H)

clf()
fig = figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()

function animate(i)
    ax1.clear()
    ax2.clear()

    ax1.plot(xpoints, abs2.(ψt[i].data), "C0")
    ax1.set_ylim(0, 0.2)
    ax1.set_xlabel(L"x")
    ax1.set_ylabel(L"|\psi(x)|^2")

    ax2.plot(xpoints, pot, "C1")
    ax2.set_ylim(0, 2)
    ax2.set_ylabel(L"V(x)")
    title("Gaussian state hitting potential step")
    tight_layout(rect=[0, 0, 1, 1])

    ax1.yaxis.label.set_color("C0")
    ax2.yaxis.label.set_color("C1")
end

frames = [1:1:size(ψt)[1];]
movie = anim.FuncAnimation(fig, animate, frames=frames, repeat=false, interval=200)
movie.save("potential_barrier.mp4")
