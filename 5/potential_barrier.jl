using QuantumOptics
using PyPlot

# rc("text", usetex=true)
# rc("text.latex", preamble= "\\usepackage{amsmath}")
# rc("font", family="serif")

pygui(true)

xmin = -300
xmax = 300
Npoints = 256
Δx = (xmax - xmin) / Npoints

b_position = PositionBasis(xmin, xmax, Npoints)
b_momentum = MomentumBasis(b_position)

V0 = 0.2 # Height of Barrier
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

x0 = -80
sigma0 = 30
# E = 0.605
E = 0.5948
m = 0.5
k = sqrt(2*m*E)
p0 = k

ψ_0 = gaussianstate(b_position, x0, p0, sigma0)
tmax = 300
timesteps = 25
T = collect(range(0.0, stop=tmax, length=timesteps))
tout, ψt = timeevolution.schroedinger(T, ψ_0, H)

idx = 25
clf()
fig = figure(figsize=(12,4))
ax1 = fig.add_subplot(131)
axp = fig.add_subplot(132)
ax_pass = fig.add_subplot(133)

ax2 = ax1.twinx()

for idx in trunc.(Int, collect(range(1, stop=25, length=5)))
    ax1.plot(xpoints, normalize(abs2.(ψt[idx].data)), "C0")
    ax1.set_ylim(0, 0.5)
    ax1.set_xlabel(L"x")
    ax1.set_ylabel(L"|\psi(x)|^2")

    ax2.plot(xpoints, pot, "C1")
    ax2.set_ylim(0, .5)
    ax2.set_ylabel(L"V(x)")

    ax1.yaxis.label.set_color("C0")
    ax2.yaxis.label.set_color("C1")

    mom = Tpx * ψt[idx]
    ppoints = samplepoints(b_momentum)
    axp.plot(ppoints, abs2.(mom.data))
    axp.set_xlabel(L"p")
    axp.set_ylabel(L"|\psi(p)|^2")
    axp.yaxis.label.set_color("C0")
end


pass_arr = []
trans_coeff_arr = []
energies = [0.2:0.01:1.3;]
for E in energies
k = sqrt(2*m*E)
p0 = k
ψ_0 = gaussianstate(b_position, x0, p0, sigma0)
tmax = 300
timesteps = 25
T = collect(range(0.0, stop=tmax, length=timesteps))
tout, ψt = timeevolution.schroedinger(T, ψ_0, H)

ψ_dens = abs2.(ψt[idx].data) ./ Δx
pass = sum(ψ_dens[trunc(Int, size(ψ_dens)[1]/2):end] .* Δx)
push!(pass_arr, pass)
q = sqrt(2 * m * (E - V0))
trans_coeff = (2*k*q)^2 / ((2*k*q)^2 + (q^2-k^2)^2*sin(q*L)^2)
push!(trans_coeff_arr, trans_coeff)
end

ax_pass.plot(energies, pass_arr, "C0", label="Simulation")
ax_pass.plot(energies, trans_coeff_arr, "C1", label="Analytic")
ax_pass.set_xlabel(L"E")
ax_pass.set_ylabel(L"T")
ax_pass.legend(bbox_to_anchor=(.45, .25), loc=2, borderaxespad=0.)


suptitle("Gaussian state hitting potential barrier")
# subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
tight_layout(rect=[0, 0, 1, .95])
gcf()
