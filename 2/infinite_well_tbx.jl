using QuantumOptics
using PyPlot
using PyCall
pygui(true)

# Parameters
L = π
Npoints = 64
V0 = 1000

b_position = PositionBasis(-L, L, Npoints)
step = 2 * L / Npoints
m = .1999 / sqrt(step)

p = momentum(b_position)

function potential(x)
    if abs(x) > L/2
        return V0
    else
        return 0
    end
end

V = potentialoperator(b_position, potential)

H_kin = p^2/2m / (π^2 / L^2)
H = H_kin + V

E, states = eigenstates((H + dagger(H))/2, 3);

xpoints = samplepoints(b_position)
clf()
for n in 1:3
    lbl = "n=$(n)"
    plot(xpoints, E[n] .+ abs2.(states[n].data./sqrt(step)), label=lbl)
end

xlabel(L"x")
ylabel(L"|\psi_n(x)|^2")
title("Infinite well wave function densities in QO Julia")
legend(loc="upper right", bbox_to_anchor=(.95, .8))

gcf()

# println("Energies QO Julia:")
# for n in 1:3
#     println(E[n])
# end
