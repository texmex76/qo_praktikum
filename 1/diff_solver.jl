using DifferentialEquations
using PyPlot

dt = .1
tend = 10.0
y0 = Float64[1,1]
m = 1.0
ω = 1.0

tvec=[0:dt:tend;]


function f(dy,y,p,t)
    dy[1] = y[2]/m
    dy[2] = -m*ω^2*y[1]
end

prob = ODEProblem(f,y0,(tvec[1],tvec[end]))

clf()
sol = solve(prob,Tsit5(),saveat=tvec)
t = sol.t
y = sol.u
x = getindex.(y,1)
p = getindex.(y,2)

k = ω^2*m
E = @. 1/2*k*x^2 + p^2/2m

figure(figsize=(6, 5))
clf()
suptitle("Harmonischer Oszillator")
subplots_adjust(top=0.85)

subplot(221)
plot(t, x)
xlabel(L"t")
ylabel(L"x(t)")

subplot(222)
plot(t, p, "C1")
xlabel(L"t")
ylabel(L"p(t)")

subplot(223)
plot(x, p, "C2")
xlabel(L"x(t)")
ylabel(L"p(t)")

subplot(224)
plot(t, E, "C3")
xlabel(L"x(t)")
ylabel(L"E(t)")

tight_layout(rect=[0, 0, 1, .95])
gcf()
