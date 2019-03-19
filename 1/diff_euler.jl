using PyPlot
pygui(false)

m = 1
ω = 1

x0 = 1
p0 = 1

Δt = .1
t_end = 10
steps = t_end / Δt

x = []
p = []
push!(x, x0)
push!(p, p0)

for i in 1:steps-1
    p_last = p[size(p)[1]]
    x_last = x[size(x)[1]]
    push!(x, x_last + p_last / m * Δt)
    p_last = p[size(p)[1]]
    x_last = x[size(x)[1]]
    push!(p, p_last -m * ω^2 * x_last * Δt)
end

figure(figsize=(6, 5))
clf()
suptitle("Harmonischer Oszillator Euler")
t = [Δt:Δt:t_end;]

k = ω^2*m
E = @. 1/2*k*x^2 + p^2/2m

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
