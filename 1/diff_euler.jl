using PyPlot

m = 1
ω = 8

x0 = 1
p0 = 1

Δt = .01
t = 2
steps = t / Δt

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
suptitle("Harmonischer Oszillator")
t = [Δt:Δt:t;]

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

tight_layout()
gcf()
