using PyPlot
# pygui(true)

x = [0:.01:2;] * pi
y = sin.(x)

figure(1)
clf()
plot(x/pi, y)
xlabel(L"x/ \pi")
ylabel(L"\sin(x)")
grid()
gcf()
