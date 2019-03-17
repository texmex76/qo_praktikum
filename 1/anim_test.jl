using PyPlot
using PyCall
using Printf
anim = pyimport("matplotlib.animation")
#pygui(true)

fig = figure(1)
clf()
x = [0:.01:2;] * pi
y = sin.(x)

function animate(k)
    k += 1
    clf()
    plot(x, y.*k)
    ylim(-105, 105)
    title("k=" *@sprintf("%03.2f", k))
    xlabel(L"x")
    ylabel(L"k \sin(x)")
    grid()
end

frames = [50:1:100;]
movie = anim.FuncAnimation(fig, animate, frames=frames, repeat=false, interval=20)
movie.save("animation2.mp4")
