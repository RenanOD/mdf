using Plots
pyplot(size=(800,800*3/4))

include("mdf.jl")

n = 7
x = linspace(0, 1, n)
xx = linspace(0, 1, 10*n)
u(x) = 120*4*x.*(1-x) + 26
ϕ = u(x)

m = 3
ɛ = 2e-3
tf = 45
ts = 3
linest = [:dash, :dot, :dashdot]

M = heat_exp(ϕ, 0, tf, x[n], ɛ, m*ts)
E = exact_heat(u, ϕ, 0, tf, x[n], ɛ, 3)

plot(u, 0, 1, c=:black, ls=:solid, label="t=0s")
plot!(xx, E[:,2], c=:blue, ls=:solid, label="t=15s")
plot!(xx, E[:,3], c=:blue, ls=:solid, label="t=15s")
plot!(xx, E[:,4], c=:blue, ls=:solid, label="t=45s")

for i = 1:m
  plot!(x, M[:,i*ts], m=(3,:dot), c=:black, ls=linest[i],
    label="t = $(tf*i/m)")
end

yaxis!("Temperatura")
xaxis!("Posição")
