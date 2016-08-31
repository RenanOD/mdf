using Plots
pyplot(size=(800,800*3/4))

include("mdf.jl")

n = 7
x = linspace(0, 1, n)
u(x) = 120*4*x.*(1-x) + 26
ϕ = u(x)

m = 3
ɛ = 2e-3
tf = 45
ts = 3
linest = [:dash, :dot, :dashdot]

M = heat_imp(ϕ, 0, tf, x[n], ɛ, m*ts)

plot(u, 0, 1, c=:black, ls=:solid, label="t=0s")

for i = 1:m
  plot!(x, M[:,i*ts], m=(3,:dot), c=:black, ls=linest[i],
    label="t = $(tf*i/m)")
end

yaxis!("Temperatura")
xaxis!("Posição")

png("imp")