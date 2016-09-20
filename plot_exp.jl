using Plots
pyplot(size=(800,800*3/4))

include("mdf.jl")

n = 7
x = linspace(0, 1, n)
xx = linspace(0, 1, 10*n)
g(x) = 120*4*x.*(1-x) + 26
ϕ = g(x)

m = 3
ɛ = 2e-3
tf = 45
ts = 3
linest = [:dash, :dot, :dashdot]

M = heat_exp(ϕ, 0, tf, x[n], ɛ, m*ts)

plot()

for i = 0:m
  F = Float64[]
  for j in xx
    push!(F, u_ex(j, i*tf/m, g, ɛ, x[n]))
  end
  t = i*tf/m
  plot!(xx, F, c=:blue, ls=:solid, label="t=$t s")
end

for i = 1:m
  plot!(x, M[:,i*ts], m=(3,:dot), c=:black, ls=linest[i],
    label="t = $(tf*i/m)")
end

yaxis!("Temperatura")
xaxis!("Posição")
