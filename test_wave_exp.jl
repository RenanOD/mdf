include("wave_exp.jl")

g = 1
dx = .05
dt = .005
xf = 10
tf = 8
a = 1
e = 0.1
h = 1
f1(x) = exp(-(x - a)^2/e)
f2(x) = exp(-(x - a)^2/e)

sol = wave_exp(f1, f2, xf, tf, dx, dt, a, e, h, g)

J = Int(xf/dx)
N = Int(tf/dt)

est = 15
passo = floor(N/(est-1))
aux = zeros(J-2, est)
j = 1
k = 1
aux[:, k] = sol[J-1:2*J-4,1]

for i = 2:N
  j = j + 1
  if j == passo
    k = k + 1
    aux[:, k] = sol[J-1:2*J-4, i]
    j = 0
    if k >= est
      break
    end
  end
end

x = Float64[]
for i = 1:J
  push!(x, i*dx)
end

a = plot(x, aux[:, 1])

for i = 2:est
  plot!(x, aux[:, i])
end
a













