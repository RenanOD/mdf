include("wave_cn.jl")

g = 1
dx = .05
dt = .005
xf = 9
tf = 6
a = 1
e = 0.1
h = 1
f1(x) = exp(-(x - a)^2/e)
f2(x) = exp(-(x - a)^2/e)

sol = wave_cn(f1, f2, xf, tf, dx, dt, a, e, h, g)

J = Int(xf/dx)
N = Int(tf/dt)

est = 180
passo = floor(N/(est))
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
for i = 1:(J-2)
  push!(x, i*dx)
end

tt = 0
anim = @animate for i = 1:est
t = (i-1)*(N/est)*dt
if t%.5 < dt
  tt = t
  plot(x, aux[:, i], line=(3), label="$t s", title = "Crank-Nicolson (dt = $dt s)", xlabel="Posicao", ylabel="Altura", ylim=(-.7,1.4))
else
  plot(x, aux[:, i], line=(3), label="$tt s", title = "Crank-Nicolson (dt = $dt s)", xlabel="Posicao", ylabel="Altura", ylim=(-.7,1.4))
end
end

gif(anim, "cn.gif", fps=30)
