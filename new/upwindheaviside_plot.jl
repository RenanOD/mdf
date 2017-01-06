using Plots
pyplot()

function upwindheaviside(f, dx, dt, U, x0, xf)

  Unum = dx/dt
  s = U/Unum
  x = linspace(x0, xf, (xf-x0)/dx + 1)
  J = length(x)
  N = J

  tf = (N-1)*dt

  sol = Array(Float64, J, N)
  sol[:,1] = f(x)

  A = speye(N, J)

  A[1,1] = 1 - s

  for i in 2:N
    A[i,i] = 1 - s
    A[i,i-1] = s
  end

  for i in 2:N
    sol[:,i] = A*sol[:,i-1]
  end
  sol[1,1] += 0*s
  for i = 2:N
    sol[1,i] += 0*s
  end
  return sol
end

function upwindheaviside_plot()

  function g(x)
    g = Float64[]
    for i in x
      push!(g, i <= 1 ? 0 : 1)
    end
    return g
  end

  dx = 0.00625
  dt = 0.005
  U = 1
  (x0, xf) = (0, 8)
  x = linspace(x0, xf, (xf-x0)/dx + 1)
  sol = upwindheaviside(g, dx, dt, U, x0, xf)

  P = plot(x, sol[:, 1], label="t=0 s", title="f(x)= x <= 1 ? 0 : 1, dx=$dx, dt=$dt, U=$U", ylim=(0,2))

  for i in 1:3
    n = i*div(xf-x0,4*dx)
    t = round(n*dt, 2)
    plot!(x, sol[:, n], label="t=$t s")
  end
  P
end