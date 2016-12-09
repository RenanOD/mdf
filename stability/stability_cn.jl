using Plots
pyplot()

  g = 1.0
  dx = .05
  dt = .0025
  xf = 9.0
  tf = 6.0
  h = 1.0

  J = Int(xf/dx)
  N = Int(tf/dt)
  P = Float64[]
  K = linspace(1, 100, 10)

  for k in K
    M1 = eye(J-2)
    M2 = zeros(Float64, J-2, J-2)
    M3 = zeros(Float64, J-2, J-2)

    for j = 1:J-3
     M2[j,j+1] = k*g*h/2
     M2[j+1,j] = -k*g*h/2
     M3[j,j+1] = k/2
     M3[j+1,j] = -k/2
    end
      
    A = [M1 M2; M3 M1]
    B = [M1 -M2; -M3 M1]

    C = B*inv(A)
    push!(P, maximum(svdvals(C)))
  end

p=plot(K, P, ylim=(0,2), xlabel="k", ylabel="|||C|||", title="Crank-Nicolson")
