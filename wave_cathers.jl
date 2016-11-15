function wave_cathers(f1, f2, xf, tf, dx, dt, a, e, h, g)

	J = Int(xf/dx)
	N = Int(tf/dt)
	sol = Array(Float64, 2*J-4, N)
	M1 = spzeros(Float64, J-2, J-2)
	M2 = spzeros(Float64, J-2, J-2)
	M3 = spzeros(Float64, J-2, J-2)

	x = Float64[]
	for i = 0:J
      push!(x, i*dx)
    end
    for i = 1:J-2
		sol[i,1] = f1(x[i+1])
		sol[i+J-2,1] = f2(x[i+1])
    end
    m = Float64[]
    for i = 1:J-1
    	push!(m, (x[i+1] - x[i])/(6*dt))
    end

	for i = 2:J-2
    	M1[i-1,i] = m[i]
    	M1[i,i-1] = m[i]
    	M2[i,i-1] = -h*g/4
    	M2[i-1,i] = h*g/4
    	M3[i,i-1] = -1/4
    	M3[i-1,i] = 1/4
	end
	for i = 1:J-2
		M1[i, i] = 2*m[i] + 2*m[i+1]
	end
	A = [M1 M2; M3 M1]
	B = [M1 -M2; -M3 M1]

	for i = 2:N
		sol[:,i] = A\(B*sol[:,i-1])
	end
	return sol
end