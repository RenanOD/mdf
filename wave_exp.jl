function wave_exp(f1, f2, xf, tf, dx, dt, a, e, H)

	g = 1.0
	J = Int(xf/dx)
	N = Int(tf/dt)
	sol = Array(Float64, 2*J-4, N)
	M1 = spzeros(Float64, J-2, J-2)
	M2 = spzeros(Float64, J-2, J-2)
	M3 = spzeros(Float64, J-2, J-2)
	k1 = dt*g*H/(2*dx)
	k2 = dt/(2*dx)

	for i = 1:J-2
		sol[i,1] = f1((i+1)*dx)
		M1[i,i] = 1
	end
	for i = 1:J-3
		M2[i,i+1] = k1
		M2[i+1,i] = -k1
		M3[i,i+1] = k2
		M3[i+1,i] = -k2
	end
	for i = J-1:2*J-4	
		sol[i,1] = f2((i+1)*dx)
	end

	B = [M1 M2; M3 M1]

	for i = 2:N
		sol[:,i] = B*sol[:,i-1]
	end
	return sol
end



