using NumInt

function u_ex(x, t, u, k, L)
	c(n) = (2/L)*clenshaw_rule(x->(u(x)-(u(L) - u(0))*x/L - u(0))*sin(n*pi*x/L), 0, L)

	r = ((u(L) - u(0))/L)*x + u(0)
	for n = 1:11
		r += c(n)*exp(-n^2*pi^2*k*t/L^2)*sin(n*pi*x/L)
	end
	return r
end

function exact_heat(u, u0, t0, tf, L, k, M)

	t = linspace(t0, tf, M+1)
	ul = length(u0)
	xx = linspace(0, L, 10*ul)
	U = Array(Float64, 10*ul, M+1)
	for j = 1:M+1
		for i = 1:10*ul
			U[i,j] = u_ex(xx[i], t[j], u, k, L)
		end
	end
	return U
end

