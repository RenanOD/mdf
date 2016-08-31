function dirichlet(f, cc1, cc2, x0, xn, k1, k2, N)
	h = (xn - x0)/N
	A = spzeros(N-1, N-1)
	F = zeros(N-1, 1)

	F[1] = k1*cc1	; F[N-1] = k1*cc2
	A[1,2] = -k1; A[N-1, N-2] = -k1

	for i = 1:N-1
		A[i, i] = 2*k1 + h^2*k2
	end
	for i = 2 : N-2
		A[i, i + 1] = -k1
		A[i, i - 1] = -k1
		F[i] = h^2*f(x0+i*h)
	end
	U = A \ F
	return U
end

function neumann(f, cc1, cc2, x0, xn, k, N)
	h = (xn - x0)/N
	A = spzeros(N+1, N+1)
	F = zeros(N+1, 1)

	F[1] = h^2*f(x0) - 2*h*cc1	; F[N+1] = h^2*f(xn) + 2*h*cc2
	A[1,2] = -2; A[N+1, N] = -2

	for i = 1:N+1
		A[i, i] = 2 + h^2*k
	end
	for i = 2 : N
		A[i, i + 1] = -1
		A[i, i - 1] = -1
		F[i] = h^2*f(x0+(i - 1)*h)
	end
	U = A \ F
	return U
end

function heat_exp(u0, t0, tf, L, k, M)
	d = (tf - t0)/M
	ul = length(u0)
	h = L/(ul - 1)
	A = spzeros(ul, ul)
	U = Array(Float64, ul, M+1)

	A[1,1] = 1; A[ul, ul] = 1

	for i = 2:ul-1
		A[i, i] = 1 - 2k*d/h^2
		A[i, i+1] = k*d/h^2
		A[i, i-1] = k*d/h^2
	end

	U[1:ul, 1] = u0

	for j = 2:M+1
		U[1:ul, j] = A*U[1:ul, j-1]
	end

	return U
end

function heat_imp(u0, t0, tf, L, k, M)
	d = (tf - t0)/M
	ul = length(u0)
	h = L/(ul - 1)
	A = spzeros(ul-2, ul-2)
	b = Array(Float64, ul-2, M)
	U = Array(Float64, ul, M+1)

	A[1,1] = -2*k-h^2/d; A[1, 2] = k;
	A[ul-2, ul-2] = -2*k-h^2/d; A[ul-2, ul-3] = k;

	for i = 2:ul-3
		A[i, i] = -2*k-h^2/d
		A[i, i+1] = k
		A[i, i-1] = k
	end

	U[:, 1] = u0
	U[1,:] = u0[1]
	U[ul,:] = u0[ul]

	for j = 2:M+1
	  b[:,j-1] = -h^2*U[2:ul-1,j-1]/d
	  b[1,j-1] += -k*u0[1]
	  b[ul-2,j-1] += -k*u0[ul]
	  U[2:ul-1,j] = A\b[:,j-1]
	end

	return U
end