using NumInt

function u_ex(x, t, g, eps, L; minh = 1e-4, maxn = 11)
  if t == 0
    return g(x)
  elseif x == 0 || x == L
    return g(x)
  end
  
  c(n) = (2/L)*clenshaw_rule(x->(g(x)-(g(L) - g(0))*x/L - g(0))*sin(n*pi*x/L), 0, L)

  h = 1.0
  n = 1
  r = ((g(L) - g(0))/L)*x + g(0)

  while n < maxn && h > minh
    cn = c(n)
    h = exp(-(n*pi/L)^2*eps*t)
    r += h * sin(n*pi*x/L) * cn
    n += 1
  end
  return r
end
