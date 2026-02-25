"""
    W(κ, μ, z)

Standard Whittaker W_{κ,μ}(z) expressed via Tricomi U:

    W_{κ,μ}(z) = exp(-z/2) * z^(μ+1/2) * U(μ-κ+1/2, 2μ+1, z)
"""
function W(κ, μ, z)
  return exp(-z/2) * z^(μ + 1/2) * U(μ - κ + 1/2, 2μ + 1, z)
end

"""
    dW(κ, μ, z)

Derivative of the standard Whittaker W_{κ,μ}(z)

    (d/dz) W_{κ,μ}(z) = (1/2 + (μ+1/2)/z) W_{κ,μ}(z) - W_{κ+1/2, μ+1/2}(z) / √z
"""
function dW(κ, μ, z)
  W0 = W(κ,     μ,     z)
  W1 = W(κ+1/2, μ+1/2, z)
  return (0.5 + (μ + 0.5)/z) * W0 - W1/sqrt(z)
end

"""
    dlogW(κ, μ, z)

Log-derivative of the standard Whittaker W_{κ,μ}(z)

  (d/dz) log(W_{κ,μ}(z)) = (1/2 + (μ+1/2)/z) - W_{κ+1/2, μ+1/2}(z) / (W_{κ,μ}(z)*√z)
"""
function dlogW(κ, μ, z)
  W0 = W(κ,     μ,     z)
  W1 = W(κ+1/2, μ+1/2, z)
  return (0.5 + (μ + 0.5)/z) - W1/(W0*sqrt(z))
end

####################################################################
"""
    M(κ, μ, z)

Standard Whittaker M_{κ,μ}(z) expressed via Kummer 1F1:

    M_{κ,μ}(z) = exp(-z/2) * z^(μ+1/2) * 1F1(μ-κ+1/2; 2μ+1; z)
"""
function M(κ, μ, z)
  return exp(-z/2) * z^(μ + 1/2) * F11(μ - κ + 1/2, 2μ + 1, z)
end

"""
    dM(κ, μ, z)

Derivative of the standard Whittaker M_{κ,μ}(z)

    (d/dz) M_{κ,μ}(z) = (-1/2 + (μ+1/2)/z) M_{κ,μ}(z) +  (a/b) * M_{κ-1/2, μ+1/2}(z) / √z
"""
function dM(κ, μ, z)
  a = μ - κ + 0.5
  b = 2μ + 1
  M0 = M(κ,       μ,       z)
  M1 = M(κ - 0.5, μ + 0.5, z)
  return (-0.5 + (μ + 0.5)/z) * M0 + (a/b) * M1 / sqrt(z)
end

"""
    dlogM(κ, μ, z)

Log-derivative of the standard Whittaker M_{κ,μ}(z)
"""
function dlogM(κ, μ, z)
  a = μ - κ + 0.5
  b = 2μ + 1
  M0 = M(κ,       μ,       z)
  M1 = M(κ - 0.5, μ + 0.5, z)
  return (-0.5 + (μ + 0.5)/z) + (a/b) * M1 / (M0*sqrt(z))
end
