#######################################################################
# Good ol whittaker functions and their derivatives
#######################################################################
@inline function _do_use_ode_W(κ, μ, z)
  a = μ - κ + 0.5
  b = 2μ + 1
  return isapprox(b, round(b); atol=1e-12) && abs(a) > 8 && z < 10
end

@inline function _W_U(κ,μ,z)
  return exp(-z/2) * z^(μ + 1/2) * U(μ - κ + 1/2, 2μ+1, z)
end

"""
    _W_all(κ,μ,z; z0 = nothing, Nasymp :: Int = 40, hmax = 0.05)

Returns W, its derivative, and its log-derivative

"""
function _W_all(κ, μ, z; z0=nothing, Nasymp :: Int = 40, hmax = 0.05)
  _do_use_ode_W(κ, μ, z) && return _W_solve_ODE(κ, μ, z; z0=z0, Nasymp=Nasymp, hmax=hmax)

  a = μ - κ + 0.5
  W0 = _W_U(κ, μ, z)
  W1 = _W_U(κ - 0.5, μ + 0.5, z)
  dWdz = (-0.5 + (μ + 0.5)/z) * W0 - a*W1 / sqrt(z)
  dlogWdz = dWdz / W0
  return W0, dWdz, dlogWdz
end

"""
    W(κ, μ, z)

Standard Whittaker W_{κ,μ}(z)
"""
function W(κ, μ, z; kwargs...)
  W0, _, _ = _W_all(κ, μ, z; kwargs...)
  return W0
end

"""
    dW(κ, μ, z)

Derivative of the standard Whittaker W_{κ,μ}(z)
"""
function dW(κ, μ, z; kwargs...)
  _, dWdz, _ = _W_all(κ, μ, z; kwargs...)
  return dWdz
end

"""
    dlogW(κ, μ, z)

Log-derivative of the standard Whittaker W_{κ,μ}(z)
"""
function dlogW(κ, μ, z; kwargs...)
  _, _, dlogWdz = _W_all(κ, μ, z; kwargs...)
  return dlogWdz
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
  return dM(κ, μ, z) / M(κ, μ, z)
end
