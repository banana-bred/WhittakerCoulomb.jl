#######################################################################
# Routines for calculating Whittaker functions via their ODE
#######################################################################

"Whittaker ODE: w'' + q(z) w = 0'"
@inline function _q_whittakerODE(κ, μ, z)
  aquart = one(z)/4
  return -aquart + κ/z + (aquart - μ*μ)/(z*z)
end

"RHS of the ODE with y(z) = dlogW = W'/W at z:  y'(z) = -y(z)² - q(z)"
@inline function _RHS_whittakerODE(κ :: Real, μ :: Real, z :: Real, y :: Real)
  @show κ, μ, z, y
  return -y*y - _q_whittakerODE(κ, μ, z)
end

@inline function _choose_z0(κ, μ, z; z0=nothing)
  z0 === nothing || return z > z0 ? z : z0
  return max(z + 50.0, 60.0, 20*abs(κ) + 40.0,  20*abs(μ) + 40.0)
end

"""
    _asymptotic_initializer(κ, μ,  z0; Nasymp :: Int = 60)

Initialize the asymptotic value of u and u' at large z0 where
W ~ e^{-z/2} z^κ u(z) -> u(z) = Σ c_n  z^{-n}
"""
function _asymptotic_initializer(κ, μ,  z0; Nasymp :: Int = 60)
  T = promote_type(typeof(κ), typeof(μ), typeof(z0))
  κ, μ, z0 = T.( (κ, μ, z0) )

  a =  μ - κ + one(T)/2
  c = -μ - κ + one(T)/2

  term = one(T)
  Sn   = term
  Snp1 = zero(T)
  invz = one(T)/z0
  for n in 0:(Nasymp-2)
    term *= -(a+n)*(c+n) * invz / (n+1)
    Sn   += term
    Snp1 += (n+1)*term
  end

  u0  = Sn
  up0 = -Snp1/z0

  return u0, up0

end

function _tricomi_initializer(κ, μ, z0)
  a= μ- κ  + 0.5
  W0 = _W_U(κ, μ, z0)
  W1 = _W_U(κ, μ, z0)

  dW0 = (-0.5 + (μ + 0.5)/z0) * W0 - a * W1 / sqrt(z0)

  scale = exp(z0/2 - κ*log(z0))
  u0 = scale*W0
  up0 = scale * (dW0 + (0.5 - κ/z0) * W0)

  return u0, up0, W0,dW0

end

"""
    _scaled_whittaker_W_RHS!(dU,U,p,z)

RHS of the ODE for the scaled Whittaker e^{z/2} z^{-κ} W_{κ,μ}(z).
U  = [U,U']
dU = [U',U'']
"""
function _scaled_whittaker_W_RHS!(dU, U, p, z)
  κ, μ = p
  u, up = U
  dU[1] = up
  dU[2]= -2*(κ/z - 1/2)*up - ( (κ-1/2)^2 -μ^2 ) / (z*z) * u
  return nothing
end

"""
    _W_solve_ODE(κ, μ, z; z0=nothing, Nasymp = 60, alg=Rodas5P(), reltol=1e-12, abstol=1e-12)

Returns W, dWdz, dlogW/dz by solving the Whittaker ODE
"""
function _W_solve_ODE(κ, μ, z; z0=nothing, Nasymp = 60, alg=Rodas5P(), reltol=1e-12, abstol=1e-12)
  @assert z>0

  T = promote_type(typeof(κ), typeof(μ), typeof(z))
  κ, μ, z = T.( (κ, μ, z) )

  # -- asymptotic starting point
  z0 = T(_choose_z0(κ, μ, z; z0))
  @assert z0 > z "Asymptotic limit ($(z0)) cannot be smaller than the target value ($(z))"

  # -- asymptotic solutions
  u0, up0 = _asymptotic_initializer(κ, μ, z0; Nasymp=Nasymp)

  # -- solve down to z
  u0up0 = T[u0, up0]
  prob = ODEProblem(_scaled_whittaker_W_RHS!, u0up0, (z0, z), (κ, μ))

  sol = solve(prob, alg; reltol=reltol, abstol=abstol)
  u, up = sol.u[end]

  # -- log-derivative
  dlogW = up/u + κ/z - T(1/2)

  # -- build the unscaled solution W_{κ,μ}(z)
  W = exp(-z/2) * z^κ * u

  # -- derivative
  dW = dlogW * W

  # -- extra goodies
  logabsW = -z/2 + κ*log(z) + log(abs(u))
  phaseW = u/abs(u)

  return W, dW, dlogW, logabsW, phaseW

end
