#######################################################################
# Routines for calculating Whittaker functions via their ODE
#######################################################################

"Whittaker ODE: w'' + q(z) w = 0'"
@inline function _q_whittakerODE(κ, μ, z)
  aquart = one(z)/4
  return -aquart + κ/z + (aquart - μ*μ)/(z*z)
end

"RHS of the ODE with y(z) = dlogW = W'/W at z:  y'(z) = -y(z)² - q(z)"
@inline function _RHS_whittakerODE(κ :: Real, μ :: Real, z :: Real, y :: Real)
  return -y*y - _q_whittakerODE(κ, μ, z)
end

"""
    _W_asymptotic

Evalutes the asymptotic form of W_{κ,μ}(z) for large z:

    W ~ e^{-z/2} z^κ Σ_{n≥0} c_n z^{-n}
    a =  μ - κ + 1/2
    c = -μ - κ + 1/2
    c₀ = 1
    c_{n+1} = c_n * -(a+n)*(c+n)/((n+1)*z)

Returns:
- logabsW = logabs(|W(z0)|)
- phase = W(z0)/|W(z0)|
- dlogW
"""
function _W_asymptotic(κ, μ, z0; Nasymp::Int = 40)
  @assert z0 > 0
  T = promote_type(typeof(κ), typeof(μ), typeof(z0))
  κ, μ, z0 = T.( (κ, μ, z0) )

  a =  μ - κ + one(T)/2
  c = -μ - κ + one(T)/2

  term = one(T)
  S    = term
  Sn   = zero(T)

  invz = one(T)/z0
  for n in 0:(Nasymp-2)
    term *= -(a+n)*(c+n) * invz / (n+1)
    S += term
    Sn += (n+1)*term
  end

  W = exp(-z0/2) * z0^κ * S
  absW = abs(W)
  logabsW = log(absW)
  phase = W/absW
  dlogW = -T(0.5) + κ/z0 - (Sn/S)/z0

  return logabsW, phase, dlogW

end

"""
    _W_solve_ODE(κ :: Real, μ :: Real, ζ :: Real; z0=nothing, Nasymp :: Int = 40, hmax :: Real = 0.05)

Computes W, dW, dW/W for real z>0 using the asymptotic expansion of W for large z
and numerical integration of the Whittaker ODE towards the target z

Returns:
- W = W_{κ,μ}(z),
- dW = dW_{κ,μ}(z)/dz,
- dlogW = dW/W at z
"""
function _W_solve_ODE(κ :: Real, μ :: Real, z :: Real; z0=nothing, Nasymp::Int=40, hmax :: Real=0.05)
  @assert z > 0 "z must be positive"

  # -- homogenize types
  T = promote_type(typeof(κ), typeof(μ), typeof(z), typeof(hmax))
  κ, μ, z, hmax = T.( (κ, μ, z, hmax) )

  # -- get z0; if not supplied, pick one that's probably big enough
  z0 = z0 === nothing ? max(z + T(50), T(60), abs(κ) + abs(μ) + T(50)) : T(z0)

  # TODO: replace this with a simple asymptotic evaluation that returns
  @assert z0 > z "cant step down"

  # -- asymptotic solution at z=z0
  # log|W|, phase, dlogW = ..
  L, phase, y = _W_asymptotic(κ, μ, z0; Nasymp=Nasymp)

  # -- step inwards down to z
  span = z0 - z
  nsteps = max(50, Int(ceil(span / hmax)))
  h = -span / nsteps

  zi = z0
  for _ in 1:nsteps
    # -- RK4 for y and log|W|
    k1y = h * _RHS_whittakerODE(κ, μ, zi, y) # <-- y
    k1L = h * y                # <-- log|W|

    z2  = zi + h/2
    y2  = y  + k1y/2
    k2y = h * _RHS_whittakerODE(κ, μ, z2, y2)
    k2L = h * y2

    y3  = y  + k2y/2
    k3y = h * _RHS_whittakerODE(κ, μ, z2, y3)
    k3L = h * y3

    z4  = zi + h
    y4  = y  + k3y
    k4y = h * _RHS_whittakerODE(κ, μ, z4, y4)
    k4L = h * y4

    y  += (k1y + 2k2y + 2k3y + k4y) / 6
    L  += (k1L + 2k2L + 2k3L + k4L) / 6
    zi += h
  end

  W = phase * exp(L)
  dWdz = y*W
  return W, dWdz, y
  #      W, dWdz, dlogWdz

end
