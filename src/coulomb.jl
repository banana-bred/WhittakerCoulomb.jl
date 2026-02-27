#########################################
############### HELPERS #################
#########################################
"ν = Z * sqrt(μred)/sqrt(-2E)"
function _ν(E; Z :: Real = 1.0, μred :: Real = 1.0)
  @assert E<0 "ν(E) requries E < 0"
  return (Z*sqrt(μred))/sqrt(-2E)
end

@inline function _prefactor_EN(ν::Real, l::Int)
  logabs, phase = _logprefactor_EN(ν, l)
  return exp(logabs)*phase
end

@inline function _logprefactor_EN(ν :: Real, l::Int)

  @assert ν > 0 "ν must be > 0"
  x1 = ν + l + 1
  x2 = ν - l

  # -- x2=ν-l will take us into the territory of poles in Γ
  if x2 <= 0 && x2 == round(x2)
    throw(DomainError(x2, "Γ(ν-1) has poles at nonpositive integers"))
  end

  logabs1, sgn1 = logabsgamma(x1)
  logabs2, sgn2 = logabsgamma(x2)

  # -- log( |ν / (Γ(x1) * Γ(x2))| ); already asserted ν=|ν| because real and positive
  logabs_dtmp = log(ν) - (logabs1 + logabs2)
  logabs_Enorm = 0.5*logabs_dtmp
  phase = (sgn1*sgn2 > 0) ? 1.0 + 0.0im : im

  return logabs_Enorm, phase

end

"The function A(ε,l) = Π_{n=0}^l (1+n²ε) from Seaton (2002)"
@inline _seatonA(ε,l) = prod(1+n^2*ε for n in 0:l)

"The function B(ε,l) (113) from Seaton (2002)"
@inline _seatonB(ε,l) = ε <=0 ? _seatonA(ε,l) : _seatonA(ε,l) / ( 1-exp(-2π/sqrt(ε)) )

"The function H(ε,l) = ε≤0 ? 0 : 1/( exp(2π/k) - 1 ) from Seaton (2002); k=sqrt(ε)"
@inline _seatonH(ε) = ε≤0 ? 0 : 1/( exp(2π/sqrt(ε)) - 1 )

####################################################################
"""
    σ(η, l::Int)

Coulomb phase shift

    σ_l = (logΓ(l+1 - iη) - logΓ(l+1 + iη)) * i;2

or its analytic continuation
"""
@inline function σl(η, l::Int)
  lp1 = l+1
  return ( loggamma(lp1 - im*η) - loggamma(lp1 + im*η) ) * im/2
end

@inline function ctemp(η, l::Int)
  σ = σl(η, l)
  return exp((π/2) * (im*l - η) - im *σ)
end

####################################################################
"""
    coulW(r, ν, l; Z, μred)

Standard unnormalized Coulomb-Whittaker decaying solution for ε < 0
"""
function coulW(r, ν :: Real, l :: Int; Z :: Real = 1.0, μred :: Real = 1.0, z0=nothing)
  μ = l + 0.5
  z = 2Z*μred * r/ν
  return W(ν, μ, z; z0=z0)
end

####################################################################
"""
    dcoulW(r, ν, l; Z, μred)

Derivative of the standard unnormalized Coulomb-Whittaker decaying solution for ε < 0
"""
function dcoulW(r, ν :: Real, l :: Int; Z :: Real = 1.0, μred :: Real = 1.0, z0=nothing)
  μ = l + 0.5
  dzdr = 2Z*μred/ν
  z = dzdr*r
  return dzdr * dW(ν, μ, z; z0=z0)
end

####################################################################
"""
    dlogcoulW(r, ν, l; Z, μred)

Logarithmic derivative of the standard unnormalized Coulomb-Whittaker decaying solution for ε < 0
"""
function dlogcoulW(r, ν :: Real, l :: Int; Z :: Real = 1.0, μred :: Real = 1.0, z0=nothing)
  μ = l + 0.5
  dzdr = 2Z*μred/ν
  z = dzdr*r
  return dzdr * dlogW(ν, μ, z; z0=z0)
end

####################################################################
"""
    coulM(r, ν, l; Z, μred)

Standard unnormalized Coulomb-Whittaker increasing solution for ε < 0
in terms of Coulomb functions
"""
function coulM(r, ν :: Real, l :: Int; Z :: Real = 1.0, μred :: Real = 1.0, z0=nothing)
  μ = l + 0.5
  z = 2Z*μred*r/ν
  return M(ν, μ, z, z0=z0)
end

####################################################################
"""
    dcoulM(r, ν, l; Z, μred)

Derivative of the standard unnormalized Coulomb-Whittaker increasing solution for ε < 0
in terms of the Coulomb functions
"""
function dcoulM(r, ν :: Real, l :: Int; Z :: Real = 1.0, μred :: Real = 1.0, z0=nothing)
  μ = l + 0.5
  dzdr = 2Z*μred/ν
  z = dzdr * r
  return dzdr * dM(ν, μ, z; z0=z0)
end

####################################################################
"""
    dlogcoulM(r, ν, l; Z, μred)

Logarithmic derivative of the standard unnormalized Coulomb-Whittaker increasing solution for ε < 0
in terms of the Coulomb functions
"""
function dlogcoulM(r, ν :: Real, l :: Int; Z :: Real = 1.0, μred :: Real = 1.0, z0=nothing)
  μ = l + 0.5
  dzdr = 2Z*μred/ν
  z = dzdr * r
  return dzdr * dlogM(ν, μ, z; z0=z0)
end

####################################################################
"""
    coulW_EN(r, ν, l; Z)

Energy normalized decaying Coulomb-Whittaker function W(r,ν,l); 2.53 in Aymar, Greene, Luc-Koeing (1996)
"""
function coulW_EN(r, ν::Real, l::Int; Z :: Real = 1.0, μred :: Real = 1.0, z0=nothing)
  pref = _prefactor_EN(ν, l)
  return pref * coulW(r, ν, l; Z = Z, μred = μred, z0=z0)
end

"""
    dcoulW_EN(r, ν, l; Z)

Derivative of the energy normalized decaying Coulomb-Whittaker function W(r,ν,l)
"""
function dcoulW_EN(r, ν::Real, l::Int; Z :: Real = 1.0, μred :: Real = 1.0, z0=nothing)
  pref = _prefactor_EN(ν, l)
  return pref * dcoulW(r, ν, l; Z = Z, μred = μred, z0=z0)
end

####################################################################
"""
    coulM_EN(r, ν, l; Z)

Energy normalized increasing Coulomb-Whittaker function M(r,ν,l)
"""
function coulM_EN(r, ν::Real, l::Int; Z::Real=1.0, μred :: Real = 1.0, z0=nothing)
  pref = _prefactor_EN(ν, l)
  return pref * coulM(r, ν, l; Z = Z, μred = μred, z0=z0)
end

"""
    dcoulM_EN(r, ν, l; Z)

Derivative of the energy normalized decaying Coulomb-Whittaker function M(r,ν,l)
"""
function dcoulM_EN(r, ν::Real, l::Int; Z::Real=1.0, μred :: Real = 1.0, z0=nothing)
  pref = _prefactor_EN(ν, l)
  return pref * dcoulM(r, ν, l; Z = Z, μred = μred, z0=z0)
end

####################################################################
"""
    seatony(κ, λ, z)

The function y(κ, λ; z) in Seaton (2002) equation (14)

    y(κ,λ;z) = κ^(λ+1/2)/Γ(2λ+1) M_{κ,λ}(z)

where M_{κ,λ}(z) is the standard Whittaker function M
"""
function seatony(κ :: Real, λ :: Int, z :: Real)
  return κ^(λ+0.5)/gamma(2λ+1) * M(κ,λ,z)
end

####################################################################
"""
    seatonf(ε, l, r)

The function f(ε, l, r) in Seaton (2002) equation (22)

    f(ε,l;r) = y(κ=1/sqrt(-ε), λ=l+1/2, z=2r/κ)

Note that Seaton's ε = 2ħE/[m(Z₁Z₂e²)²] (2E in atomic units)
"""
function seatonf(ε, l, r)
  κ = 1/sqrt(-ε)
  z = 2r/κ
  λ = l + 0.5
  return seatony(κ, λ, z)
end

####################################################################
"""
    dseatonf(ε, l, r)

The derivative of the function f(ε, l, r) in Seaton (2002) equation (22)
in terms of contiguous functions

        l*f'(l) =               f(l-1) - [l²/r - 1]  * f(l)    (79)
    (l+1)*f'(l) = [l+1)²/r - 1]*f(l)   - [1+(l+1)²ε] * f(l+1)  (80)

Note that Seaton's ε = 2ħE/[m(Z₁Z₂e²)²] (2E in atomic units)
"""
function dseatonf(ε, l, r)
  @assert l >=0 "l ($(l)) must be non-negative"
  f(l) = seatonf(ε, l, r)
  if l == 0
    return ( ((l+1)^2/r - 1)*f(l) - (1+(l+1)^2*ε)*f(l+1) ) / (l+1)
  end
  return (f(l-1) - (l^2/r - 1)*f(l)) / l
end

"""
    seatonh(ε, l, r)

The function h(ε, l, r) in Seaton (2002) equation (49)

    h(ε,l;r) = Γ(l+1-κ)/(πκ^l) * W_{κ,l+1/2}(2r/κ)
             + [iH(ε) + cot(πκ)] A(ε,l) f(ε,l;r)   (r>0)
"""
function seatonh(ε, l, r)
  @assert r>0 "h(ε,l,r) is only defined here for r > 0"
  κ=1/sqrt(-ε)
  λ=l+0.5
  z=2r/κ
  return gamma(l+1-κ)/(π*κ^l) * W(κ,λ,z) +
  (im*_seatonH(ε) + cot(π*κ)) * _seatonA(ε,l)*seatonf(ε,l,r)
end

####################################################################
"""
    dseatonh(ε, l, r)

The derivative of the function h(ε, l, r) in Seaton (2002) equation (49)
in terms of contiguous functions

        l*h'(l) =        (l+l²ε)*h(l-1) - (l²/r-1)*h(l)    (81)
    (l+1)*f'(l) = ((l+1)²/r - 1)*h(l)   -          h(l+1)  (82)

Note that Seaton's ε = 2ħE/[m(Z₁Z₂e²)²] (2E in atomic units)
"""
function dseatonh(ε, l, r)
  @assert l >=0 "l ($(l)) must be non-negative"
  h(l) = seatonh(ε, l, r)
  if l == 0
    return ( ((l+1)^2/r - 1)*h(l) - h(l+1) ) / (l+1)
  end
  return ((1+l^2*ε)*h(l-1) - (l^2/r-1)*h(l)) / l
end

# ####################################################################
# """
#     coulM_alt(r, ν, l; Z)

# Alternate normalization of the increasing Coulomb-Whittaker function M(r,ν,l) where

#     2M(r,ν,l; Z) = f/sin(β) + h/cos(β)

# where β=π(ν-l)
# """
# function coulM_alt(r, E::Real, l::Int; Z::Real=1.0, μred :: Real = 1.0)
#   β = π*(v(E,Z=Z,μred=μred)-l)
#   ε = 2E
#   return seatonf(ε,l,r)/sin(β) + seatonh(ε,l,r)/cos(β)
# end

# ####################################################################
# """
#     dcoulM_alt(r, ν, l; Z)

# Derivative of the alternate normalization of the increasing Coulomb-Whittaker function M(r,ν,l) where

#     2M(r,ν,l; Z) = f/sin(β) + h/cos(β)

# where β=π(ν-l)
# """
# function dcoulM_alt(r, E::Real, l::Int; Z::Real=1.0, μred :: Real = 1.0)
#   β = π*(v(E,Z=Z,μred=μred)-l)
#   ε = 2E
#   return dseatonf(ε,l,r)/sin(β) + dseatonh(ε,l,r)/cos(β)
# end

# ####################################################################
# """
#     coulW_alt(r, ν, l; Z)

# Alternate normalization of the decreasing Coulomb-Whittaker function W(r,ν,l) where

#     2W(r,ν,l; Z) = - f/cos(β) + h/sin(β)

# where β=π(ν-l)
# """
# function coulM_alt(r, E::Real, l::Int; Z::Real=1.0, μred :: Real = 1.0)
#   β = π*(v(E,Z=Z,μred=μred)-l)
#   ε = 2E
#   return -seatonf(ε,l,r)/cos(β) + seatonh(ε,l,r)/sin(β)
# end

# ####################################################################
# """
#     coulW_alt(r, ν, l; Z)

# Derivative of the alternately normlzed decreasing Coulomb-Whittaker function W(r,ν,l) where

#     2W(r,ν,l; Z) = - f/cos(β) + h/sin(β)

# where β=π(ν-l)
# """
# function dcoulM_alt(r, E::Real, l::Int; Z::Real=1.0, μred :: Real = 1.0)
#   β = π*(v(E,Z=Z,μred=μred)-l)
#   ε = 2E
#   return -dseatonf(ε,l,r)/cos(β) + dseatonh(ε,l,r)/sin(β)
# end
