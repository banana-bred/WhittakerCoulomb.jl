module WhittakerCoulomb

  using HypergeometricFunctions: _₁F₁ as F11, U
  using SpecialFunctions: gamma

  export ν
  export W, dW, dlogW
  export M, dM, dlogM
  export coulW,    dcoulW
  export coulW_EN, dcoulW_EN
  export coulM,    dcoulM
  export coulM_EN, dcoulM_EN

  ####################################################################
  # -- Helper  Functions

  "ν = Z * sqrt(μred)/sqrt(-2E)"
  function ν(E; Z :: Real = 1.0, μred :: Real = 1.0)
    @assert E<0 "ν(E) requries E < 0"
    return (Z*sqrt(μred))/sqrt(-2E)
  end

  @inline function _prefactor(ν::Real, l::Int)
    return ν^0.5 / sqrt(gamma(ν+l+1)*gamma(ν-l))
  end

  ####################################################################

  # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  # vvvvvvv EXPORT DEFINITIONS vvvvvvv
  # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  ####################################################################
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

  ####################################################################
  """
      coulW(r, ν, l; Z, μred)

  Standard unnormalized Coulomb-Whittaker decaying solution for ε < 0

      W_{ν, l+1/2}(z)

  with

     z = 2κr
     κ = sqrt(-2ε) = Z/ν

  """
  function coulW(r, ν :: Real, l :: Int; Z :: Real = 1.0, μred :: Real = 1.0)
    μ = l + 0.5
    z = 2Z*μred *  r/ν
    return W(ν, μ, z)
  end

  ####################################################################
  """
      dcoulW(r, ν, l; Z, μred)

  Derivative of the standard unnormalized Coulomb-Whittaker decaying solution for ε < 0

      W_{ν, l+1/2}(z)
  """
  function dcoulW(r, ν :: Real, l :: Int; Z :: Real = 1.0, μred :: Real = 1.0)
    μ = l + 0.5
    dzdr = 2Z*μred /ν
    z = dzdr * r
    return dzdr * dW(ν, μ, z)
  end

  ####################################################################
  """
      coulM(r, ν, l; Z, μred)

  Standard unnormalized Coulomb-Whittaker increasing solution for ε < 0

      M_{ν, l+1/2}(z)

  with

     z = 2κr
     κ = sqrt(-2ε) = Z/ν

  """
  function coulM(r, ν :: Real, l :: Int; Z :: Real = 1.0, μred :: Real = 1.0)
    μ = l + 0.5
    z = 2Z*μred * r/ν
    return M(ν, μ, z)
  end

  ####################################################################
  """
      dcoulM(r, ν, l; Z, μred)

  Derivative of the standard unnormalized Coulomb-Whittaker increasing solution for ε < 0

      M_{ν, l+1/2}(z)
  """
  function dcoulM(r, ν::Real, l::Int; Z :: Real = 1.0, μred :: Real = 1.0)
    μ = l + 0.5
    dzdr = 2Z*μred /ν
    z = dzdr * r
    return dzdr * dM(ν, μ, z)
  end

  ####################################################################
  """
      coulW_EN(r, ν, l; Z)

  Energy normalized decaying Coulomb-Whittaker function W(r,ν,l); 2.53 in Aymar, Greene, Luc-Koeing (1996)
  """
  function coulW_EN(r, ν::Real, l::Int; Z :: Real = 1.0, μred :: Real = 1.0)
    pref = _prefactor(ν, l)
    return pref * coulW(r, ν, l; Z = Z, μred = μred)
  end

  """
      dcoulW_EN(r, ν, l; Z)

  Derivative of the energy normalized decaying Coulomb-Whittaker function W(r,ν,l)
  """
  function dcoulW_EN(r, ν::Real, l::Int; Z :: Real = 1.0, μred :: Real = 1.0)
    pref = _prefactor(ν, l)
    return pref * dcoulW(r, ν, l; Z = Z, μred = μred)
  end

  ####################################################################
  """
      coulM_EN(r, ν, l; Z)

  Energy normalized increasing Coulomb-Whittaker function M(r,ν,l)
  """
  function coulM_EN(r, ν::Real, l::Int; Z::Real=1.0, μred :: Real = 1.0)
    pref = _prefactor(ν, l)
    return pref * coulM(r, ν, l; Z = Z, μred = μred)
  end

  """
      dcoulM_EN(r, ν, l; Z)

  Derivative of the energy normalized decaying Coulomb-Whittaker function M(r,ν,l)
  """
  function dcoulM_EN(r, ν::Real, l::Int; Z::Real=1.0, μred :: Real = 1.0)
    pref = _prefacto(ν, l)r
    return pref * dcoulM(r, ν, l; Z = Z, μred = μred)
  end

end # module WhittakerCoulomb
