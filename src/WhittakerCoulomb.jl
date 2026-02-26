module WhittakerCoulomb

  using HypergeometricFunctions: _₁F₁ as F11, U
  using CoulombFunctions: coulombs
  using SpecialFunctions: gamma, logabsgamma

  # export ν
  export M, dM, dlogM
  export W, dW, dlogW
  export coulM,    dcoulM, dlogcoulM
  export coulW,    dcoulW, dlogcoulW
  export coulM_EN, dcoulM_EN
  export coulW_EN, dcoulW_EN
  export seatony
  export seatonf, dseatonf
  export seatonh, dseatonh


  include("helpers.jl")
  include("whittaker.jl")
  include("coulomb.jl")

end # module WhittakerCoulomb
