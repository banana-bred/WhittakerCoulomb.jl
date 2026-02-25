module WhittakerCoulomb

  using HypergeometricFunctions: _₁F₁ as F11, U
  using CoulombFunctions: coulombs
  using SpecialFunctions: gamma

  export ν
  export W, dW, dlogW
  export M, dM, dlogM
  export coulW,    dcoulW
  export coulW_EN, dcoulW_EN
  export coulM,    dcoulM
  export coulM_EN, dcoulM_EN
  export seatony
  export seatonf, dseatonf
  export seatonh, dseatonh


  include("helpers.jl")
  include("whittaker.jl")
  include("coulomb.jl")

end # module WhittakerCoulomb
