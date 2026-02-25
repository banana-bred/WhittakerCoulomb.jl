"ν = Z * sqrt(μred)/sqrt(-2E)"
function ν(E; Z :: Real = 1.0, μred :: Real = 1.0)
  @assert E<0 "ν(E) requries E < 0"
  return (Z*sqrt(μred))/sqrt(-2E)
end

@inline function _prefactor(ν::Real, l::Int)
  return ν^0.5 / sqrt(gamma(ν+l+1)*gamma(ν-l) + 0im)
end

"The function A(ε,l) = Π_{n=0}^l (1+n²ε) from Seaton (2002)"
@inline _seatonA(ε,l) = prod(1+n^2*ε for n in 0:l)

"The function B(ε,l) (113) from Seaton (2002)"
@inline _seatonB(ε,l) = ε <=0 ? A(ε,l) : A(ε,l) / ( 1-exp(-2π/sqrt(ε)) )

"The function H(ε,l) = ε≤0 ? 0 : 1/( exp(2π/k) - 1 ) from Seaton (2002); k=sqrt(ε)"
@inline _seatonH(ε) = ε≤0 ? 0 : 1/( exp(2π/sqrt(ε)) - 1 )
