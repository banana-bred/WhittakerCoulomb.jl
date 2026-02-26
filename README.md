Simple procedures to calculate the Whittaker functions, the Whittaker Coulomb functions, the energy-normalized Whittaker Coulomb functions, and their derivatives.
WIP

### Whittaker functions
The Whittaker differential equation is 

$$
\frac{d^2w}{dz^2} + \left( -\frac{1}{4} + \frac{\kappa}{z} + \frac{1/4 - \mu^2}{z^2} \right) w = 0.
$$

Two solutions to the above are the Whittaker functions $W_{\kappa,\mu}(z)$ and $M_{\kapa,\mu}(z), can be defined in terms of the confluent hypergeometric function $_1F_1(a;b;z)$ and the Tricomi confluent hypergeometric function $U(a,b,z)$

$$
W_{\kappa,\mu}(z) = e^{-\frac{z}{2}}  z^{\mu+\frac{1}{2}}  U\left( \mu - \kappa + \frac{1}{2}; 1 + 2\mu; z \right)
$$
$$
M_{\kappa,\mu}(z) = e^{-\frac{z}{2}}  z^{\mu+\frac{1}{2}}  {}_1F_1\left( \mu - \kappa + \frac{1}{2}; 1 + 2\mu; z \right)
$$

### Whittaker Coulomb functions
These are just the regular Whittaker functions above with a specific parameterization in terms of the orbital angular momentum $l$, the effective quantum number $ν = 1/\sqrt{-2E}$ in atomic units, and a distance $r$ (assuming a target charge of 1 and a reduced mass of 1).
- $\kappa\to ν$
- $\mu\to l+1/2$
- $z\to 2r/ν$

### Energy-normalized Whittaker functions
Energy normalization is a common choice for functions like these, especially at positive energies

$$
\int\limits_0^\infty dzW_{\kappa_1,\mu}(z) W_{\kappa_2,\mu}(z) = \delta(\kappa_1 - \kappa_2)
$$

which results in the multiplicative prefactor $\sqrt{ \nu /\Gamma(\nu+l+1)\Gamma(\nu-l) }$ to the Whittaker Coulomb functions
