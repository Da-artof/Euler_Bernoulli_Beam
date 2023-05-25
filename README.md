# Euler_Bernoulli_Beam
A C program to find the displacement of a cantilevered beam from equilibrium subject to various force profiles.
We find time - stationary solutions for the displacement $w\left(x\right)$ to the Euler-Lagrange equation for a cantilevered beam with a forcing $q\left(x\right)\sin\ \omega t$:

$$\mu\left(x\right)\omega^2w\left(x\right)= -\frac{\partial^2}{\partial x^2}\left(K\left(x\right)\frac{\partial^2}{\partial x^2}\right)w\left(x\right) + q\left(x\right)$$

where $\mu\left(x\right)$ is the mass density and $K\left(x\right)$ is the flexural rigidity of the beam. 

The boundary conditions are:

at $x = 0$: $w\left(x\right) = 0$ and $\frac{\partial w\left(x\right)}{\partial x} = 0$

at $x = L$: $\frac{\partial^2w\left(x\right)}{\partial x^2} = 0$ and $\frac{\partial^3w\left(x\right)}{\partial x^3} = 0$

This is a 4th order ordinary differential equation, as the displacement here is only a function of $x$. We use an explicit finite difference method to form a banded matrix $A$ of coefficients, using the boundary conditions to assign values at the endpoints, then solve $A\boldsymbol{w} = \boldsymbol{q}$. More information on the underlying physics can be found at: https://en.wikipedia.org/wiki/Euler–Bernoulli_beam_theory.
