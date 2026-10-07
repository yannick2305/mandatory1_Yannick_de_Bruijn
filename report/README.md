# Mandatory assignment 1

Yannick de Bruijn

Answers to 1.2.3 and 1.2.4. Notation: $k_x = m_x\pi$, $k_y = m_y\pi$, $h = 1/N$, $C = c\Delta t/h$.

## 1.2.3 Exact solution

Claim: $u(t,x,y) = e^{\imath(k_x x + k_y y - \omega t)}$ solves $u_{tt} = c^2 \nabla^2 u$.

Differentiating the exponential:

```math
u_{tt} = (-\imath\omega)^2 u = -\omega^2 u, \qquad
u_{xx} = (\imath k_x)^2 u = -k_x^2 u, \qquad
u_{yy} = (\imath k_y)^2 u = -k_y^2 u
\quad\Rightarrow\quad
\nabla^2 u = -(k_x^2 + k_y^2)\,u
```

Insert in the wave equation:

```math
-\omega^2 u = -c^2 (k_x^2 + k_y^2)\,u
\quad\Longleftrightarrow\quad
\omega^2 = c^2 (k_x^2 + k_y^2)
\quad\Longleftrightarrow\quad
\omega = c\sqrt{k_x^2 + k_y^2}
```

So (1.6) satisfies the wave equation exactly, for any $k_x, k_y$, provided $\omega$ is given by this dispersion relation.

## 1.2.4 Dispersion coefficient

$k_x = k_y = k$. Insert $u^n_{ij} = e^{\imath(kh(i+j) - \tilde{\omega}\, n\Delta t)}$ into (1.3). Shifting an index is a phase factor:

```math
u^{n\pm 1}_{ij} = e^{\mp\imath\tilde{\omega}\Delta t}\, u^n_{ij}, \qquad
u^n_{i\pm 1,j} = e^{\pm\imath kh}\, u^n_{ij}, \qquad
u^n_{i,j\pm 1} = e^{\pm\imath kh}\, u^n_{ij}
```

Using $e^{\imath\theta} - 2 + e^{-\imath\theta} = 2\cos\theta - 2 = -4\sin^2(\theta/2)$.

LHS of (1.3):

```math
\frac{u^{n+1}_{ij} - 2u^n_{ij} + u^{n-1}_{ij}}{\Delta t^2}
= \frac{e^{-\imath\tilde{\omega}\Delta t} - 2 + e^{\imath\tilde{\omega}\Delta t}}{\Delta t^2}\, u^n_{ij}
= -\frac{4}{\Delta t^2}\sin^2\!\Big(\frac{\tilde{\omega}\Delta t}{2}\Big)\, u^n_{ij}
```

RHS of (1.3): each of the two spatial differences in the same way (identical since $k_x = k_y$):

```math
\frac{u^n_{i+1,j} - 2u^n_{ij} + u^n_{i-1,j}}{h^2}
= \frac{e^{\imath kh} - 2 + e^{-\imath kh}}{h^2}\, u^n_{ij}
= -\frac{4}{h^2}\sin^2\!\Big(\frac{kh}{2}\Big)\, u^n_{ij}
\quad\Rightarrow\quad
\text{RHS} = -\frac{8c^2}{h^2}\sin^2\!\Big(\frac{kh}{2}\Big)\, u^n_{ij}
```

LHS = RHS, divide by $-4u^n_{ij}$:

```math
\frac{1}{\Delta t^2}\sin^2\!\Big(\frac{\tilde{\omega}\Delta t}{2}\Big)
= \frac{2c^2}{h^2}\sin^2\!\Big(\frac{kh}{2}\Big)
\quad\Longleftrightarrow\quad \sin^2\!\Big(\frac{\tilde{\omega}\Delta t}{2}\Big) = 2C^2 \sin^2\!\Big(\frac{kh}{2}\Big)
```
Now $C = 1/\sqrt{2}$ $\Rightarrow$ $2C^2 = 1$:

```math
\sin^2\!\Big(\frac{\tilde{\omega}\Delta t}{2}\Big) = \sin^2\!\Big(\frac{kh}{2}\Big)
\quad\Rightarrow\quad
\tilde{\omega}\Delta t = kh
```

and with $\Delta t = Ch/c = h/(\sqrt{2}\,c)$:

```math
\tilde{\omega} = \frac{kh}{\Delta t} = \sqrt{2}\,c\,k = c\sqrt{k^2 + k^2} = c\sqrt{k_x^2 + k_y^2} = \omega.
```

## 1.2.5 Genrate the Gif


