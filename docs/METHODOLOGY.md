# Mathematical and computational methodology

## Model and assumptions

The computational domain is the unit square $\Omega=(0,1)^2$. The main model is

$$
\frac{\partial u}{\partial t}-\alpha\Delta u=f(x,y),\qquad
u(x,y,0)=0,\qquad u|_{\partial\Omega}=0.
$$

The source $f$ is spatially varying but constant in time. The principal experiments use $\alpha=0.1$ and final time $T=0.1$. The PDE is linear, and experiments use synthetic observations, not instrument readings.

## From Taylor series to the FTCS update

Write $U_{i,j}^n\approx u(x_i,y_j,t_n)$ and let $h_x,h_y,\Delta t$ be the mesh widths. Taylor expansions give

$$
\frac{u(t+\Delta t)-u(t)}{\Delta t}=u_t+O(\Delta t),
$$

$$
\frac{u(x+h_x)-2u(x)+u(x-h_x)}{h_x^2}=u_{xx}+O(h_x^2),
$$

with an analogous centered difference in $y$. Substitution into the heat equation gives

$$
U_{i,j}^{n+1}=U_{i,j}^n+r_x(U_{i+1,j}^n-2U_{i,j}^n+U_{i-1,j}^n)
+r_y(U_{i,j+1}^n-2U_{i,j}^n+U_{i,j-1}^n)+\Delta t F_{i,j},
$$

where $r_x=\alpha\Delta t/h_x^2$ and $r_y=\alpha\Delta t/h_y^2$. Its stability restriction is $r_x+r_y\le 1/2$; for $h_x=h_y=h$, $\alpha\Delta t/h^2\le 1/4$.

The formal consistency order for smooth solutions is $O(\Delta t+h_x^2+h_y^2)$. The independent sinusoidal verification in Notebook 09 tests spatial and temporal orders separately. This is not a universal convergence guarantee for nonsmooth sources.

## Parameterized source and observation operator

Approximate the source by $f_\theta=\sum_{k=1}^{p}\theta_k\phi_k$, with Gaussian source functions $\phi_k$. At the interior spatial nodes, write $F=B\theta$. Then the state recurrence is

$$
U^{n+1}=MU^n+\Delta t B\theta,\qquad M=I+\alpha\Delta t L_h.
$$

Extract temperatures at sensor positions and selected times; stacking these produces the predicted observation vector

$$
\widehat d=A\theta.
$$

**Notation:** $d$ is the measured/noisy observation vector; $\widehat d$ is the vector predicted from trial or reconstructed coefficients. A saved variable named `d_predicted` is the same kind of quantity as $\widehat d$; it is not a second mathematical model.

## Reconstruction and regularization

For $m$ scalar observations, the quadratic objective used in the core inverse experiments is

$$
J_\lambda(\theta)=\frac{1}{2m}\lVert A\theta-d\rVert_2^2+\frac{\lambda}{2}\lVert\theta\rVert_2^2.
$$

The minimizer satisfies $(A^TA+m\lambda I)\widehat\theta=A^Td$. For $\lambda>0$, the finite-dimensional objective is strictly convex; that does not imply the original physical source has been identified uniquely from the measurements.

**Different errors:** coefficient error $\lVert\widehat\theta-\theta^*\rVert_2/\lVert\theta^*\rVert_2$ is not source-field error $\lVert B\widehat\theta-F^*\rVert_2/\lVert F^*\rVert_2$. These are reported separately.

## Discrete adjoint

For the discrete recurrence and the objective evaluated at the chosen observation times, define $g_n$ as the transposed observation residual injection (zero at unobserved times). With terminal convention $p^{N_t+1}=0$,

$$
p^n=M^Tp^{n+1}+g_n\quad(n=N_t,\ldots,1),
$$

and the source-parameter gradient is

$$
\nabla J_\lambda(\theta)=\Delta t B^T\sum_{n=1}^{N_t}p^n+\lambda\theta.
$$

The original FTCS Laplacian with homogeneous Dirichlet boundaries is symmetric. The gradient was checked using both central numerical differences and the independent explicit response matrix. See [RESULTS_SUMMARY.md](RESULTS_SUMMARY.md) and Notebook 08.

## Scope

This repository demonstrates numerical verification of established techniques. It does **not** claim a novel reconstruction algorithm, proof of general sensor optimality, or a universal adjoint speedup. See the full technical report for the longer derivations and literature.
