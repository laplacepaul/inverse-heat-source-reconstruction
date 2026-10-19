---
title: "Numerical Analysis and PDE-Constrained Optimization for Inverse Heat-Source Reconstruction"
subtitle: "An independent master's-level computational mathematics investigation"
date: "19 October 2026"
lang: en-US
geometry: margin=1.0in
fontsize: 12pt
colorlinks: true
linkcolor: blue
---

# Abstract

Inverse heat-source reconstruction estimates an unknown spatial heating distribution using temperature measurements governed by a diffusion equation. The forward model may be accurate and stable even when the inverse reconstruction is poorly conditioned or nonunique. This study develops a two-dimensional heat equation solver using the Forward-Time Central-Space (FTCS) scheme and reconstructs time-independent sources represented by Gaussian basis functions. The discrete source-to-observation operator gives least-squares and Tikhonov-regularized estimation problems; a discrete adjoint is derived as an alternative gradient-evaluation strategy. Nine computational notebooks investigate analytic forward verification, superposition, observation matrices, reconstruction accuracy, noise and sensor sensitivity, larger underdetermined systems, source-model mismatch, adjoint gradients, and independent convergence and computational benchmarks. Spatial and temporal reference tests exhibit approximately second- and first-order accuracy, respectively. The well-conditioned three-coefficient model attains 0.5684% source-field error, while a 25-parameter model exhibits severe noise sensitivity without regularization. Source mismatch persists over the tested mesh-refinement range. Adjoint gradients agree with the explicit matrix gradient to approximately $7.98\times10^{-16}$ relative error; nevertheless, a direct matrix workflow is faster in the tested small linear systems. The work provides a documented and reproducible investigation of established numerical techniques, not a claim of novel algorithmic theory.

**Keywords:** inverse problems; heat equation; finite differences; Tikhonov regularization; PDE-constrained optimization; discrete adjoints; numerical conditioning; computational verification.

# Reading Guide and Mathematical Notation

This report distinguishes three levels of description: the continuous heat equation, its discrete numerical approximation, and the observation-based inverse problem. A symbol is not an additional algorithm merely because it is introduced in the derivation. For example, $d$ is the vector of measured temperatures and $\widehat d=A\widehat\theta$ is the model's predicted temperature vector after reconstruction; Notebook 4 uses the computational expression `A @ theta_ridge` for this prediction.

| Quantity | Mathematical meaning | Typical computational form |
|---|---|---|
| $u(x,y,t)$ | Continuous temperature | Analytic solution or PDE reference |
| $U^n$ | Interior temperature array at step $n$ | Forward-solver state |
| $f(x,y)$, $F$ | Continuous and grid-sampled source | Source function or array |
| $B$ | Gaussian source dictionary on the numerical grid | Basis-field matrix |
| $\theta$, $\widehat\theta$ | Trial and reconstructed source coefficients | `theta`, `theta_ridge` |
| $d$ | Measured (possibly noisy) temperature samples | Observation data vector |
| $\widehat d$ | Temperature samples predicted from fitted coefficients | `A @ theta_ridge` |
| $A$ | Linear observation-response matrix | `A` |
| $m,p,q$ | Scalar measurements, coefficients, interior state nodes | Integer dimensions |
| $H_n$ | Sampling matrix at observation step $n$ | Sensor selection / interpolation |
| $p^n$ | Adjoint *vector* at time index $n$ | Backward-propagated state |

The integer $p$ (basis dimension) and the vector $p^n$ (adjoint state) use customary mathematical notation but represent different objects; dimensions and context distinguish them. In the publication version, all displayed mathematics is rendered by LaTeX rather than pasted as special Unicode math characters.

# 1. Introduction

## 1.1 Motivation and problem statement

Heat conduction is often modeled as a forward problem: a heat source and boundary/initial conditions are prescribed, and a temperature field is predicted. In an inverse source problem, temperature measurements are supplied and the spatial source must be inferred. Such identification is relevant to thermal diagnostics and monitoring, but diffusion smooths spatial structure and reduces sensitivity to certain source features. Sparse or noisy measurements can therefore be compatible with many distinct source distributions. This distinction between forward numerical solvability and inverse identifiability motivates the present study.

The project considers a bounded square domain with fixed zero boundary temperature, zero initial temperature, constant thermal diffusivity, and a source independent of time. Synthetic measurements permit direct comparison with the known generating source. The scope is intentionally focused: we examine established finite-difference, regularization, and optimization methods with transparent numerical validation. This is neither a new PDE solver nor a physical experiment involving real sensors.

## 1.2 Aim and research objectives

The aim is to determine how the accuracy of a discretized heat equation, the choice of source basis, observations, regularization, and optimization strategy interact in finite-dimensional heat-source reconstruction. The specific objectives are to: derive and validate the FTCS forward scheme; construct Gaussian source models and the discrete observation-response matrix; compare unregularized and regularized inversion; investigate singular values, sensor placement, and noisy data; examine higher-dimensional and underdetermined cases; isolate discretization from model-mismatch effects; derive and verify a discrete adjoint; and compare complete matrix and adjoint computational workflows.

## 1.3 Scope and structure

The mathematics and related literature are reviewed in Chapter 2. Chapter 3 states the continuous problem; Chapter 4 derives the numerical scheme and examines stability and convergence. Chapters 5 and 6 develop matrix-based inversion and adjoint gradients. Chapter 7 reports experiment-specific results, Chapter 8 interprets them, and Chapter 9 concludes. Appendix A documents reproducibility and dataset provenance. Throughout the report $N$ denotes the number of *intervals* per spatial coordinate, not the number of grid points.

# 2. Literature Review and Mathematical Background

## 2.1 Parabolic diffusion and source identification

For a sufficiently regular spatial source, the heat equation defines a linear source-to-state map. Writing $\Delta_D$ for the Dirichlet Laplacian and taking zero initial temperature, the solution can be represented formally as

$$
u(t)=\int_0^t \exp\!\bigl(\alpha(t-s)\Delta_D\bigr)f\,ds. \tag{2.1}
$$

This integral conveys why source identification is distinct from direct temperature estimation: measured temperatures encode a diffused, temporally accumulated version of the source. LeVeque [1] provides the finite-difference framework for the heat equation and related evolution problems.

## 2.2 Well-posedness, singular values, and regularization

A classical well-posed problem has existence, uniqueness, and continuous dependence on data. Inverse problems may fail uniqueness or stability, although *not every* finite-dimensional or infinite-dimensional source-identification configuration does so. In the chosen basis, the key discrete question is whether the observation matrix has independent, adequately observable columns. With $A=U\Sigma V^{\mathsf T}$, a minimum-norm least-squares solution has the expansion

$$\widehat\theta=A^\dagger d=\sum_{j=1}^r\frac{u_j^{\mathsf T}d}{\sigma_j}v_j. \tag{2.2}$$

Small positive singular values amplify the relevant noise components; a null space introduces exact nonuniqueness. Foundational discussions of ill-posedness and numerical regularization are given by Engl, Hanke and Neubauer [2], Hansen [3], and Tikhonov and Arsenin [4]. Vogel [5] treats computational aspects, including regularization-parameter selection.

## 2.3 Discretization, source bases, and sensor design

Finite differences approximate derivatives using nearby sampled values. Taylor expansion establishes local approximation orders; the stability of the time integrator determines whether discretization error can be controlled over finite time intervals [1]. Source expansion in a finite dictionary reduces a functional inversion problem to parameter estimation. Its representational adequacy and the geometry of the observation operator must be considered together: adding coefficients increases flexibility but may introduce poorly observable directions.

## 2.4 PDE-constrained optimization and adjoints

PDE-constrained optimization treats the state equation as a constraint on an objective involving observations and source parameters. Foundational treatments include Lions [6] and Hinze et al. [7]. Giles and Pierce [8] develop adjoint sensitivities from both linear algebra and PDE viewpoints. Particularly relevant is Adavani and Biros [9], who study spatial source identification constrained by linear parabolic PDEs, $L^2$ Tikhonov regularization, and multigrid algorithms. Their work is more advanced and focuses on scalable solvers; our investigation instead emphasizes transparent construction, gradient validation, and small-scale empirical comparison. Standard gradient and quasi-Newton optimization principles follow Nocedal and Wright [10].

## 2.5 Position of the present study

This project does not claim a previously unknown connection between PDEs and regularization. Its contribution is a *controlled computational synthesis* of established theory, with preserved experiments that separate spatial and temporal verification, report errors with explicit definitions, and compare direct and adjoint methods on the *complete* reconstruction task.

# 3. Mathematical Formulation

## 3.1 Continuous forward problem

Let $\Omega=(0,1)\times(0,1)$, $0<t\le T$, and let $u(x,y,t)$ be the temperature relative to a prescribed reference. We solve

$$u_t-\alpha\Delta u=f(x,y),\qquad \alpha>0,\quad \Delta=\partial_{xx}+\partial_{yy}. \tag{3.1}$$

The initial and boundary conditions are

$$u(x,y,0)=0,\qquad u|_{\partial\Omega}=0. \tag{3.2}$$

The common benchmark parameters are $\alpha=0.1$ and $T=0.1$. The source is time-independent, and its units are those implied by the normalized model; no calibrated physical thermal units are asserted.

## 3.2 Analytical sinusoidal benchmark

For $f(x,y)=\sin(\pi x)\sin(\pi y)$, set $u=a(t)\sin(\pi x)\sin(\pi y)$. Since $\Delta f=-2\pi^2 f$, the PDE reduces to

$$a'(t)+2\alpha\pi^2 a(t)=1,\qquad a(0)=0. \tag{3.3}$$

Hence

$$u(x,y,t)=\frac{1-e^{-2\alpha\pi^2t}}{2\alpha\pi^2}\sin(\pi x)\sin(\pi y). \tag{3.4}$$

This solution is a verification reference, not the source model used in every subsequent experiment.

## 3.3 Gaussian source and measurements

Approximate the source by $p$ Gaussian basis functions:

$$f_\theta(x,y)=\sum_{\ell=1}^{p}\theta_\ell\phi_\ell(x,y),\quad
\phi_\ell(x,y)=\exp\!\left[-\frac{(x-a_\ell)^2+(y-b_\ell)^2}{2s_\ell^2}\right]. \tag{3.5}$$

The observation operator samples temperature at selected locations and times. Unless otherwise specified, sensors select interior *grid nodes*; the cross-grid experiments in Chapter 7 instead interpolate observations between different grids. Let $m$ be the number of scalar observations, $p$ the source coefficient count, and $q$ the number of interior state degrees of freedom. The finite inverse task is to recover $\theta\in\mathbb R^p$ from $d\in\mathbb R^m$.

# 4. Numerical Discretization, Stability, and Convergence

## 4.1 Space–time mesh

Choose $N_x,N_y$ equal spatial intervals, with $h_x=1/N_x$, $h_y=1/N_y$ and nodes $x_i=ih_x$, $y_j=jh_y$. With $N_t$ time steps, $\Delta t=T/N_t$ and $t_n=n\Delta t$. The approximation $U_{ij}^{n}\approx u(x_i,y_j,t_n)$ is computed at interior nodes, with $q=(N_x-1)(N_y-1)$. The source samples are $F_{ij}=f(x_i,y_j)$.

## 4.2 Taylor-series derivation of the finite differences

For a sufficiently smooth function of time, Taylor expansion around $t_n$ gives

$$u(t_n+\Delta t)=u(t_n)+\Delta t\,u_t(t_n)+\tfrac12\Delta t^2u_{tt}(t_n)+O(\Delta t^3),$$

therefore

$$u_t(t_n)=\frac{u(t_n+\Delta t)-u(t_n)}{\Delta t}+O(\Delta t). \tag{4.1}$$

Similarly, expansion around $x_i$ in the positive and negative directions yields

$$\begin{aligned}
u(x_i+h_x)&=u_i+h_xu_x+\tfrac12h_x^2u_{xx}+\tfrac16h_x^3u_{xxx}+\tfrac1{24}h_x^4u_{xxxx}+O(h_x^5),\\
u(x_i-h_x)&=u_i-h_xu_x+\tfrac12h_x^2u_{xx}-\tfrac16h_x^3u_{xxx}+\tfrac1{24}h_x^4u_{xxxx}+O(h_x^5).
\end{aligned}$$

Adding cancels the odd-order terms and gives

$$u_{xx}(x_i)=\frac{u_{i+1}-2u_i+u_{i-1}}{h_x^2}-\frac{h_x^2}{12}u_{xxxx}(x_i)+O(h_x^4). \tag{4.2}$$

The same construction in $y$ gives the centered second difference with error $O(h_y^2)$. Substituting the temporal and spatial differences directly into the PDE yields

$$\begin{aligned}
\frac{U_{ij}^{n+1}-U_{ij}^n}{\Delta t}={}&\alpha\Bigl[\frac{U_{i+1,j}^n-2U_{ij}^n+U_{i-1,j}^n}{h_x^2}\\
&+\frac{U_{i,j+1}^n-2U_{ij}^n+U_{i,j-1}^n}{h_y^2}\Bigr]+F_{ij}.
\end{aligned} \tag{4.3}$$

Defining $r_x=\alpha\Delta t/h_x^2$ and $r_y=\alpha\Delta t/h_y^2$, isolate the new state:

$$\boxed{U_{ij}^{n+1}=(1-2r_x-2r_y)U_{ij}^n+r_x(U_{i-1,j}^n+U_{i+1,j}^n)+r_y(U_{i,j-1}^n+U_{i,j+1}^n)+\Delta t F_{ij}.} \tag{4.4}$$

For a square grid $h_x=h_y=h$ and $r=\alpha\Delta t/h^2$, the center coefficient is $1-4r$ and each neighbour has coefficient $r$. The initial values are zero and all boundary nodes are maintained at zero; the update applies only to the interior. With the Dirichlet discrete Laplacian $L_h$, the vector form is

$$U^{n+1}=MU^n+\Delta t F,\qquad M=I+\alpha\Delta t L_h. \tag{4.5}$$

## 4.2.1 Deriving the forward time difference with an explicit remainder

For fixed $(x_i,y_j)$, Taylor's theorem for a sufficiently smooth temperature field gives

$$
 u(x_i,y_j,t_n+\Delta t)=u(x_i,y_j,t_n)+\Delta t\,u_t(x_i,y_j,t_n)+\frac{\Delta t^2}{2}u_{tt}(x_i,y_j,t_n)+O(\Delta t^3).
$$

Dividing the difference of successive temperature values by $\Delta t$ yields

$$
\frac{u(x_i,y_j,t_n+\Delta t)-u(x_i,y_j,t_n)}{\Delta t}=u_t+\frac{\Delta t}{2}u_{tt}+O(\Delta t^2).
$$

Thus the forward difference is not *equal* to $u_t$: it differs from $u_t$ by a first-order truncation term. We introduce the numerical approximation only after discarding that controlled remainder. The approximation is evaluated at $t_n$, while the unknown value is at the future time $t_{n+1}$. This asymmetry is why the temporal method has first-order, rather than second-order, accuracy.

## 4.2.2 Deriving both centered second differences

At fixed $y_j,t_n$, the expansions in the positive and negative $x$ directions are

$$
\begin{aligned}
 u(x_i+h_x,y_j,t_n)&=u_{ij}^n+h_xu_x+\frac{h_x^2}{2}u_{xx}+\frac{h_x^3}{6}u_{xxx}+\frac{h_x^4}{24}u_{xxxx}+O(h_x^5),\\
 u(x_i-h_x,y_j,t_n)&=u_{ij}^n-h_xu_x+\frac{h_x^2}{2}u_{xx}-\frac{h_x^3}{6}u_{xxx}+\frac{h_x^4}{24}u_{xxxx}+O(h_x^5).
\end{aligned}
$$

Adding eliminates the first and third derivatives. Assuming sufficient smoothness to continue the even-power expansion,

$$
 u(x_i+h_x)-2u(x_i)+u(x_i-h_x)=h_x^2u_{xx}+\frac{h_x^4}{12}u_{xxxx}+O(h_x^6).
$$

After division by $h_x^2$,

$$
 u_{xx}=\frac{u_{i+1,j}^n-2u_{ij}^n+u_{i-1,j}^n}{h_x^2}-\frac{h_x^2}{12}u_{xxxx}+O(h_x^4).
$$

The same reasoning in $y$ gives

$$
 u_{yy}=\frac{u_{i,j+1}^n-2u_{ij}^n+u_{i,j-1}^n}{h_y^2}-\frac{h_y^2}{12}u_{yyyy}+O(h_y^4).
$$

These formulas explain the phrase *central space*: both neighbours of an interior grid point enter symmetrically. The discrete Laplacian therefore uses five points in two dimensions: the centre and its four nearest axial neighbours.

## 4.2.3 Discrete equation, rearrangement, and algorithmic interpretation

Evaluate the PDE at $(x_i,y_j,t_n)$, replace the continuous derivatives by the approximations above, and write the discrete source value as $F_{ij}=f(x_i,y_j)$. Before rearrangement, the numerical equation is

$$
\begin{aligned}
\frac{U_{ij}^{n+1}-U_{ij}^n}{\Delta t}
=\alpha\Bigl[
\frac{U_{i+1,j}^n-2U_{ij}^n+U_{i-1,j}^n}{h_x^2}
+\frac{U_{i,j+1}^n-2U_{ij}^n+U_{i,j-1}^n}{h_y^2}
\Bigr]+F_{ij}.
\end{aligned}
$$

Multiplying through by $\Delta t$ and moving $U_{ij}^n$ to the right produces the explicit update (4.4). The scheme is *explicit* because all temperatures on the right-hand side are known at time $n$; solving a simultaneous system of linear equations is not required for each time step. The source contribution is $\Delta tF_{ij}$ and accumulates over successive steps, while the diffusion contribution averages information from neighbouring grid values.

At the edges of the square, the Dirichlet data enforce $U=0$. The update should be computed only for indices $1\le i\le N_x-1$ and $1\le j\le N_y-1$; applying the interior formula indiscriminately to boundary nodes would solve a different problem.

## 4.3 Local truncation error and consistency

Substituting the *exact smooth solution* into the discrete residual gives

$$\tau_{ij}^n=\tfrac{\Delta t}{2}u_{tt}-\tfrac{\alpha h_x^2}{12}u_{xxxx}-\tfrac{\alpha h_y^2}{12}u_{yyyy}+O(\Delta t^2+h_x^4+h_y^4). \tag{4.6}$$

Thus the residual is $O(\Delta t+h_x^2+h_y^2)$: first-order in time and second-order in space. This establishes consistency, not convergence by itself.

## 4.4 Von Neumann stability condition

For the homogeneous error recurrence insert $e_{ij}^n=G^n\exp[\mathrm i(i\xi+j\eta)]$. The amplification factor is

$$G=1-4r_x\sin^2(\xi/2)-4r_y\sin^2(\eta/2). \tag{4.7}$$

Requiring $|G|\le1$ for all Fourier modes gives the standard sufficient restriction

$$\boxed{r_x+r_y\le\tfrac12.} \tag{4.8}$$

For $h_x=h_y=h$, this becomes $r\le1/4$, or $\Delta t\le h^2/(4\alpha)$. The finite homogeneous Dirichlet implementation satisfies the corresponding stability requirement. Refining $h$ by a factor of two generally requires a timestep about four times smaller if one remains near the explicit stability limit.

### 4.4.1 Fourier-mode calculation and the mesh-dependent timestep

Starting from the homogeneous error equation, substitute $e_{ij}^n=G^n e^{\mathrm{i}(i\xi+j\eta)}$. Dividing by the nonzero common factor gives

$$G=1-2r_x-2r_y+r_x(e^{\mathrm{i}\xi}+e^{-\mathrm{i}\xi})+r_y(e^{\mathrm{i}\eta}+e^{-\mathrm{i}\eta}).$$

Euler's identity and $1-\cos z=2\sin^2(z/2)$ give the amplification factor shown in (4.7). Since each squared sine lies in $[0,1]$, the largest amplification factor is $1$ and the smallest is $1-4r_x-4r_y$. Requiring that all Fourier components remain bounded in amplitude is equivalent to

$$-1\le 1-4(r_x+r_y)\le 1,\qquad r_x+r_y\le\tfrac12.$$

For equal grid spacing, the restriction reads $\Delta t\le h^2/(4\alpha)$. This explains why the temporal resolution cannot be selected independently of the spatial resolution when using an explicit two-dimensional heat solver. A grid refinement that halves $h$ reduces the maximum stable step by a factor of four. These statements concern the error amplification of the homogeneous scheme; source forcing and boundary effects must be interpreted within the same fixed finite observation interval.

### 4.4.2 Consistency, stability, and convergence are different statements

Consistency follows from substituting the exact smooth solution into the discrete stencil: the local residual vanishes as $h_x,h_y,\Delta t\to0$. Stability controls whether errors are amplified by time stepping. Convergence concerns whether the *computed numerical solution* tends to the PDE solution. For the linear initial-boundary-value problem and regularity assumptions used in this study, consistency together with suitable stability yields the expected global error bound of order $O(\Delta t+h_x^2+h_y^2)$. The separate experimental rates in Notebook 9 corroborate that expectation for a smooth sinusoidal mode; they do not replace the assumptions of the mathematical estimate.

## 4.5 Analytical verification and separate convergence tests

The first notebook compares the numerical solution with (3.4); joint mesh–time refinement yields relative errors 0.01182374, 0.00322413, and 0.00080473 on grids with 12, 24, and 48 subdivisions, respectively. Joint observed orders are 1.875 and 2.003. Because $\Delta t$ changes with $h$, these do *not* isolate temporal order.

Notebook 9 separates the two contributions. The sampled sinusoidal mode $S_{ij}=\sin(\pi x_i)\sin(\pi y_j)$ satisfies $L_hS=-\mu_hS$ with

$$\mu_h=\frac8{h^2}\sin^2\left(\frac{\pi h}{2}\right),\qquad a_h(t)=\frac{1-e^{-\alpha\mu_h t}}{\alpha\mu_h}. \tag{4.9}$$

Its exact semi-discrete amplitude $a_h$ is compared with the continuous amplitude $a$ from (3.4) to measure **pure spatial error**, not the 2000-step FTCS approximation. For the temporal study, $N=16$ is held fixed and FTCS is compared against $a_h(T)$; four stable configurations with 32–256 steps are used. Table 7.7 supplies the measured values. For this smooth benchmark, the results support spatial order approximately 2 and temporal order approximately 1; they do not constitute a theorem for arbitrary discontinuous forcing.

# 5. Inverse Reconstruction and Regularization

## 5.1 From numerical state to observations

Let $B\in\mathbb R^{q\times p}$ contain the sampled Gaussian basis fields as columns, so $F=B\theta$. The recurrence becomes

$$U^{n+1}=MU^n+\Delta t B\theta,\qquad U^0=0. \tag{5.1}$$

By iteration,

$$U^n=\Delta t\sum_{k=0}^{n-1}M^kB\theta. \tag{5.2}$$

If $C_n$ selects the measured grid locations at time $n$, stack the time-indexed measurements to obtain

$$\widehat d=A\theta,\qquad A=\begin{bmatrix}
C_{n_1}\Delta t\sum_{k=0}^{n_1-1}M^kB\\
\vdots\\
C_{n_o}\Delta t\sum_{k=0}^{n_o-1}M^kB
\end{bmatrix}. \tag{5.3}$$

Here the number of observation *times* is $n_o$ in this display; the total number of scalar observations remains $m$. Operationally, each column of $A$ is formed by running the forward solver with one basis function as the source, then stacking temperatures at specified observation locations and times. Thus explicit matrix powers are unnecessary. For the initial 3-basis experiment there are 16 spatial sensors and three times, hence $A\in\mathbb R^{48\times3}$. The measured numerical rank is 3 and $\kappa_2(A)\approx1.806359$. Direct forward measurements and $A\theta^\star$ agree to about $2.08\times10^{-17}$ in the reported comparison.

## 5.1.1 Why the response matrix can be constructed column by column

The basis representation writes the unknown spatial source as $F=B\theta=\sum_{\ell=1}^{p}\theta_\ell B_{:,\ell}$. Because the discrete heat recurrence is linear in $F$ and has zero initial state, its prediction for a sum of sources is the same sum of the individual predictions. Therefore, to construct column $\ell$ of $A$, set $\theta=e_\ell$, solve the PDE once, and collect the sensor measurements at the required times. Stacking those observations gives $A_{:,\ell}$. Repeating this procedure for each basis function gives $A$ without explicitly computing $M^k$ or forming the full state-to-state matrix powers. Notebook 3 independently checks this numerical superposition identity against a full forward simulation of a combined source.

The observation operator requires care in experiments where measurements are generated on a finer grid and compared with a different inversion grid. There, measurement interpolation is part of the observation map; the sensor-selection matrix $C_n$ used in the simplest derivation must be replaced by the appropriate grid/interpolation operator.

## 5.1.2 Observed data versus predicted observations

The observation vector $d$ holds the temperatures supplied to the inverse solver. When synthetic noise is present, $d=A\theta^\star+\varepsilon$ in a matched forward model. For a candidate parameter vector $\theta$, the corresponding prediction is $d(\theta)=A\theta$. After fitting, define

$$\boxed{\widehat d=A\widehat\theta.}$$

This $\widehat d$ is the *predicted observation vector*, not an extra measurement or source term. It is the mathematical equivalent of Notebook 4's expression `A @ theta_ridge`. The discrepancy between prediction and observation is $\widehat d-d$, whose size enters the inverse objective and residual diagnostic. In a model-mismatch experiment, the identity $d=A\theta^\star+\varepsilon$ is not sufficient: additional discrepancy terms are required to represent an unmodelled source feature or a different forward discretization.

## 5.2 Noise model and error definitions

Synthetic data have the form $d=A\theta^\star+\varepsilon$ with Gaussian noise $\varepsilon\sim\mathcal N(0,\sigma^2 I)$. For several experiments, $\sigma=\eta\|d_{\mathrm{clean}}\|_\infty$ for prescribed $\eta$; this *is not the same* as fixing $\|\varepsilon\|_2/\|d_{\mathrm{clean}}\|_2=\eta$. The experiment-specific seed and scaling are recorded in the datasets and notebook code.

The two principal error metrics are distinct:

$$E_\theta=\frac{\|\widehat\theta-\theta^\star\|_2}{\|\theta^\star\|_2},\qquad E_f=\frac{\|B\widehat\theta-F^\star\|_2}{\|F^\star\|_2}. \tag{5.4}$$

Field errors use the same sampled spatial grid in numerator and denominator; a common cell-area factor cancels on uniform grids. Observation fit is separately measured by $\|A\widehat\theta-d\|_2/\|d\|_2$. Small data residual is not proof of small source error.

## 5.3 Least squares and regularization

Unregularized least squares minimizes $\|A\theta-d\|_2^2/(2m)$. Setting its gradient to zero gives $A^{\mathsf T}A\theta=A^{\mathsf T}d$. For ill-conditioned problems, one should prefer QR/SVD-based `lstsq` to explicitly inverting $A^{\mathsf T}A$, since forming normal equations can square the spectral condition number.

Zero-order Tikhonov regularization solves

$$\boxed{\widehat\theta_\lambda=\arg\min_\theta\left[\frac{\|A\theta-d\|_2^2}{2m}+\frac\lambda2\|\theta\|_2^2\right].} \tag{5.5}$$

The optimality system is

$$\left(A^{\mathsf T}A+m\lambda I\right)\widehat\theta_\lambda=A^{\mathsf T}d. \tag{5.6}$$

For $\lambda>0$ this coefficient-space objective is strictly convex. Its SVD filter factors are $\sigma_j^2/(\sigma_j^2+m\lambda)$, damping weakly observed directions. Critically, $\|\theta\|_2^2$ is a **coefficient norm**, not a spatial smoothness norm $\|B\theta\|_2^2$. The effect of the penalty therefore depends on the basis representation and scaling [2–5].

## 5.4 Rank, sparse observations, and mismatch

For $m<p$, rank$(A)\le m<p$ and a nontrivial null space exists; the data alone cannot identify all coefficients. The 18-by-25 configuration has rank 18 and nullity 7. A ratio using its positive singular values may still be finite, but must not be misrepresented as unique invertibility. Source mismatch can be represented as $F^\star=B\theta^\star+\rho$, where $\rho$ is not represented by the chosen basis. In cross-grid experiments, observations are generated on a distinct mesh and interpolated, so discretization and representational discrepancies need careful separation.

# 6. Discrete Adjoint Optimization

## 6.1 Constrained objective

For observation times $\mathcal O$, state matrices $H_n$ and residual data $d_n$, consider

$$J(\theta,U)=\frac1{2m}\sum_{n\in\mathcal O}\|H_nU^n-d_n\|_2^2+\frac\lambda2\|\theta\|_2^2, \tag{6.1}$$

subject to $U^n-MU^{n-1}-\Delta t B\theta=0$, $n=1,\dots,N_t$, with $U^0=0$. To derive the **discrete** adjoint, introduce multipliers $p^n$ and the sign convention

$$\mathcal L=J-\sum_{n=1}^{N_t}(p^n)^{\mathsf T}(U^n-MU^{n-1}-\Delta t B\theta). \tag{6.2}$$

Let $g_n=H_n^{\mathsf T}(H_nU^n-d_n)/m$ when $n\in\mathcal O$ and zero otherwise. Stationarity with respect to $U^n$ gives the backward recursion

$$\boxed{p^{N_t+1}=0,\qquad p^n=M^{\mathsf T}p^{n+1}+g_n\quad(n=N_t,\ldots,1).} \tag{6.3}$$

Differentiating with respect to $\theta$ yields

$$\boxed{\nabla J(\theta)=\Delta t\,B^{\mathsf T}\sum_{n=1}^{N_t}p^n+\lambda\theta.} \tag{6.4}$$

The implemented five-point Dirichlet heat propagator is symmetric, so $M^{\mathsf T}=M$; this simplification is *not* generally justified for other discretizations. The implementation injects residuals at the correct observation steps before accumulating the adjoint fields.

## 6.1.1 Derivation of the backward adjoint recurrence

The discrete state constraint at step $n$ is $U^n-MU^{n-1}-\Delta t B\theta=0$. With the Lagrangian sign convention in (6.2), a temperature state $U^n$ appears in the observation misfit, negatively in its own constraint, and positively through $MU^n$ in the next constraint. Define

$$g_n=\frac1m H_n^{\mathsf T}(H_n U^n-d_n)$$

at observed time levels and $g_n=0$ otherwise. For $1\le n<N_t$, stationarity with respect to $U^n$ gives

$$0=\frac{\partial\mathcal L}{\partial U^n}=g_n-p^n+M^{\mathsf T}p^{n+1}.$$

At the last time, no later constraint exists, so $p^{N_t}=g_{N_t}$. Defining $p^{N_t+1}=0$ permits a single recurrence for every time index, exactly as in (6.3). The adjoint propagates backward because the current adjoint depends on the next adjoint time level.

For the derivative with respect to the source coefficients, the direct contribution from the penalty is $\lambda\theta$, and each constrained time step contributes $+\Delta t B^{\mathsf T}p^n$. Summing these terms gives (6.4). Thus the gradient of *all* $p$ source coefficients is assembled after one forward pass and one backward adjoint pass at a trial iterate; it is not obtained by finite-differencing each coefficient in the production optimizer.

## 6.1.2 Discretize-then-optimize versus optimize-then-discretize

The present implementation differentiates the discrete recurrence *actually used by the solver*. This is known as a discretize-then-optimize construction. Deriving a continuous adjoint PDE first and subsequently discretizing it may give a different discrete gradient unless the discretizations and boundary treatments are carefully matched. The near-machine-precision agreement between the response-matrix derivative and the coded adjoint is therefore meaningful evidence about consistency with our particular FTCS implementation. It is not a statement that every continuous and discrete adjoint construction would automatically coincide.

## 6.2 Three independent gradient calculations

The analytic discrete matrix gradient is

$$\nabla J_{A}(\theta)=A^{\mathsf T}(A\theta-d)/m+\lambda\theta. \tag{6.5}$$

A central-difference check uses

$$\left[\nabla J_{\mathrm{FD}}\right]_\ell\approx\frac{J(\theta+\epsilon e_\ell)-J(\theta-\epsilon e_\ell)}{2\epsilon}. \tag{6.6}$$

Notebook 8 records relative adjoint/finite-difference agreement $5.51575697\times10^{-11}$, adjoint/matrix agreement $7.98040339\times10^{-16}$, and maximum absolute finite-difference component discrepancy $6.0116\times10^{-14}$. The fitted 25-coefficient solution differs from the direct coefficient minimizer by about $5.09057\times10^{-5}$ relative norm. These checks support code correctness for the *same discrete objective*, not a guarantee of physical source accuracy.

## 6.3 Optimization and efficiency

The 25- and 100-coefficient adjoint optimizations converged successfully, using 179 and 211 iterations. Task 8 reports 4.9278s and 8.4803s in those runs; its one-gradient cost (100 coefficients) is 0.031145s and explicit matrix build cost 0.875801s. These *individual-operation* times cannot be interpreted as end-to-end speedups. Notebook 9 supplies separate median complete-workflow benchmarks under its own configurations. A direct coefficient solve is particularly competitive for the small linear models considered here. L-BFGS-B and relevant optimization theory are described in [10].

# 7. Numerical Experiments and Results

## 7.1 Data provenance and metrics

The canonical results are the uploaded corrected notebooks and files in `data/`, rather than the older initial ZIP. Task 5 reports **coefficient errors**. Tasks 6, 7, 8, and 9 principally report **source-field errors**. Timing results are specific to their benchmark protocol. Four actual PNG figures supplied with the project are embedded below. Displayed percentages are rounded only after extracting the underlying saved arrays.

## 7.2 Small response operator and reconstruction

The three-basis reference coefficients are $\theta^\star=(1.2,0.85,1.05)$. Notebook 3's $48\times3$ observation matrix has rank 3, condition number $1.806359$, spatial grid with 32 **subdivisions** ($33\times33$ nodes), and 46 time steps with $\Delta t\approx0.002173913$. Three observation step indices are 12, 23, and 46, hence the measurement times are approximately 0.026087, 0.050000, and 0.100000. The stored noisy observation standard deviation is about $0.00045105$ and realized relative noise norm about 2.19%.

| Coefficient | True | Ordinary LS | Tikhonov, $\lambda=10^{-7}$ |
|---|---:|---:|---:|
| $\theta_1$ | 1.200000 | 1.203937 | 1.202577 |
| $\theta_2$ | 0.850000 | 0.853090 | 0.852704 |
| $\theta_3$ | 1.050000 | 1.062738 | 1.060681 |

The regularized **coefficient error** is 0.6262%, while **source-field error** is 0.5684%. Agreement with the direct solution is corroborated by BFGS minimization. The small errors pertain to this well-conditioned and matched-basis example.

## 7.3 Noise and sensor robustness (Notebook 5)

Task 5 stored 150 independent trials for each of five noise-standard-deviation levels, with Gaussian SD defined relative to maximum clean observation magnitude. The canonical mean **coefficient errors** calculated from the uploaded NPZ are:

| Prescribed noise SD / max signal | LS mean coefficient error | Ridge mean coefficient error |
|---:|---:|---:|
| 0% | 0.0000% | 0.1368% |
| 1% | 0.6567% | 0.6691% |
| 3% | 1.9225% | 1.9214% |
| 5% | 3.4625% | 3.4695% |
| 10% | 6.8161% | 6.8060% |

Ridge has a small penalty bias at nominal zero noise and provides similar performance to ordinary least squares for this well-conditioned 3-basis case. More challenging sensitivity appears when the sensors are reduced:

| Sensors | Positive-singular-value condition ratio | LS mean coefficient error | Ridge mean coefficient error |
|---:|---:|---:|---:|
| 3 | 120.676 | 322.68% | 60.16% |
| 4 | 27.126 | 65.23% | 52.47% |
| 6 | 35.829 | 57.83% | 40.78% |
| 8 | 2.004 | 4.45% | 4.44% |
| 12 | 1.376 | 3.65% | 3.65% |
| 16 | 1.806 | 3.03% | 3.03% |

The condition ratio need not be monotone in sensor count because the spatial layout changes. These outcomes are from particular sensor subsets and trial distributions, not an optimal-design claim.

![Figure 7.1. Task 5 robustness: noise, sensor count, regularization sensitivity, and conditioning.](figures/phase5_robustness_experiments.png){width=88%}

## 7.4 Ill-conditioned and underdetermined 25-coefficient inverse systems

Notebook 6 uses two response matrices, $A\in\mathbb R^{48\times25}$ and $A_{\mathrm{sparse}}\in\mathbb R^{18\times25}$. The full matrix has rank 25 and condition number about 1,025; the reduced matrix has rank 18 and nullity 7. Over 100 noisy trials the reported mean **source-field errors** were approximately 435.34% (unregularized) versus 31.23% (ridge) for 48 observations, and 2168.42% versus 48.17% for 18 observations. In such regimes, quadratic regularization stabilizes coefficient estimation while leaving irretrievable null-space information unresolved.

## 7.5 Source-model mismatch and grid variation (Notebook 7)

For inversion-grid subdivisions $N\in\{18,24,36,48\}$, reference observations use a finer 72-subdivision forward model. The cross-grid setup includes interpolation; changes in resolution also alter time discretization. Source-model mismatch adds a localized feature outside the assumed basis. The saved JSON gives:

| Grid subdivisions | Matched in-grid clean field error | Cross-grid noisy field error | Mismatched noisy field error |
|---:|---:|---:|---:|
| 18 | 4.06% | 4.17% | 7.92% |
| 24 | 4.13% | 4.19% | 8.08% |
| 36 | 4.13% | 4.16% | 8.07% |
| 48 | 4.14% | 4.16% | 8.07% |

The model-mismatch field error remains near 8% across the tested grid range, whereas matched cases are near 4%. This does **not** prove that all mesh refinement is ineffective; it shows that the tested grid changes do not eliminate this particular representation error. In the 36-subdivision case, observation residuals stay close to 1.6–1.7% despite a much larger difference between the source errors.

![Figure 7.2. Task 7 model mismatch: true and recovered source fields, differences, grid sensitivity, and observation residuals.](figures/phase7_discretization_model_mismatch.png){width=88%}

## 7.6 Discrete adjoint validation and reconstruction (Notebook 8)

The 24-subdivision, 64-step configuration has $\Delta t=0.0015625$ and total stability quantity $4r=0.36$. It uses 48 observations. The gradient errors are $5.5158\times10^{-11}$ against central finite differences and $7.9804\times10^{-16}$ against the explicit matrix gradient. The 25-parameter optimization requires 179 iterations and returns source-field error $9.4302\%$; the 100-parameter optimization requires 211 iterations. Both signal successful numerical optimization, not guaranteed unique physical recovery.

![Figure 7.3. Task 8 adjoint verification and reconstruction; the cost panel compares *one gradient* with *constructing the response matrix*.](figures/phase8_adjoint_results.png){width=88%}

## 7.7 Independent convergence and random sensor layouts (Notebook 9)

The spatial experiment compares an **exact semi-discrete sinusoidal amplitude** to the exact continuous amplitude. It is not an error estimate of a fixed-2000-step FTCS evolution. It gives:

| Spatial subdivisions | Exact semi-discrete vs continuous relative error | Observed spatial order |
|---:|---:|---:|
| 12 | 0.00054413 | — |
| 24 | 0.00013623 | 1.99793 |
| 48 | 0.00003407 | 1.99948 |

The temporal experiment fixes $N=16$ and compares FTCS to the exact semi-discrete amplitude. The unstable 16-step case was removed; all reported configurations satisfy $r_x+r_y\le0.5$:

| Time steps | $\Delta t$ | Relative time-discretization error | Observed temporal order | Stability $r_x+r_y$ |
|---:|---:|---:|---:|---:|
| 32 | 0.003125000 | 0.00279242 | — | 0.32 |
| 64 | 0.001562500 | 0.00139355 | 1.00275 | 0.16 |
| 128 | 0.000781250 | 0.00069611 | 1.00137 | 0.08 |
| 256 | 0.000390625 | 0.00034789 | 1.00069 | 0.04 |

Over 90 randomized layouts per sensor count, the mean source-field error is 42.15% for six sensors versus 20.96% for sixteen sensors; medians are 40.25% and 19.67%. Six sensors at three times yield an 18-by-25 rank-deficient inverse system; sixteen sensors give a full-column-rank 48-by-25 response matrix in the recorded layouts. The reported median ratio of largest to smallest *positive* singular values is about 5006 for the sparse configurations and 1017 for the full ones.

## 7.8 End-to-end benchmark (Notebook 9)

For each $p\in\{25,100\}$, the benchmark compares (i) build-response-matrix-plus-direct-regularized-solve with (ii) L-BFGS-B using the discrete adjoint gradient. It reports medians of three sequential elapsed-time measurements **in the runtime of that session**:

| Parameters | Full matrix solve, median | Full adjoint optimization, median | Relative coefficient gap |
|---:|---:|---:|---:|
| 25 | 0.126 s | 1.582 s | 3.129e-04 |
| 100 | 0.445 s | 2.810 s | 1.266e-04 |

Thus the complete matrix-based method is faster in these small linear trials. Different Task 8 timings and the cost of a *single* gradient must not be conflated with these benchmark totals. CPU model, SciPy and Matplotlib versions were not retained in the saved numerical outputs; the hardware-independent scaling threshold remains unestablished.

![Figure 7.4. Task 9 corrected convergence experiments, randomized-sensor error distributions, and complete-method timing.](figures/task9_validation_figure.png){width=88%}

# 8. Discussion and Limitations

## 8.1 Forward accuracy is not inversion accuracy

The FTCS solver reproduces expected orders for a smooth benchmark and passes superposition tests. Yet its correct implementation is only one component of successful source recovery. Diffusion reduces sensitivity to spatial features, especially when sensors are sparse. The distinction is visible in the model-mismatch test: a small observation residual does not certify spatial source accuracy. Source reconstruction should therefore be evaluated using source-field diagnostics when a known synthetic truth exists, and via independent observations or uncertainty estimates when it does not.

## 8.2 Rank, conditioning, and measurement design

The well-conditioned 3-coefficient model recovers its known source accurately. Expanding to 25 parameters produces weakly observed directions; reducing observations to 18 necessarily introduces a seven-dimensional null space. Rank deficiency is a failure of uniqueness, while a large finite condition number indicates sensitivity along observed directions. Random sensor experiments show that *placement* matters as well as count. These experiments motivate, but do not solve, an optimal sensor-placement problem.

## 8.3 Regularization and basis dependence

Tikhonov regularization suppresses amplification along small singular values. It yields substantial benefit in the noisy 25-coefficient experiment but introduces a mild bias in the well-conditioned noiseless 3-coefficient problem. Coefficient-norm penalties depend on the scaling and overlap of basis functions. A larger basis gives representational freedom while also increasing potential ill-conditioning. Future studies should compare coefficient penalties with source-field smoothness penalties and explore regularization-parameter selection without using the unavailable true source.

## 8.4 Mismatch and interpretation of convergence

The Task 7 mismatch case contains a source component outside the Gaussian dictionary. The observed residual field error remains across tested grids because numerical discretization and model representation address different questions. The Task 9 semi-discrete spatial reference is especially valuable here: it shows how to isolate a numerical error component rather than attribute all differences to the PDE solver. A separate asymptotic theory for the inverse reconstruction under mesh refinement is outside the scope of this study.

## 8.5 Adjoint correctness versus computational efficiency

The adjoint gradient is validated to near machine precision against the explicit matrix gradient, and the optimizer agrees closely with a direct minimizer. Nevertheless, the *complete* 25- and 100-parameter matrix workflows finish sooner in the recorded environment. Avoiding an explicit sensitivity matrix is a structural advantage that may become relevant at larger scales, but the experiment does not establish a speedup crossover or general hardware-independent superiority.

## 8.6 Boundaries of inference

The study uses synthetic data, a unit-square domain, constant diffusivity, zero initial and Dirichlet boundary values, time-independent sources, fixed Gaussian dictionaries, and a modest set of source dimensions. Random trials quantify empirical variation for specified noise models; they do not provide calibrated uncertainty bands for real instruments. The 90 randomized layouts are finite samples, not a comprehensive geometry search. Task 9 wall-clock timings depend on execution environment. Finally, the current archival package was inspected statically and had been run in Colab by the user, but a new full clean-kernel run under a locked environment has not been performed as part of this manuscript consolidation.

# 9. Conclusions and Future Work

The project derived and verified a two-dimensional FTCS heat solver, connected Gaussian source coefficients to sparse observations, and tested least-squares, Tikhonov, and discrete-adjoint reconstruction. Independent semi-discrete spatial comparisons produced orders close to 2, while the corrected stable temporal experiment produced orders close to 1. In a three-parameter matched-basis problem, source-field reconstruction error was 0.5684%; larger, poorly conditioned and rank-deficient problems were much more noise sensitive. Regularization substantially reduced errors in the large inverse tests but could not supply measurements that were never collected. The mismatch experiment demonstrated a persistent error component from an inadequate dictionary across the tested grid resolutions. The discrete adjoint gradient was accurately reproduced by both finite differences and the explicit matrix gradient, although direct matrix reconstruction was faster in the small tested workflows.

The broader conclusion is that reliable source reconstruction requires simultaneous attention to the PDE approximation, parameterization, observation geometry, regularization, and optimization. A numerically accurate forward solver, low observation residual, or accurate optimization gradient cannot independently guarantee an accurate recovered source.

Potential extensions include adaptive or multiscale basis dictionaries, sensor-location and observation-time optimization, spatially informed regularization and parameter-selection rules, uncertainty quantification, heterogeneous diffusivity and more general boundary conditions, and larger-scale matrix-free PDE-constrained optimization. These are research directions suggested by the computational results, not established improvements of the present work.

# References

[1] LeVeque, R. J. (2007). *Finite Difference Methods for Ordinary and Partial Differential Equations: Steady-State and Time-Dependent Problems*. SIAM. https://doi.org/10.1137/1.9780898717839

[2] Engl, H. W., Hanke, M., & Neubauer, A. (1996). *Regularization of Inverse Problems*. Kluwer Academic Publishers.

[3] Hansen, P. C. (2010). *Discrete Inverse Problems: Insight and Algorithms*. SIAM. https://doi.org/10.1137/1.9780898718836

[4] Tikhonov, A. N., & Arsenin, V. Y. (1977). *Solutions of Ill-Posed Problems*. V. H. Winston & Sons.

[5] Vogel, C. R. (2002). *Computational Methods for Inverse Problems*. SIAM. https://doi.org/10.1137/1.9780898717570

[6] Lions, J.-L. (1971). *Optimal Control of Systems Governed by Partial Differential Equations*. Springer.

[7] Hinze, M., Pinnau, R., Ulbrich, M., & Ulbrich, S. (2009). *Optimization with PDE Constraints*. Springer. https://doi.org/10.1007/978-1-4020-8839-1

[8] Giles, M. B., & Pierce, N. A. (2000). An introduction to the adjoint approach to design. *Flow, Turbulence and Combustion*, 65, 393–415. https://doi.org/10.1023/A:1011430410075

[9] Adavani, S. S., & Biros, G. (2008). Multigrid algorithms for inverse problems with linear parabolic PDE constraints. *SIAM Journal on Scientific Computing*, 31(1), 369–397. https://doi.org/10.1137/070687426

[10] Nocedal, J., & Wright, S. J. (2006). *Numerical Optimization* (2nd ed.). Springer. https://doi.org/10.1007/978-0-387-40065-5

# Appendix A. Reproducibility and Data Provenance

## A.1 Project structure

The companion package includes nine numbered notebooks in `notebooks/`, datasets and machine-readable metrics in `data/`, four submitted PNG figures in `figures/`, and README/dependency documentation. The manuscript is based on the corrected Notebook 9, not the older initial archive. Tasks 3–9 export numeric and JSON artifacts; not every notebook independently exports an NPZ file.

## A.2 Run order and dependencies

Run notebooks in ascending order from 01 to 09 in Python with NumPy, SciPy, Matplotlib, and Jupyter/Colab. Notebook 4 uses the saved Task 3 observation NPZ. Preserve the random seeds and noise-generation formulas from each notebook. The accompanying `requirements.txt` names the packages but is not a dependency lockfile; the Task 9 summary records NumPy version 2.1.3 and seed 20261018, while complete CPU and package metadata were not preserved. Benchmarks should be rerun and hardware documented if cross-platform performance is an intended claim.

## A.3 Artifact handling

Load only numeric keys from `phase8_results.npz` using default `np.load(..., allow_pickle=False)`. The `metrics` field is an object array: use the separate plain JSON `phase8_metrics.json` instead of enabling pickle for an untrusted archive. The four figures are: `phase5_robustness_experiments.png`, `phase7_discretization_model_mismatch.png`, `phase8_adjoint_results.png`, and `task9_validation_figure.png`. Corrected Task 9 exports contain temporal steps 32, 64, 128, and 256 exclusively.

## A.4 Verification status

The nine notebooks and associated saved artifacts were inspected for source consistency and data availability; the user previously executed them in Colab. A clean rerun of the assembled archive was not part of this manuscript-production pass. Numerical results remain experiment-specific; a future reproduction should archive Python, NumPy, SciPy, Matplotlib, platform, CPU, runtime configuration, and any package versions used for timing.

## A.5 Verification checklist for future readers

A reproducibility run should: verify that the forward temperature remains zero under zero source and homogeneous initial/boundary conditions; confirm the timestep satisfies $r_x+r_y\le1/2$; check a sinusoidal source against the analytic solution; test superposition by comparing the direct mixed-source run with $A\theta$; verify the original three-coefficient observation matrix has shape $48\times3$; check that the corrected temporal study includes only $32,64,128,256$ steps; load source and prediction arrays with the correct definition of relative error; verify the adjoint gradient against both central differences and the explicit matrix gradient; and compare complete reconstruction timings without equating one gradient evaluation with an entire optimization run.

For a future public archive, record exact Python, NumPy, SciPy and Matplotlib versions, runtime platform, machine type, random generator and seed, notebook execution order, and checksums of the archived datasets. The present files support the stated numerical findings, but the most recent combined package has not been certified by a single clean, scripted end-to-end rerun.

# Appendix B. Symbols and Notation

| Symbol | Meaning |
|---|---|
| $u$, $U_{ij}^n$ | Continuous temperature and its discrete approximation |
| $f$, $F$ | Continuous and sampled spatial heat source |
| $\alpha$ | Constant diffusivity |
| $h_x,h_y,\Delta t$ | Spatial mesh widths and timestep |
| $N_x,N_y,N_t$ | Numbers of space intervals and time steps |
| $q,p,m$ | Interior temperature unknowns, source coefficients, scalar observations |
| $B$ | Discrete source basis matrix, $q\times p$ |
| $M$ | Explicit forward heat propagator, $q\times q$ |
| $A$ | Discrete source-to-observation response matrix, $m\times p$ |
| $\theta$ | Unknown source coefficient vector |
| $\lambda$ | Coefficient-norm Tikhonov parameter |
| $p^n$ | Time-indexed discrete adjoint state (distinct from integer $p$) |
| $E_\theta,E_f$ | Relative coefficient error and source-field error |


# Appendix C. Reconstructing the Observation Operator

The numerical response matrix can be constructed without matrix exponentiation. For each $\ell=1,\ldots,p$, create a coefficient vector with entry $\ell$ equal to one and all other entries zero. Compute the corresponding Gaussian source field $B e_\ell$, initialize the interior temperature to zero, and advance the stable FTCS solver. At each chosen observation step, sample the temperature at the prescribed sensor nodes, keeping the same stacking order across all columns. Assemble these observation vectors as the columns of $A$.

Given a reference coefficient vector $\theta^\star$, one may calculate the direct forward observations from the combined source $B\theta^\star$ and independently calculate $A\theta^\star$. Their agreement within floating-point precision checks the source-basis mapping, forward linearity, sensor indexing, and observation stacking. Once a noisy measurement vector $d$ is generated, the reconstruction is a separate operation: determine $\widehat\theta$ by least squares or regularized optimization, then report both $\widehat d=A\widehat\theta$ and the reconstructed spatial source $B\widehat\theta$. Comparing $\widehat d$ with $d$ assesses data fit, whereas comparing $B\widehat\theta$ with a known synthetic source assesses source accuracy.

The column construction is possible because the governing state equation and observation map are linear. For nonlinear PDEs, a fixed response matrix generally no longer represents the complete source-to-data map, although local sensitivity matrices and adjoint gradients may still be used.
