---
title: Richards Equation
description: ''
date: 2017-12-22T00:00:00.000Z
authors:
  - name: Rowan Cockett
    userId: vKndfPAZO7WeFxLH1GQcpnXPzfH3
    orcid: 0000-0002-7859-8394
    corresponding: true
    email: rcockett@eoas.ubc.ca
    roles: null
    affiliations:
      - University of British Columbia
---

+++

```{admonition} Preface
This chapter presents a computationally scalable algorithm for solving inverse problems for hydraulic parameters in vadose zone flow using the Richards equation. This work has been submitted for peer review and the preprint is available on _arXiv_ {cite:p}`Cockett2017`; preliminary versions of this research were presented at two conferences {cite:p}`Cockett2013, Cockett2013a`.
```

---

# Introduction

Studying the processes that occur in the vadose zone, the region between the earth's surface and the fully saturated zone, is of critical importance for understanding our groundwater resources. Fluid flow in the vadose zone is described by the Richards equation and parameterized by hydraulic conductivity, which is a nonlinear function of pressure head {cite:p}`Richards1931, Celia1990`. Typically, hydraulic conductivity is heterogeneous and can have a large dynamic range. In any site characterization, the spatial estimation of the hydraulic conductivity function is an important step. Achieving this, however, requires the ability to efficiently solve and optimize the nonlinear, time-domain Richards equation. Rather than working with a full, implicit, 3D time-domain system of equations, simplifications are consistently used to avert the conceptual, practical, and computational difficulties inherent in the parameterization and inversion of the Richards equation. These simplifications typically parameterize the conductivity and assume that it is a simple function in space, often adopting a homogeneous or one dimensional layered soil profile (cf. {cite:p}`Binley2002, Deiana2007, Hinnell2010, Liang2014`). Due to the lack of constraining hydrologic data, such assumptions are often satisfactory for fitting observed measurements, especially in two and three dimensions as well as in time. However, as more data become available, through spatially extensive surveys and time-lapse proxy measurements (e.g. direct current resistivity surveys and distributed temperature sensing), extracting more information about subsurface hydrogeologic parameters becomes a possibility. The proxy data can be directly incorporated through an empirical relation (e.g. {cite:p}`Archie1942`) or time-lapse estimations can be structurally incorporated through some sort of regularization technique {cite:p}`HaberHoltzman2013, ho, Hinnell2010`. Recent advances have been made for the forward simulation of the Richards equation in a computationally-scalable manner {cite:p}`RichardsFOAM`. However, the inverse problem is non-trivial, especially in three-dimensions {cite:p}`Towara2015`, and must be considered using modern numerical techniques that allow for spatial estimation of hydraulic parameters. However, this is especially intricate to both derive and implement due to the nonlinear, time-dependent forward simulation and potential model dependence in many aspects of the Richards equation (e.g. multiple empirical relations, boundary/initial conditions). To our knowledge, there has been no large-scale inversion for distributed hydraulic parameters in three dimensions using the Richards equation as the forward simulation.

Inverse problems in space and time are often referred to as history matching problems (see {cite:t}`DeanChen2011, OliverBook2008, SarmaDurlofskyAziz2007, Oliver01, hydrusCalibration2012` and reference within). Inversions use a flow simulation model, combined with some a-priori information, in order to estimate a spatially variable hydraulic conductivity function that approximately yields the observed data. The literature shows a variety of approaches for this inverse problem, including trial-and-error, stochastic methods, and various gradient based methods {cite:p}`Bitterlich2004, Binley2002, Carrick2010, Durner1994, Finsterle2011c, Mualem1976, Simunek1996`. The way in which the computational complexity of the inverse method scales becomes important as problem size increases {cite:p}`Towara2015`. Computational memory and time often become a bottleneck for solving the inverse problem, both when the problem is solved in 2D and, particularly, when it is solved in 3D {cite:p}`hao`. To solve the inverse problem, stochastic methods are often employed, which have an advantage in that they can examine the full parameter space and give insights into non-uniqueness {cite:p}`Finsterle2011`. However, as the number of parameters we seek to recover in an inversion increases, these stochastic methods require that the forward problem be solved many times, which often makes these methods impractical. This scalability, especially in the context of hydrogeophysics has been explicitly noted in the literature (cf. {cite:t}`Binley2002, Deiana2007, Towara2015, Linde2016`).

Derivative-based optimization techniques become a practical alternative when the forward problem is computationally expensive or when there are many parameters to estimate (i.e. thousands to millions). Inverse problems are ill-posed and thus to pose a solvable optimization problem, an appropriate regularization is combined with a measure of the data misfit to state a deterministic optimization problem {cite:p}`tikhonov1977`. Alternatively, if prior information can be formulated using a statistical framework, we can use Bayesian techniques to obtain an estimator through the Maximum A Posteriori model (MAP) {cite:p}`somersallo`. In the context of Bayesian estimation, gradient based methods are also important, as they can be used to efficiently sample the posterior {cite:p}`BuiThanhGhattas2015`.

A number of authors have sought solutions for the inverse problem, where the forward problem is the Richards equation (cf. {cite:p}`Bitterlich2002, Iden2007, hydrusCalibration2012` and references within). The discretization of the Richards equation is commonly completed by an implicit method in time and a finite volume or finite element method in space. Most work uses a Newton-like method for the resulting nonlinear system, which arises from the discretization of the forward problem. For the deterministic inverse problem using the Richards equation, previous work uses some version of a Gauss-Newton method (e.g. Levenberg-Marquardt), with a **direct** calculation of the sensitivity matrix {cite:p}`Finsterle2011, Simunek1996, Bitterlich2002`. However, while these approaches allow for inversions of moderate scale, they have one major drawback: the sensitivity matrix is large and dense; its computation requires dense linear algebra and a non-trivial amount of memory (cf. {cite:p}`Towara2015`). Previous work used either external numerical differentiation (e.g. PEST) or automatic differentiation in order to directly compute the sensitivity matrix {cite:p}`Finsterle2011c, Bitterlich2002, Doherty2004, Towara2015`. Finite difference can generate inaccuracies in the sensitivity matrix and, consequently, tarry the convergence of the optimization algorithm. Furthermore, external numerical differentiation is computationally intensive and limits the number of model parameters that can be estimated.

The goal of this chapter is to suggest a modern numerical formulation that allows the inverse problem to be solved **without explicit** computation of the sensitivity matrix by using **exact** derivatives of the discrete formulation {cite:p}`hao`. Our technique is based on the discretize-then-optimize approach, which discretizes the forward problem first and then uses a deterministic optimization algorithm to solve the inverse problem {cite:p}`Gunz03`. To this end, we require the discretization of the forward problem. Similar to the work of {cite:p}`Celia1990`, we use an implicit Euler method in time and finite volume in space. Given the discrete form, we show that we can analytically compute the derivatives of the forward problem with respect to distributed hydraulic parameters and, as a result, obtain an implicit formula for the sensitivity. The formula involves the solution of a linear time-dependent problem; we avoid computing and storing the sensitivity matrix directly and, rather, suggest a method to efficiently compute the product of the sensitivity matrix and its adjoint times a vector. Equipped with this formulation, we can use a standard inexact Gauss-Newton method to solve the inverse problem for distributed hydraulic parameters in 3D. This large-scale distributed parameter estimation becomes computationally tractable with the technique presented in this chapter and can be employed with any iterative Gauss-Newton-like optimization technique.

This chapter is structured as follows: in Section \ref{sec:richards-forward}, we discuss the discretization of the forward problem on a staggered mesh in space and backward Euler in time; in Section \ref{sec:richards-inverse}, we formulate the inverse problem and construct the implicit functions used for computations of the Jacobian-vector product. In Section \ref{sec:richards-validation}, we demonstrate the validity of the implementation of the forward problem and sensitivity calculation. In Section \ref{sec:richards-results}, we validate the numerical implementation and compare to the literature. Chapter \ref{ch:applications} will expand upon the techniques introduced in this chapter to show the effectiveness of the implicit sensitivity algorithm in comparison to existing numerical techniques.

## Attribution and dissemination

To accelerate both the development and dissemination of this approach, we have built these tools on top of an open source framework for organizing simulation and inverse problems in geophysics {cite:p}`simpeg2015`. We have released our numerical implementation under the permissive MIT license. Our implementation of the implicit sensitivity calculation for the Richards equation and associated inversion implementation is provided and tested to support 1D, 2D, and 3D forward and inverse simulations with respect to custom empirical relations and sensitivity to parameters within these functions. The source code can be found at <https://github.com/simpeg/simpeg> and may be a helpful resource for researchers looking to use or extend our implementation. I have presented early versions of this work at two international conferences {cite:p}`Cockett2013a, Cockett2013` and have submitted a version of this manuscript for peer review {cite:p}`Cockett2017`.

% (sec:richards-forward)=

# Forward problem

In this section, we describe the Richards equations and its discretization {cite:p}`Richards1931`. The Richards equation is a nonlinear parabolic partial differential equation (PDE) and we follow the so-called mixed formulation presented in {cite:p}`Celia1990` with some modifications. In the derivation of the discretization, we give special attention to the details used to efficiently calculate the effect of the sensitivity on a vector, which is needed in any derivative based optimization algorithm.

## Richards equation

The parameters that control groundwater flow depend on the effective saturation of the media, which leads to a nonlinear problem. The groundwater flow equation has a diffusion term and an advection term which is related to gravity and only acts in the $z$-direction. There are two different forms of the Richards equation; they differ in how they deal with the nonlinearity in the time-stepping term. Here, we use the most fundamental form, referred to as the 'mixed'-form of the Richards equation {cite:p}`Celia1990`:

% equations/richards/richards-mixed

```{math}
:label: eq:richards-mixed
\frac{\partial \theta(\psi)}{\partial t} - \nabla \cdot k(\psi) \nabla \psi - \frac{\partial k(\psi)}{\partial z} = 0
\quad \psi \in \Omega
```

where $\psi$ is pressure head, $\theta(\psi)$ is volumetric water content, and $k(\psi)$ is hydraulic conductivity. This formulation of the Richards equation is called the 'mixed'-form because the equation is parameterized in $\psi$ but the time-stepping is in terms of $\theta$. The hydraulic conductivity, $k(\psi)$, is a heterogeneous and potentially anisotropic function that is assumed to be known when solving the forward problem. In this chapter, we assume that $k$ is isotropic, but the extension to anisotropy is straightforward {cite:p}`simpeg2015, fvtutorial`. The equation is solved in a domain, $\Omega$, equipped with boundary conditions on $\partial \Omega$ and initial conditions, which are problem-dependent.

% For this study, we assume the so-called "natural" or no-flux boundary conditions on $\psi$
%% TODO: Why introduce BC here..?
% $$ \GRAD \psi \cdot {\vec n} = 0 \quad \psi \in \partial \Omega $$

An important aspect of unsaturated flow is noticing that both water content, $\theta$, and hydraulic conductivity, $k$, are functions of pressure head, $\psi$. There are many empirical relations used to relate these parameters, including the Brooks-Corey model {cite:p}`Brooks1964` and the van Genuchten-Mualem model {cite:p}`Mualem1976, VanGenuchten1980`. The van Genuchten model is written as:

% equations/richards/van-genuchten

% \begin{subequations}
% \label{eq:van-genuchten}

% :label: eq:van-genuchten-water-retention
% :label: eq:van-genuchten-hydraulic-conductivity

```{math}
:label: eq:van-genuchten
\begin{align*}
  \theta(\psi) &=
  \left\{\begin{aligned}
      \theta_r& + \frac{\theta_s- \theta_r}{(1+|\alpha \psi|^n)^m}  & \psi < 0 \\
      \theta_s& & \psi \ge 0
  \end{aligned}\right.
  \\
  k(\psi) &=
  \left\{\begin{aligned}
      K_s & \theta_e(\psi)^l(1-(1- \theta_e(\psi)^{-m})^m)^2 & \psi < 0 \\
      K_s& & \psi \ge 0
  \end{aligned}\right.
\end{align*}
```

% \end{subequations}

where

```{math}
:label: eq:van-genuchten-params
\theta_e(\psi) = \frac{\theta(\psi) - \theta_r}{\theta_s - \theta_r},
\qquad
m=1- \frac{1}{n},
\qquad
n > 1
```

Here, $\theta_r$ and $\theta_s$ are the residual and saturated water contents, $K_s$ is the saturated hydraulic conductivity, $\alpha$ and $n$ are fitting parameters, and, $\theta_e(\psi) \in [0,1]$ is the effective saturation. The pore connectivity parameter, $l$, is often taken to be $\frac{1}{2}$, as determined by {cite:t}`Mualem1976`. {numref}`Figure %s <fig:van-genuchten>` shows the functions over a range of negative pressure head values for four soil types (sand, loam, sandy clay, and clay). The pressure head varies over the domain $\psi \in (-\infty, 0)$. When the value is close to zero (the left hand side), the soil behaves most like a saturated soil where $\theta = \theta_s$ and $k = K_s$. As the pressure head becomes more negative, the soil begins to dry, which the water retention curve shows as the function moving towards the residual water content ($\theta_r$). Small changes in pressure head can change the hydraulic conductivity by several orders of magnitude; as such, $k(\psi)$ is a highly nonlinear function, making the Richards equation a nonlinear PDE.

```{figure} images/van-genuchten.png
:name: fig:van-genuchten
The water retention curve and the hydraulic conductivity function for four canonical soil types of sand, loam, sandy clay, and clay.
```

## Discretization

The Richards equation is parameterized in terms of pressure head, $\psi$. Here, we describe simulating the Richards equation in one, two, and three dimensions. We start by discretizing in space and then we discretize in time. This process yields a discrete, nonlinear system of equations; for its solution, we discuss a variation of Newton's method.

### Spatial Discretization

In order to conservatively discretize the Richards equation, we introduce the flux ${\vec f}$ and rewrite the equation as a first order system of the form:

% equations/richards/richards-mixed-first-order
% \begin{subequations}

% \label{eq:richards-mixed-first-order-a}
% \label{eq:richards-mixed-first-order-b}

```{math}
:label: eq:richards-mixed-first-order
\begin{align*}
\frac{\partial \theta(\psi)}{\partial t} - \nabla \cdot {\vec f} - \frac{\partial k(\psi)}{\partial z} & = 0 \\
k(\psi)^{-1} {\vec f} & =  \nabla \psi
\end{align*}
```

% \end{subequations}

We then discretize the system using a standard staggered finite volume discretization (cf. {cite:t}`Ascher2008, haber2015computational, fvtutorial`, and Appendix \ref{ch:discretize}). This discretization is a natural extension of mass-conservation in a volume where the balance of fluxes into and out of a volume are conserved {cite:p}`LipnikovMisitas2013`. Here, it is natural to assign the entire cell one hydraulic conductivity value, $k$, which is located at the cell center. Such assigning leads to a piecewise constant approximation for the hydraulic conductivity and allows for discontinuities between adjacent cells. From a geologic perspective, discontinuities are prevalent, as it is possible to have large differences in hydraulic properties between geologic layers in the ground. The pressure head, $\psi$, is also located at the cell centers and the fluxes are located on cell faces, which lead to the usual staggered mesh or Marker and Cell (MAC) discretization in space {cite:p}`fletcher`. We demonstrate the discretization in 1D, 2D and 3D on the tensor mesh in {numref}`Figure %s <fig:richards-finite-volume>`. We discretize the function, $\psi$, on a cell-centered grid, which results in a grid function, $\bfpsi$. We use bold letters to indicate other grid functions.

```{figure} images/finite-volume.png
:name: fig:richards-finite-volume
Discretization of unknowns in 1D, 2D and 3D space. Red circles are the locations of the discrete hydraulic conductivity $K$ and the pressure head $\psi$.
The arrows are the locations of the discretized flux $\vec f$ on each cell face. Modified after {cite:t}`fvtutorial`.
```

The discretization of a diffusion-like equation on an orthogonal mesh is well-known (see {cite:p}`Haber2001,fletcher,HaberHeldmannAscher07,Ascher2011` and reference within). We discretize the differential operators by using the usual mass balance consideration and the elimination of the flux, $\bff$ [^2]. This spatial discretization leads to the following discrete nonlinear system of ordinary differential equations (assuming homogeneous Dirichlet boundary conditions):

[^2]: Here we assume an isotropic conductivity that leads to a diagonal mass matrix and this yields easy elimination of the fluxes.

% equations/richards/richards-mixed-discrete

```{math}
:label: eq:richards-mixed:discrete
\frac{d \boldsymbol{\theta}(\boldsymbol{\psi})}{d t}
- \mathbf{D}
    \text{ diag}
    \left(
        \mathbf{k}_{Av}(\boldsymbol{\psi}^{n+1})
    \right)
\mathbf{G} \boldsymbol{\psi}
- \mathbf{G}_z
    \left(
        \mathbf{k}_{Av}(\boldsymbol{\psi}^{n+1})
    \right)
=0
```

Here, $\bfD$ is the discrete divergence operator and $\mathbf{G}$ is the discrete gradient operator. The discrete derivative in the $z$-direction is written as $\mathbf{G}_z$. The values of $\psi$ and $k(\psi)$ are known on the cell-centers and must be averaged to the cell-faces, which we complete through harmonic averaging {cite:p}`Haber2001`.

% equations/richards/mass-matrix-kav

```{math}
:label: eq:mass-matrix-kav
\mathbf{k}_{Av}(\boldsymbol{\psi}) = \frac{1}{\mathbf{A}_v\frac{1}{\mathbf{k}(\boldsymbol{\psi})}}
```

where $\bfA_{v}$ is a matrix that averages from cell-centers to faces and the division of the vector is done pointwise; that is, we use the vector notation, $(1/{\bfv})_{i} = 1/\bfv_{i}$. We incorporate boundary conditions using a ghost-point outside of the mesh {cite:p}`tos`.

### Time discretization and stepping

The Richards equation is often used to simulate water infiltrating an initially dry soil. At early times in an infiltration experiment, the pressure head, $\psi$, can be close to discontinuous. These large changes in $\psi$ are also reflected in the nonlinear terms $k(\psi)$ and $\theta(\psi)$; as such, the initial conditions imposed require that an appropriate time discretization be chosen. Hydrogeologists are often interested in the complete evolutionary process, until steady-state is achieved, which may take many time-steps. Here, we describe the implementation of a fully-implicit backward Euler numerical scheme. Higher-order implicit methods are not considered here because the uncertainty associated with boundary conditions and the fitting parameters in the Van Genuchten models (eq. {eq}`eq:van-genuchten`) have much more effect than the order of the numerical method used.

The discretized approximation to the mixed-form of the Richards equation, using fully-implicit backward Euler, reads:

% equations/richards/richards-mixed-discrete-time

```{math}
:label: eq:richards-mixed-discrete-time
    F(\boldsymbol{\psi}^{n+1},\boldsymbol{\psi}^n) = \frac{
    \boldsymbol{\theta}(\boldsymbol{\psi}^{n+1}) - \boldsymbol{\theta}(\boldsymbol{\psi}^n)
    }{\Delta t}
    - \mathbf{D}
        \text{ diag}
        \left(
            \mathbf{k}_{Av}(\boldsymbol{\psi}^{n+1})
        \right)
    \mathbf{G} \boldsymbol{\psi}^{n+1}
    - \mathbf{G}_z
        \left(
            \mathbf{k}_{Av}(\boldsymbol{\psi}^{n+1})
        \right)
    =
    0
```

This is a nonlinear system of equations for $\bfpsi\nn$ that needs to be solved numerically by some iterative process. Either a Picard iteration (as in {cite:t}`Celia1990`) or a Newton root-finding iteration with a step length control can be used to solve the system. Note that to deal with dependence of $\theta$ with respect to $\psi$ in Newton's method, we require the computation of $\deriv{\bftheta}{\bfpsi}$. We can complete this computation by using the analytic form of the hydraulic conductivity and water content functions (e.g. derivatives of eq. {eq}`eq:van-genuchten`). We note that a similar approach can be used for any smooth curve, even when the connection between $\theta$ and $\psi$ are determined empirically (for example, when $\theta(\psi)$ is given by a spline interpolation of field data).

## Solving the nonlinear equations

Regardless of the empirical relation chosen, we must solve {eq}`eq:richards-mixed-discrete-time` using an iterative root-finding technique. Newton's method iterates over $m=1,2,\dots$ until a satisfactory estimation of $\bfpsi^{n+1}$ is obtained. Given $\bfpsi\nnm$, we approximate $\FF(\bfpsi^{n+1},\bfpsi\n)$ as:

% equations/richards/richards-newton

```{math}
:label: eq:richards-newton
F(\boldsymbol{\psi}^{n+1},\boldsymbol{\psi}^n)
\approx
F(\boldsymbol{\psi}^{n,m},\boldsymbol{\psi}^n) +
\mathbf{J}_{\boldsymbol{\psi}^{n+1,m}} \delta \boldsymbol{\psi}
```

where the Jacobian for iteration, $m$, is:

% equations/richards/richards-newton-deriv-setup

```{math}
:label: eq:richards-newton-deriv-setup
\mathbf{J}_{\psi^{n+1,m}}  = \left. {\frac {\partial F(\boldsymbol{\psi},\boldsymbol{\psi}^n)}{\partial \boldsymbol{\psi}}} \right|_{\boldsymbol{\psi}^{n+1,m}}
```

The Jacobian is a large dense matrix, and its computation necessitates the computation of the derivatives of $\FF(\bfpsi\nnm,\bfpsi\n)$. We can use numerical differentiation in order to evaluate the Jacobian (or its product with a vector). However, in the context of the inverse problem, an exact expression is preferred. Given the discrete forward problem, we obtain that:

% equations/richards/richards-newton-deriv

```{math}
:label: eq:richards-newton-deriv
    \bfJ_{ \bfpsi\nnm} =
    \frac{1}{\Delta t} \deriv{\bftheta(\bfpsi\nnm)}{\bfpsi\nnm}
    -
    \deriv{}{\bfpsi\nnm}
    \left(
    \bfD
        \diag{\bfk_{Av}(\bfpsi\nnm)}
    \mathbf{G}\bfpsi\nnm
    \right)
    -
    \mathbf{G}_z
    \deriv{\bfk_{Av}(\bfpsi\nnm)}{\bfpsi\nnm}
```

Here, recall that $\bfk_{Av}$ is harmonically averaged and its derivative can be obtained by the chain rule:

% equations/richards/mass-matrix-kav-deriv

```{math}
:label: eq:mass-matrix-kav-deriv
\deriv{\mathbf{k}_{Av}(\boldsymbol{\psi})}{\boldsymbol{\psi}}
=
\text{diag}\left(
    (\mathbf{A}_v \mathbf{k}^{-1}(\boldsymbol{\psi}))^{-2} \right)
\mathbf{A}_v
\text{diag} \left(
    \mathbf{k}^{-2}(\boldsymbol{\psi}) \right)
\deriv{\mathbf{k}(\boldsymbol{\psi})}{\boldsymbol{\psi}}
```

Similarly, for the second term in ({eq}`eq:richards-newton-deriv`) we obtain:

% equations/richards/richards-newton-deriv-detail

```{math}
:label: eq:richards-newton-deriv-detail
\frac{\partial}{\partial\boldsymbol{\psi}}
\left(
\mathbf{D}
    \text{ diag}\left(\mathbf{k}_{Av}(\boldsymbol{\psi})\right)
\mathbf{G} \boldsymbol{\psi}
\right)
=
\mathbf{D}
    \text{ diag}\left(\mathbf{k}_{Av}(\boldsymbol{\psi})\right)
\mathbf{G}
+
\mathbf{D}
    \text{ diag}\left(\mathbf{G}\boldsymbol{\psi}\right)
\frac{\partial\mathbf{k}_{Av}(\boldsymbol{\psi})}{\partial\boldsymbol{\psi}}
```

Here the notation $\nnm$ has been dropped for brevity. For the computations above, we need the derivatives of functions $\bfk(\bfpsi)$ and $\bftheta(\bfpsi)$; note that, since the relations are assumed local (point wise in space) given the vector, $\bfpsi$, these derivatives are diagonal matrices. For Newton's method, we solve the linear system:

% equations/richards/richards-newton-update

```{math}
:label: eq:richards-newton-update
\mathbf{J}_{\psi^{n+1,m}}\, \delta \boldsymbol{\psi}\ =\ -  F(\boldsymbol{\psi}^{n+1,m},\boldsymbol{\psi}^n)
```

For small-scale problems, we can solve the linear system using direct methods; however, for large-scale problems, iterative methods are more commonly used. The existence of an advection term in the PDE results in a non-symmetric linear system in the Newton solve. Thus, when using iterative techniques to solve this system, an appropriate iterative method, such as {sc}`bicgstab` or {sc}`gmres` {cite:p}`saad,templates`, must be used. For a discussion on solver choices in the context of the Richards equation please see {cite:t}`RichardsFOAM`.

At this point, it is interesting to note the difference between the Newton iteration and the Picard iteration suggested in {cite:p}`Celia1990`. We can verify that the Picard iteration uses an approximation to the Jacobian $\bfJ_{\psi\nnm}\, \delta \bfpsi$ given by dropping the second term from {eq}`eq:richards-newton-deriv-detail`. This term can have negative eigenvalues and dropping it is typically done when considering the lagged diffusivity method {cite:p}`vogelBook`. However, as discussed in {cite:p}`vogelBook`, ignoring this term can slow convergence.

Finally, a new iterate is computed by adding the Newton update to the last iterate:

$$\bfpsi\nnmm = \bfpsi\nnm+ \alpha \delta \psi$$
where $\alpha$ is a parameter that guarantees that

$$ \|\FF(\bfpsi^{n,m+1},\bfpsi\n) \| < \| \FF(\bfpsi^{n,m},\bfpsi\n) \| $$
To obtain $\alpha$, we perform an Armijo line search {cite:p}`nw`. In our numerical experiments, we have found that this method can fail when the hydraulic conductivity is strongly discontinuous and changes rapidly. In such cases, Newton's method yields a poor descent direction. Therefore, if the Newton iteration fails to converge to a solution, the update is performed with the mixed-form Picard iteration. Note that Picard iteration can be used, even when Newton’s method fails, because Picard iteration always yields a descent direction {cite:p}`vogelBook`.

At this point, we have discretized the Richards equation in both time and space while devoting special attention to the derivatives necessary in Newton's method and the Picard iteration as described in {cite:p}`Celia1990`. The exact derivatives of the discrete problem will be used in the following two sections, which outline the implicit formula for the sensitivity and its incorporation into a general inversion algorithm. The implementation is provided as a part of the open source {sc}`SimPEG` project ({cite:t}`simpeg2015`, <https://simpeg.xyz>).

% (sec:richards-inverse)=

# Inverse Problem

The location and spatial variability of, for example, an infiltration front over time is inherently dependent on the hydraulic properties of the soil column. As such, direct or proxy measurements of the water content or pressure head at various depths along a soil profile contain information about the soil properties. We pose the inverse problem, which is the estimation of distributed hydraulic parameters, given either water content or pressure data. We frame this problem under the assumption that we wish to estimate hundreds of thousands to millions of distributed model parameters. Due to the large number of model parameters that we aim to estimate in this inverse problem, Bayesian techniques or external numerical differentiation, such as the popular PEST toolbox {cite:p}`Doherty2004`, are not computationally feasible. Instead, we will employ a direct method by calculating the exact derivatives of the discrete the Richards equation and solving the sensitivity implicitly. For brevity, we show the derivation of the sensitivity for an inversion model of only saturated hydraulic conductivity, $K_s$, from pressure head data, $\bfdo$. This derivation can be readily extended to include the use of water content data and inverted for other distributed parameters in the heterogeneous hydraulic conductivity function. We will demonstrate the sensitivity calculation for multiple distributed parameters in the numerical examples (Section \ref{sec:richards-examples}).

The Richards equation simulation produces a pressure head field at all points in space as well as through time. Data can be predicted, $\bfdp$, from these fields and compared to observed data, $\bfdo$. To be more specific, we let $\bfPsi = [(\bfpsi^{1})^{\top},\ldots,(\bfpsi^{n_{t}})^{\top}]^{\top}$ be the (discrete) pressure field for all space and $n_{t}$ time steps. When measuring pressure head recorded only in specific locations and times, we define the predicted data, $\bfdp$, as $\bfdp = \bfP \bfPsi(\bfm)$. Here, the vector $\mathbf{m}$ is the vector containing all of the parameters which we are inverting for (e.g. $K_s, \alpha, n, \theta_r$, or $\theta_s$ when using the van Genuchten empirical relation). The matrix, $\bfP$, interpolates the pressure head field, $\bfPsi$, to the locations and times of the measurements. Since we are using a simple finite volume approach and backward Euler in time, we use linear interpolation in both space and time to compute $\bfdp$ from $\bfPsi$. Thus, the entries of the matrix $\bfP$ contain the interpolation weights. For linear interpolation in 3D, $\bfP$ is a sparse matrix that contains up to eight non-zero entries in each row. Note that the time and location of the data measurement is independent and decoupled from the numerical discretization used in the forward problem. A water retention curve, such as the van Genuchten model, can be used for computation of predicted water content data, which requires another nonlinear transformation, $\bfdp^\theta = \bfP \bftheta(\bfPsi(\bfm), \bfm)$. Note here that the transformation to water content data, in general, depends on the model to be estimated in the inversion, which will be addressed in the numerical examples. For brevity in the derivation that follows, we will make two simplifying assumptions: (1) that the data are pressure head measurements, which requires a linear interpolation that is not dependent on the model; and, (2) that the model vector, $\bfm$, describes only distributed saturated hydraulic conductivity. Our software implementation _does not_ make these assumptions; our numerical examples will use water content data, a variety of empirical relations, and calculate the sensitivity to multiple heterogeneous empirical parameters.

We can now formulate the discrete inverse problem to estimate saturated hydraulic conductivity, $\bfm$, from the observed pressure head data, $\bfdo$. We frame the inversion as an optimization problem, which minimizes a data misfit and a regularization term. Chapter \ref{ch:framework} showed an approach for geophysical inversions where hundreds of thousands to millions of distributed parameters are commonly estimated in a deterministic inversion {cite:p}`tikhonov1977, DougTutorial, Constable1987, haber2015computational`. Please refer to the previous chapter for the details of this inversion methodology. The hydrogeologic literature also shows the use of these techniques; however, there is also a large community advancing stochastic inversion techniques and geologic realism (cf. {cite:t}`Linde2015`). Regardless of the inversion algorithm used, an efficient method to calculate the sensitivity is crucial; this method is the focus of our work.

% (sec:richards-derivs)=

## Implicit sensitivity calculation

The optimization problem requires the derivative of the pressure head with respect to the model parameters, $\frac{\partial\bfPsi}{\partial\bfm}$. We can obtain an approximation of the sensitivity matrix through a finite difference method on the forward problem {cite:p}`Simunek1996, Finsterle2011, Finsterle2011c`. One forward problem, or two, when using central differences, must be completed for each column in the Jacobian at every iteration of the optimization algorithm. This style of differentiation proves advantageous in that it can be applied to any forward problem; however, it is highly inefficient and introduces errors into the inversion that may slow the convergence of the scheme {cite:p}`Doherty2004`. Automatic differentiation (AD) can also be used {cite:p}`nw`. However, AD does not take the structure of the problem into consideration and often requires that the dense Jacobian be explicitly formed. {cite:t}`Bitterlich2002` presents three algorithms (finite difference, adjoint, and direct) to directly compute the elements of the dense sensitivity matrix for the Richards equation. As problem size increases, the memory required to store this dense matrix often becomes a practical computational limitation {cite:p}`hao2, Towara2015`. As we show next, it is possible to explicitly write the derivatives of the Jacobian and evaluate their products with vectors using only sparse matrix operations. This algorithm is much more efficient than finite differencing, especially for large-scale simulations, since it does not require explicitly forming and storing a large dense matrix. Rather, the algorithm efficiently computes matrix-vector and adjoint matrix-vector products with sensitivity. We can use these products for the solution of the Gauss-Newton system when using the conjugate gradient method, which bypasses the need for the direct calculation of the sensitivity matrix and makes solving large-scale inverse problems possible. Other geophysical inverse problems have used this idea extensively, especially in large-scale electromagnetics (cf. {cite:t}`hao`). The challenge in both the derivation and implementation for the Richards equation lies in differentiating the nonlinear time-dependent forward simulation with respect to multiple distributed hydraulic parameters.

The approach to implicitly constructing the derivative of the Richards equation in time involves writing the whole time-stepping process as a block bi-diagonal matrix system. The discrete Richards equation can be written as a function of the model. For a single time-step, the equation is written:

% equations/richards/richards-timestep

```{math}
:label: eq:richards-timestep
\begin{align*}
F(\boldsymbol{\psi}^{n+1}(\mathbf{m}),\boldsymbol{\psi}^n(\mathbf{m}),\mathbf{m}) =
\frac{\boldsymbol{\theta}^{n+1}(\boldsymbol{\psi}^{n+1}) - \boldsymbol{\theta}^n(\boldsymbol{\psi}\n)}{\Delta t}
\nonumber\\-
\ \mathbf{D}
    \text{ diag}\left(\mathbf{k}_{Av}(\boldsymbol{\psi}^{n+1},\mathbf{m})\right)
    \mathbf{G} \boldsymbol{\psi}^{n+1}
-
\mathbf{G}_{z} \mathbf{k}_{Av}(\boldsymbol{\psi}^{n+1},\mathbf{m})
= 0
\end{align*}
```

In this case, $\bfm$ is a vector that contains all the parameters of interest. Note that $\bfpsi\nn$ and $\bfpsi\n$ are also functions of $\bfm$. In general, $\bftheta\nn$ and $\bftheta\n$ are also dependent on the model; however, for brevity, we will omit these derivatives. The derivatives of $\FF$ to the change in the parameters $\bfm$ can be written as:

% equations/richards/richards-timestep-chain

```{math}
:label: eq:richards-timestep-chain
\nabla_\mathbf{m}  F(\boldsymbol{\psi}^n,\boldsymbol{\psi}^{n+1},\mathbf{m})
=
\frac{\partial F}{\partial \mathbf{k}_{Av}}
\frac{\partial\mathbf{k}_{Av}}{\partial\mathbf{m}}
+ \frac{\partial F}{\partial \boldsymbol{\psi}^n}\frac{\partial\boldsymbol{\psi}^n}{\partial\mathbf{m}}
+ \frac{\partial F}{\partial \boldsymbol{\psi}^{n+1}}\frac{\partial\boldsymbol{\psi}^{n+1}}{\partial\mathbf{m}}
=0
```

or, in more detail:

% equations/richards/richards-timestep-deriv

```{math}
:label: eq:richards-timestep-deriv
\begin{align*}
\frac{1}{\Delta t}
\left(
    \frac{\partial \boldsymbol{\theta}^{n+1}}{\partial\boldsymbol{\psi}^{n+1}}
    \frac{\partial \boldsymbol{\psi}^{n+1}}{\partial\mathbf{m}}
    -
    \frac{\partial \boldsymbol{\theta}^n}{\partial\boldsymbol{\psi}^n}
    \frac{\partial \boldsymbol{\psi}^n}{\partial\mathbf{m}}
\right)
-
\mathbf{D}
    \text{ diag}\left( \mathbf{G} \boldsymbol{\psi}^{n+1} \right)
    \left(
        \frac{\partial \mathbf{k}_{Av}}{\partial\mathbf{m}} +
        \frac{\partial \mathbf{k}_{Av}}{\partial\boldsymbol{\psi}^{n+1}}
        \frac{\partial \boldsymbol{\psi}^{n+1}}{\partial\mathbf{m}}
    \right)
\nonumber\\
-\
\mathbf{D}
    \text{ diag}\left( \mathbf{k}_{Av}(\boldsymbol{\psi}^{n+1}) \right)
\mathbf{G}
\frac{\partial \boldsymbol{\psi}^{n+1}}{\partial\mathbf{m}}
-
\mathbf{G}_{z}
\left(
    \frac{\partial \mathbf{k}_{Av}}{\partial\mathbf{m}} +
    \frac{\partial \mathbf{k}_{Av}}{\partial\boldsymbol{\psi}^{n+1}}
    \frac{\partial \boldsymbol{\psi}^{n+1}}{\partial\mathbf{m}}
\right)
& = 0
\end{align*}
```

The above equation is a linear system of equations and, to solve for $\deriv{\bfPsi}{\bfm}$, we rearrange the block-matrix equation:

% equations/richards/richards-timestep-deriv-blocks

```{math}
:label: eq:richards-timestep-deriv-blocks
\begin{align*}
\overbrace{
    \left[
        \frac{1}{\Delta t}
        \frac{\partial \boldsymbol{\theta}^{n+1}}{\partial\boldsymbol{\psi}^{n+1}}
        -\mathbf{D}
        \text{ diag}\left( \mathbf{G} \boldsymbol{\psi}^{n+1} \right)
        \frac{\partial \mathbf{k}_{Av}}{\partial\boldsymbol{\psi}^{n+1}}
        -\mathbf{D}
        \text{ diag}\left( \mathbf{k}_{Av}(\boldsymbol{\psi}^{n+1},\mathbf{m}) \right)
        \mathbf{G}
        - \mathbf{G}_{z}
        \frac{\partial \mathbf{k}_{Av}}{\partial\boldsymbol{\psi}^{n+1}}
    \right]
}^{\mathbf{A}_0(\boldsymbol{\psi}^{n+1})}
\frac{\partial \boldsymbol{\psi}^{n+1}}{\partial\mathbf{m}}
\nonumber\\
+
\underbrace{
    \left[
        -\frac{1}{\Delta t}
        \frac{\partial \boldsymbol{\theta}^n}{\partial\boldsymbol{\psi}^n}
    \right]
}_{\mathbf{A}_{-1}(\boldsymbol{\psi}^n)}
\frac{\partial \boldsymbol{\psi}^n}{\partial\mathbf{m}}
=
\underbrace{
\left[
    -\mathbf{D}
    \text{ diag}\left( \mathbf{G} \boldsymbol{\psi}^{n+1} \right)
    \frac{\partial \mathbf{k}_{Av}}{\partial\mathbf{m}}
    -\mathbf{G}_{z}
    \frac{\partial \mathbf{k}_{Av}}{\partial\mathbf{m}}
\right]
}_{\mathbf{B}(\psi^{n+1})}&
\end{align*}
```

Here, we use the subscript notation of $\bfA_0(\bfpsi\nn)$ and $\bfA_{-1}(\bfpsi\n)$ to represent two block-diagonals of the large sparse matrix $\bfA({\bfPsi},\bfm)$. Note that all of the terms in these matrices are already evaluated when computing the Jacobian of the Richards equations in Section \ref{sec:richards-forward} and that they contain only basic sparse linear algebra manipulations without the inversion of any matrix. If $\bfpsi_0$ does not depend on the model, meaning the initial conditions are independent, then we can formulate the block system as:

% equations/richards/richards-timestep-deriv-matrix

```{math}
:label: eq:richards-timestep-deriv-matrix
\overbrace{
    \left[
    \begin{array}{cccc}
        \mathbf{A}_0(\boldsymbol{\psi}_1)\\
        \mathbf{A}_{-1}(\boldsymbol{\psi}_1)&\mathbf{A}_0(\boldsymbol{\psi}_2)\\
        &\mathbf{A}_{-1}(\boldsymbol{\psi}_2)&\mathbf{A}_0(\boldsymbol{\psi}_3)\\
        &\qquad\ddots&\qquad\ddots\\
        &&\mathbf{A}_{-1}(\boldsymbol{\psi}_{n_{t}-1})&\mathbf{A}_0(\boldsymbol{\psi}_{n_{t}})\\
    \end{array}
    \right]
}^{\mathbf{A}(\boldsymbol{\Psi},\mathbf{m})}
\overbrace{
    \left[
    \begin{array}{c}
        \frac{\partial \boldsymbol{\psi}_1}{\partial \mathbf{m}}\\
        \frac{\partial \boldsymbol{\psi}_2}{\partial \mathbf{m}}\\
        \vdots\\
        \frac{\partial \boldsymbol{\psi}_{n_{t}-1}}{\partial \mathbf{m}}\\
        \frac{\partial \boldsymbol{\psi}_{n_{t}}}{\partial \mathbf{m}}\\
    \end{array}
    \right]
}^{\frac{\partial\boldsymbol{\Psi}}{\partial\mathbf{m}}}
    =
\overbrace{
    \left[
    \begin{array}{c}
        \mathbf{B_1}(\boldsymbol{\psi}_1)\\
        \mathbf{B_2}(\boldsymbol{\psi}_2)\\
        \vdots\\
        \mathbf{B_{n-1}}(\boldsymbol{\psi}_{n_{t}-1})\\
        \mathbf{B_n}(\boldsymbol{\psi}_{n_{t}})\\
    \end{array}
    \right]
}^{\mathbf{B}(\boldsymbol{\Psi},\mathbf{m})}
```

This is a block matrix equation; both the storage and solve will be expensive if it is explicitly computed. Indeed, its direct computation is equivalent to the adjoint method {cite:p}`Bitterlich2002, DeanChen2011`.

Since only matrix vector products are needed for the inexact Gauss-Newton optimization method, the matrix $\bfJ$ is never needed explicitly and only the products of the form $\bfJ \bfv$ and $\bfJ^\top \bfz$ are needed for arbitrary vectors $\bfv$ and $\bfz$. Projecting the full sensitivity matrix onto the data-space using $\bfP$ results in the following equations for the Jacobian:

% equations/richards/richards-timestep-deriv-mult
%\begin{subequations}

```{math}
:label: eq:richards-timestep-deriv-mult
\mathbf{J} = \mathbf{P} \mathbf{A}(\boldsymbol{\Psi},\mathbf{m})^{-1} \mathbf{B}(\boldsymbol{\Psi},\mathbf{m})
```

```{math}
:label: eq:richards-timestep-deriv-mult-adjoint
\mathbf{J}^\top =   \mathbf{B}(\boldsymbol{\Psi},\mathbf{m})^\top \mathbf{A}(\boldsymbol{\Psi},\mathbf{m})^{-\top} \mathbf{P}^\top
```

%\end{subequations}

In these equations, we are careful to not write $\deriv{\bfPsi}{\bfm}$, as it is a large dense matrix which we do not want to explicitly compute or store. Additionally, the matrices $\bf A(\bfPsi,\bfm)$ and $\bf B(\bfPsi,\bfm)$ do not even need to be explicitly formed because the matrix $\bf A(\bfPsi,\bfm)$ is a triangular block-system, which we can solve using forward or backward substitution with only one block-row being solved at a time (this is equivalent to a single time step). To compute the matrix vector product, $\bfJ \bfv$, we use a simple algorithm:

1. Given the vector $\bfv$ calculate $\bfy = \mathbf{Bv}$
2. Solve the linear system $\mathbf{Aw} = \bfy$ for the vector $\bfw$
3. Set $\bfJ \bfv = \bfP \bfw$

Here, we note that we complete steps (1) and (2) using a for-loop with only one block-row being computed and stored at a time. As such, only the full solution, $\bfPsi$, is stored and all other block-entries may be computed as needed. There is a complication here if data is in terms of water content or effective saturation, as the data projection is no longer linear and may have model dependence. These complications can be dealt with using the chain rule during step (3). Similarly, to compute the adjoint $\bfJ^\top \bfz$ involves the intermediate solve for $\bfy$ in $\bfA^\top \bfy = \bfP^\top \bfz$ and then computation of $\bfB^\top \bfy$. Again, we solve the block-matrix via backward substitution with all block matrix entries being computed as needed. Note that the backward substitution algorithm can be viewed as time stepping, which means that it moves from the final time back to the initial time. This time stepping is equivalent to the adjoint method that is discussed in {cite:t}`DeanChen2011` and references within. The main difference between our approach and the classical adjoint-based method is that our approach yields the **exact** gradient of the discrete system; no such guarantee is given for adjoint-based methods.

The above algorithm and the computations of all of the different derivatives summarizes the technical details of the computations of the sensitivities. Equipped with this "machinery", we now demonstrate that validity of our implementation.

% (sec:richards-results)=

# Numerical results

The focus of this section is to validate and compare our algorithm and implementation to the literature. The following chapter will focus on applications of this work as well as demonstrate computational scalability of the algorithm for realistic field examples (Chapter \ref{ch:applications}).

% (sec:richards-validation)=

## Validation

### Forward problem

The Richards equation has no analytic solution, which makes testing the code more involved. Here we have chosen to use a fictitious source experiment to rigorously test the code. In this experiment, we approximate an infiltration front by an arctangent function in one dimension, which is centered over the highly nonlinear part of the van Genuchten curves, with $\psi\in[-60,-20]$ centimeters. The arctangent curve advects into the soil column with time. The advantage of using an analytic function is that the derivative is known explicitly and can be calculated at all times. However, it should be noted that the Richards equation does not satisfy the analytic solution exactly, but differs by a function, $S(x,t)$. We refer to this function as the fictitious source term. The analytic function that we used has similar boundary conditions and shape to an example in {cite:t}`Celia1990` and is considered over the domain $x\in[0, 1]$.

```{math}
:label: eq:fictitiousSource
\Psi_{\text{true}}(x,t)=-20\arctan(20((x-0.25)-t))-40
```

This analytic function is shown at times 0 and 0.5 in {numref}`Figure %s <fig:richards-validation-source>` and has a pressure head range of $\psi\in[-60,-20]$. We can compare these values to the van Genuchten curves in {numref}`Figure %s <fig:van-genuchten>`. We can then put the known pressure head into the Richards equation {eq}`eq:richards-mixed` and calculate the analytic derivatives and equate them to a source term, $S(x,t)$. Knowing this source term and the analytic boundary conditions, we can solve discretized form of the Richards equation, which should reproduce the analytic function in Equation {eq}`eq:fictitiousSource`. {numref}`Table %s <table:richards-source>` shows the results of the fictitious source test when the number of mesh-cells is doubled and the time-discretization is both fixed and equivalent to the mesh size (i.e. $k=h$). In this case, we expect that the backward-Euler time discretization, which is $\mathcal{O}(\delta)$, will limit the order of accuracy. The final column of {numref}`Table %s <table:richards-source>` indeed shows that the order of accuracy is $\mathcal{O}(\delta)$. The higher errors in the coarse discretization are due to high discontinuities and changes in the source term, which the coarse discretization does not resolve. We can complete a similar procedure in two and three dimensions and these tests show similar results of convergence. The rigorous testing of the code presented provides confidence in the forward simulation that is used throughout the following sections of this chapter.

```{figure} images/richards-validation-source.png
:name: fig:richards-validation-source
Fictitious source test in 1D showing the analytic function $\Psi_{\text{true}}$ at times 0.0 and 0.5 and the numerical solution $\Psi(x,0.5)$ using the mixed-form Newton iteration.
```

% tables/richards/source

```{list-table} Fictitious source test for Richards equation in 1D using the mixed-form Newton iteration.
:name: table:richards-source
:header-rows: 1
* - Mesh Size (n)
  - $||\Psi(x,0.5) - \Psi_{\text{true}}(x,0.5)||_\infty$
  - Order Decrease, $\mathcal{O}(\delta)$ \\[0.5em]
* - 64
  - 5.485569e+00
  - \-
* - 128
  - 2.952912e+00
  - 0.894
* - 256
  - 1.556827e+00
  - 0.924
* - 512
  - 8.035072e-01
  - 0.954
* - 1024
  - 4.086729e-01
  - 0.975
* - 2048
  - 2.060448e-01
  - 0.988
* - 4096
  - 1.034566e-01
  - 0.994
* - 8192
  - 5.184507e-02
  - 0.997
```

### Inverse problem

In order to test the implicit sensitivity calculation, we employ derivative and adjoint tests as described in {cite:t}`haber2015computational`. Given that the Taylor expansion of a function $f(\mathbf{m} + h \Delta \mathbf{m})$ is

```{math}
f(\mathbf{m}+ h  \Delta \mathbf{m}) = f(\mathbf{m}) + h  \mathbf{J} \Delta \mathbf{m} + \mathcal{O}(h^2),
```

for any of the model parameters considered, we see that our approximation of $f(\mathbf{m}+ h  \Delta \mathbf{m})$ by $ f(\mathbf{m}) + h \mathbf{J} \Delta \mathbf{m}$ should converge as $\mathcal{O}(h^2)$ as $h$ is reduced. This allows us to verify our calculation of $\mathbf{J}\mathbf{v}$. To verify the adjoint, $\mathbf{J}^\top\mathbf{v}$, we check that

```{math}
\mathbf{w}^\top\mathbf{J}\mathbf{v} = \mathbf{v}^\top \mathbf{J}^\top \mathbf{w}
```

for any two random vectors, $\mathbf{w}$ and $\mathbf{v}$. These tests are run for all of the parameters considered in an inversion of the Richards equation. Within our implementation, both the derivative and adjoint tests are included as unit tests which are run on any updates to the implementation (<https://travis-ci.org/simpeg/simpeg>)[^1].

[^1]: Testing has since moved to GitHub Actions and Azure Cloud.

## Comparison to literature

Code-to-code comparisons have been completed for comparison to {cite:t}`Celia1990`, which can be found in {cite:t}`richardsceliacomparison`. The following results are a direct comparison to the results produced by {cite:t}`Celia1990` for the Picard iteration only; {cite:t}`Celia1990` did not implement the Newton iteration. The direct comparison is completed: (a) to give confidence to the numerical results; (b) to compare the Newton iteration to the Picard iteration of the mixed formulation for a well-known example; and (c) demonstrate the use of multiple empirical models. Here, we use the Haverkamp model {cite:p}`Haverkamp1977` (rather than the classically used van Genuchten model) for the water retention and hydraulic conductivity functions.

% equations/richards/haverkamp

```{math}
\begin{align*}
    \theta(\psi) &= \frac{\alpha (\theta_s - \theta_r)}{\alpha + |\psi|^\beta} + \theta_r \\
    K(\psi) &= K_s \frac{A}{A+|\psi|^\gamma}
\end{align*}
```

We used parameters of $\alpha=1.611 \times 10^6$, $ \theta_s = 0.287$, $ \theta_r = 0.075$, $ \beta = 3.96$, $ A = 1.175 \times 10^6$, $ \gamma = 4.74$, and $K_s = 9.44 \times 10^{-5}$ m/s, which are the same as in {cite:t}`Celia1990`. The 40 cm high 1D soil column has initially dry conditions with a pressure head $\psi_0(x,0)=-61.5$cm. The boundary conditions applied are inhomogeneous Dirichlet with the top of the soil column, $\psi(40\text{cm},t)=-20.7$cm, and the bottom of the soil column, $\psi(0\text{cm},t)=-61.5$cm. The initial conditions are not consistent with the boundary condition at the top of the soil profile. This inconsistency leads to a boundary layer and a steep gradient in the pressure head at early times; as such, we anticipate that the Newton iteration will converge slowly at these times. The spatial grid is regular and has a spacing of 1.0cm, while the time-stepping, $\Delta t$, is manipulated. The three iterative methods described in Section \ref{sec:richards-forward} are implemented and compared at 360s: (1) head-based form Picard; (2) mixed-form Picard; and, (3) mixed-form Newton. {numref}`Figure %s <fig:richards-celia>` shows the solution obtained with the three iterative methods. Comparing the head-based formulation to the mixed-formulation for a large time-step (e.g. $\Delta t = 120$s) shows the degradation of the head-based method. Not only is the infiltration front smoothed, there is underestimate of the front location ({numref}`Figure %sa <fig:richards-celia>`). The underestimate of the infiltration front location is due to a loss of mass, which can be traced back the initial formulation of the head-based method {cite:p}`Celia1990`. The mixed-formulation, solved with either a Picard iteration or a Newton method, conserves mass and correctly identifies the spatial location of the infiltration front ({numref}`Figure %sb <fig:richards-celia>`). The results obtained here show excellent agreement with {cite:t}`Celia1990`.

```{figure} images/richards-validation-celia.png
:name: fig:richards-celia
Comparison of results to {cite:t}`Celia1990` showing the differences in the
(a) head-based and (b) mixed formulations for $t$=360s.
```

# Conclusions

The number of parameters that are estimated in the Richards equation inversions has grown and will continue to grow as time-lapse data and geophysical data integration become standard in site characterizations. The increase in data quantity and quality provides the opportunity to estimate spatially distributed hydraulic parameters from the Richards equation; doing so requires efficient simulation and inversion strategies. In this chapter, we have shown a derivative-based optimization algorithm that does not store the Jacobian, but rather computes its effect on a vector (i.e. $\bf Jv$ or $\bf J^\top z$). By not storing the Jacobian, the size of the problem that we can invert becomes much larger. We have presented efficient methods to compute the Jacobian that can be used for all empirical hydraulic parameters, even if the functional relationship between parameters is obtained from laboratory experiments.

Our technique allows a deterministic inversion, which includes regularization, to be formulated and solved for any of the empirical parameters within the Richards equation. For a full 3D simulation, as many as ten spatially distributed parameters may be needed, resulting in a highly non-unique problem; as such, we may not be able to reasonably estimate all hydraulic parameters. Depending on the setting, amount of a-priori knowledge, quality and quantity of data, the selection of which parameters to invert for may vary. Our methodology enables practitioners to experiment in 1D, 2D and 3D with full simulations and inversions, in order to explore the parameters that are important to a particular dataset. Our numerical implementation is provided in an open source repository (<https://simpeg.xyz>) and is integrated into the framework presented in Chapter \ref{ch:framework}. The goal of Chapter \ref{ch:applications} is to work with these techniques in an experiment that represents a field infiltration experiment. The following chapter will also further document the scalability of this approach when moving to 3D inversions.
