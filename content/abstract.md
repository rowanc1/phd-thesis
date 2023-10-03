---
title: Abstract
description: The development of a geophysical framework requires considering a number of disciplines and geophysical problems to ensure generality as well as extensibility.
numbering: false
---

+++

Inverse modeling is a powerful tool for extracting information about the subsurface from geophysical and hydrologic data. Geophysical inverse problems are inherently multidisciplinary, requiring elements from the relevant physics, numerical simulation, and optimization, as well as knowledge of the geologic setting, hydrologic processes, and a comprehension of the interplay between all of these elements. Increasingly geoscientists are tackling complex problems that require integration of multiple types of information in order to better characterize the subsurface. However, many of the sub-fields of geophysics are developing simulation and inversion approaches, algorithms, and supporting software in isolation. This isolation is a barrier to quantitative integration and leads to inefficiencies in advancing interdisciplinary research. Greater efficiencies, and higher quality outcomes, could be achieved if (hydro)geophysicists had a common framework to accelerate an integrated approach.

The development of a geophysical framework requires considering a number of disciplines and geophysical problems to ensure generality as well as extensibility. As such, I have worked with diverse collaborators on geophysical problems (e.g. electromagnetics and potential fields) to inform and research the structure of this framework. However, the goal is also to have the framework work outside of geophysics and most notably in hydrogeology; I have focused on vadose zone fluid flow as a model problem. Fluid flow in the vadose zone is governed by the Richards equation; it is parameterized by hydraulic conductivity, which is a nonlinear function of pressure head. Currently the computational scalability of the Richards equation inversion is a barrier for three dimensional inversions in hydrogeophysics. Existing methods explicitly calculate the sensitivity matrix using finite difference or automatic differentiation, however, for large-scale problems these methods are constrained by computation and memory, respectively. This dissertation provides an implicit sensitivity algorithm that enables large-scale inversion problems for distributed parameters in the Richards equation to become tractable on modest computational resources. The main goal of my thesis is to organize the components of (hydro)geophysical simulations and inverse problems, and synthesize these into a comprehensive, modular framework.