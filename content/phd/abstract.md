---
title: Abstract
description: ''
date: 2021-08-25T20:06:05.379Z
authors:
  - name: Rowan Cockett
    userId: vKndfPAZO7WeFxLH1GQcpnXPzfH3
    orcid: 0000-0002-7859-8394
    corresponding: null
    email: null
    roles: null
    affiliations: null
name: blank
oxa:
---

+++

```{admonition} Lay Summary
Geophysical methods gather data remotely to enable insights into subsurface structure and processes (e.g. locating economic resources or monitoring environmental changes). The information derived from geophysical methods is of crucial importance in resource exploration, environmental remediation, and the study of deep-earth processes. Interpretation of geophysical data requires a combination of numerical simulation and inversion. Inversion is a procedure for using data to estimate an image or model of the earth (this is similar to medical imaging). Increasingly, geoscientists are tackling complex problems that require integration of multiple types of information in order to better characterize the subsurface. In hydrogeology and geophysics, this quantitative integration requires advances in both disciplines, as well as a framework for this collaboration. The objective of this dissertation is to identify and refine a computational framework that enables and encourages sustained cross-disciplinary communication, which is a necessary step in integrated geophysical simulation research.
```

+++

## A framework for geophysical inversions with application to vadose zone parameter estimation

Inverse modeling is a powerful tool for extracting information about the subsurface from geophysical and hydrologic data. Geophysical inverse problems are inherently multidisciplinary, requiring elements from the relevant physics, numerical simulation, and optimization, as well as knowledge of the geologic setting, hydrologic processes, and a comprehension of the interplay between all of these elements. Increasingly geoscientists are tackling complex problems that require integration of multiple types of information in order to better characterize the subsurface. However, many of the sub-fields of geophysics are developing simulation and inversion approaches, algorithms, and supporting software in isolation. This isolation is a barrier to quantitative integration and leads to inefficiencies in advancing interdisciplinary research. Greater efficiencies, and higher quality outcomes, could be achieved if (hydro)geophysicists had a common framework to accelerate an integrated approach.

The development of a geophysical framework requires considering a number of disciplines and geophysical problems to ensure generality as well as extensibility. As such, I have worked with diverse collaborators on geophysical problems (e.g. electromagnetics and potential fields) to inform and research the structure of this framework. However, the goal is also to have the framework work outside of geophysics and most notably in hydrogeology; I have focused on vadose zone fluid flow as a model problem. Fluid flow in the vadose zone is governed by the Richards equation; it is parameterized by hydraulic conductivity, which is a nonlinear function of pressure head. Currently the computational scalability of the Richards equation inversion is a barrier for three dimensional inversions in hydrogeophysics. Existing methods explicitly calculate the sensitivity matrix using finite difference or automatic differentiation, however, for large-scale problems these methods are constrained by computation and memory, respectively. This dissertation provides an implicit sensitivity algorithm that enables large-scale inversion problems for distributed parameters in the Richards equation to become tractable on modest computational resources. The main goal of my thesis is to organize the components of (hydro)geophysical simulations and inverse problems, and synthesize these into a comprehensive, modular framework.
