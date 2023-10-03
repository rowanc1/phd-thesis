---
title: Conclusions
description: ''
subject: Chapter 5
numbering:
  enumerator: 5.%s
---

+++

Forward and inverse modelling in geophysics requires solving and optimizing large-scale partial differential equations. Thus, many components including linear algebra, optimization routines, discretizations, and model conceptualizations, are required to interact. Advances in instrumentation, the monitoring of time-lapse processes, and acquisition of multiple geophysical and hydrologic data types are driving a need for a more integrated geoscience approach. This includes integrating information through multiphysics simulations and coupled geophysical inversions. In addition to recovering 3D models, geophysics is increasingly being used in time-lapse imaging of fluid flow processes, requiring both computational scalability to 4D inverse problems and interdisciplinary collaborations. There have been significant advances in the past three decades in computational geophysics; researchers are now able to both simulate and invert almost any geophysical method in three dimensions (cf. {cite:t}`Oldenburg2016`). One of the major challenges ahead of geophysics as a discipline is how to systematically improve the quantitative interfaces and integrations between hydrological, geophysical, and geological information and processes. The research necessary to address these challenges will require the interdisciplinary community to build upon, as well as augment, standard practices; this presupposes that researchers have access to consistent methodologies that can be extended, adapted, and combined.

Adapting interdisciplinary methodologies to geophysical simulations and inversions inherently requires that a diverse suite of methods and applications be considered across hydrogeology, geophysics, and geology. Throughout this work, my colleagues and I have tried to summarize and/or reproduce many methodologies in the geophysical inversion literature. Other communities, such as astrophysics and machine learning {cite:p}`astropy, scikit-learn` have organized these efforts and research communities around open, accessible, and actionable ideas. Adapting these learnings to geoscience, I strive to complete all of my research such that it is immediately reproducible and openly available.

Much of my dissertation required the development of software, which is implementation and engineering by its nature. However, the software, although extremely useful, is not the aim of my research. If the goal is a framework for quantitative geoscience integration, a simplistic, high level conceptualization is easy to present. However, a picture or a paragraph describing a framework cannot be tested, interrogated, nor used beyond its static form. Software is the means by which I test, extend, organize, and abstract inherently computational ideas in a rigorous, scientific way. The aim of my research is to identify, explore, and formalize a framework for simulation and parameter estimation in geophysics. I applied this framework to vadose zone flow and, with the help of collaborators, to several other geophysical methodologies.

# Contributions and dissemination

The conclusions from each component of my thesis are contained within each chapter and each appendix. However, the central aim of my dissertation was to develop a framework for geophysical inversions: this has been disseminated through three publications {cite:p}`Cockett2017, Heagy2016, simpeg2015`, four extended abstracts {cite:p}`HeagySEG2014, Kang2015, Heagy2015, HeagyEM62017`, over twenty conference presentations, and a dedicated international workshop[^1]. Furthermore, this framework has demonstrated value in several new research areas, methodologies, and case studies (cf. {cite:t}`Kang2017, Miller2017, Kang2016, Rosenkjaer2016`). Overall, the contributions of my thesis are twofold:

[^1]: At the Banff International Research Station, see <https://www.birs.ca/events/2016/2-day-workshops/16w2695>

1. a conceptual organization and synthesis of geophysical simulations and inversions into a framework that has been rigorously, numerically tested; and
2. an algorithm for large-scale vadose zone parameter estimation for any distributed hydraulic parameter, regardless of the empirical relationship used.

One outcome of a framework approach is the accelerated transfer of ideas from one discipline to another. For example, the implicit sensitivity calculation for the Richards equation {cite:p}`Cockett2017` was heavily inspired by research completed in time-domain electromagnetics simulations and inversions {cite:p}`Heagy2016`. The refinement and application of this algorithm to hydrology significantly improved numerical scalability for the 3D inverse problem. Chapter \ref{ch:applications} showed significant improvements in memory over explicitly forming the sensitivity matrix by over two orders of magnitude for the example shown, bringing this inversion into the range of possibility on modest computational resources. Additionally, the complexities of the Richards equation were generalized and synthesized to improve other geophysical methods in the framework. These improvements were especially demonstrated with regard to dealing with multiple physical properties that may or may not be estimated and occur in distributed empirical relations throughout the forward simulation framework.

The framework presented is designed to decouple concerns and expose well-defined interfaces between the many components necessary in geophysical simulations and inversions. Chapter \ref{ch:applications} and Appendix \ref{ch:extensions} showed many demonstrations of these ideas, for example: (a) exposing model conceptualization that are decoupled from the sensitivity calculation allows custom, parameterized inversions to be completed by combining various, predefined mappings; (b) the declarative interface of differential operators and derivatives is decoupled from the structure and type of mesh, which allows the physics to be written once and used across many types of meshes; or (c) the decoupling of the physics from the definition of the geophysical survey, which allows many types of geophysical surveys to be combined with a single physical problem (e.g. in electromagnetics) or conversely different approximations of a physical problem (e.g. dimensionality) to be combined with a single survey. Given the number of choices in geophysical simulations and inversions this type of combinatorial, decoupled approach could provide significant acceleration to the unique geoscience integration problems of the future.

## Software

The conceptual organization developed could not have been created without the aid of software. Furthermore, even if it were, it would be of limited utility, difficult, or impossible to validate, and would not make significant progress towards sustained, reproducible, quantitative geoscience integrations. It is not until a framework is implemented and tested from a number of non-overlapping geoscience perspectives that the assumptions, inconsistencies, or redundancies come to light and are available to interrogation. The main software package that was, and continues to be, developed is {sc}`SimPEG` (<https://github.com/simpeg/simpeg>), which defines the framework and hosts a collection of other geophysical methods written by many collaborators across six universities. Currently there are methods for: vadose zone flow {cite:p}`Cockett2017`, direct current resistivity and induced polarization {cite:p}`Kang2016`, time-domain and frequency-domain electromagnetics {cite:p}`Heagy2016`, magnetotellurics {cite:p}`Rosenkjaer2016`, magnetics and gravity {cite:p}`Fournier2016, Miller2017`, and a number of example linear inverse problems. Many of these geophysical methods also have different formulations (e.g. integral equation, differential equation, etc.), dimensionalities (e.g. 1D, 2.5D, 3D), and survey components (e.g. sources and receivers). All software is disseminated with the MIT license to encourage permissionless innovation.

# Outlook and continuing work

This thesis constructed a preliminary organization and synthesis of simulations and deterministic inversions in a few subdisciplines of geophysics and hydrogeology. This conceptual framework and computational implementation has demonstrated utility in advancing current and future research, however, caveats and qualifiers abound.

Some of the more intricate parallel algorithms (e.g. stochastic optimization) would require updates to the implementation and framework; a lack of these scalable parallel algorithms is limiting when tackling problems with many sources (e.g. airborne electromagnetics). There is still significant work to do in coupling and nesting various geophysical problems in a robust way. Due to the current focus on enabling researchers, the framework offers little in high-level interfaces to inversions that are typical of _black box_ industry use (i.e. data in, model out). Regardless of these shortcomings, many of my colleagues rely upon and extend this framework and implementation in their day to day research. My goal was to accelerate their work and connect their research to a community of collaborators who are explicitly working towards common goals.

This thesis is positioned from a perspective of looking out to a future of multidisciplinary, multi-data-type, quantitatively integrated geoscience simulations and inversions. A future where joint, cooperative, coupled, parameterized, multiphysics inversions and simulations are the norm rather than the exception. Where multiple existing, robust, and computationally-efficient methodologies are combined to extract all possible information from disparate geoscience datasets. To realize this sort of ubiquitous, quantitative communication between disciplines and methodologies requires an organized and integrated community that can effectively work together. My research is _aimed_ here. This future will not be realized by one person nor by one research group. In the field of machine learning, {cite:t}`olah2017research` note that "[t]he maintainable size of the field is controlled by how its members trade off the energy between communicating and understanding." The curation of ideas is just as important as their creation. There is a _research debt_ created by an exclusive focus on research novelty; future advances also require distillation, synthesis, and explicit communication. There is a significant amount of effort ahead of us to achieve effective communication and collaboration with our geoscience peers in geology and hydrogeology. This communication and quantitative integration is the webbing on which the future of our _geoscience_ field depends. My approach, therefore, has been to research and disseminate a numerical framework that attempts to support and enable a number of these diverse interdisciplinary collaborations. I hope that a lasting contribution of my work is the open, modular approach that I have curated and the community that I have helped to seed around these ideas.