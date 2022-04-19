---
title: Dissemination
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
subject: Appendix D
journal: PhD Thesis
open_access: true
license: CC-BY-4.0
doi: 10.14288/1.0362383
numbering:
  enumerator: D.%s
  heading_1: true
  heading_2: true
  heading_3: false
  heading_4: false
  heading_5: false
  heading_6: false
---

+++

Throughout the development of {sc}`SimPEG`, we have focused on building both a framework ({numref}`Figure %s <fig:classOutline>`) and a toolbox that is flexible and extensible. The toolbox includes utilities that we use repeatedly in our research, for example, visualization routines, mesh generation, and synthetic modeling. Additionally, when writing new code for differential operators and PDE systems, several functions are available to help test and verify results including: (a) checking derivatives and expected order of convergence, (b) comparing to analytics, and (c) adjoint tests for the sensitivity operators (cf. {cite:t}`haber2015computational`). These tests can be written in-line in an interactive development paradigm and then rapidly transfered and incorporated as unit-tests to ensure future code changes do not change functionality. Currently the entire {sc}`SimPEG` project has upwards of 80\% test coverage, with all core functionality (e.g. discretization, optimization) being extensively tested. All changes are combined using the practice of continuous integration, supported by freely available tools for open source projects such as TravisCI and GitHub. Additionally, we have focused on writing and maintaining documentation for core functionality (<https://docs.simpeg.xyz>). The documentation is also practicing continuous integration and is updated with all changes to {sc}`SimPEG` and is built and hosted through ReadTheDocs {cite:p}`RTFD`. We are hopeful that our efforts and adherence to best practices in open source code development (cf. {cite:t}`Wilson2014`) will encourage a community to exercise reproducible research around these scientific tools (cf. {cite:t}`Fomel2009`).
