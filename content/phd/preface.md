---
title: Preface
description: The research for this dissertation was completed while studying at the University of British Columbia. This research has resulted in three peer reviewed publications, three expanded conference abstracts, and several auxiliary works.
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

The research for this dissertation was completed while studying at the University of British Columbia. This research has resulted in three peer reviewed publications, three expanded conference abstracts, and several auxiliary works. The main focus of my thesis is on a framework for geophysical simulations and inversions that increases quantitative geoscience communication. In 2016, Dr. Oldenburg, Dr. Pidlisecky, Lindsey Heagy and I organized an international conference around this work that was sponsored by the Banff International Research Station; excerpts from the introduction of my thesis were used in the conference proposal.

Chapter \ref{ch:framework} presents a framework for simulation and parameter estimation for geophysical applications. An earlier version of which was published in {cite:t}`simpeg2015`, and ideas from this chapter also have been presented at several international conferences (cf. {cite:t}`Cockett2014c, Cockett2015c, Cockett2015a`).

Chapter \ref{ch:richards} presents a computationally scalable algorithm for solving inverse problems for hydraulic parameters in vadose zone flow using the Richards equation. This work has been submitted for peer review and the preprint is available on _arXiv_ {cite:p}`Cockett2017`; preliminary versions of this research were presented at two conferences {cite:p}`Cockett2013, Cockett2013a`.

Chapter \ref{ch:applications} involves several numerical examples, which were inspired by work from my undergraduate thesis, of which two papers were published during the course of my graduate research {cite:p}`Pidlisecky2013, Cockett2014`. The forward simulation framework for multi-parameter simulations and inversions in time-domain physical problems used in this chapter was derived from collaborative work between electromagnetics and vadose zone flow {cite:p}`Heagy2016`. One of the numerical examples in Chapter \ref{ch:applications} has previously been published in {cite:p}`Cockett2017`.

Two of the appendices contain supporting materials on finite volume and several numerical examples and case studies. Appendix \ref{ch:discretize} on finite volume contains work and figures that have been published in a computational tutorial on finite volume {cite:p}`fvtutorial`. Additionally, much of this work is supported by course material and instruction from Dr. Eldad Haber, Dr. Uri Ascher, and Dr. Chen Grief {cite:p}`haber2015computational, Ascher2011`. Appendix \ref{ch:extensions} presents an adaptation of the forward simulation framework published in {cite:t}`Heagy2016` for the Richards equation. This appendix also summarizes conclusions and insights from three extended conference abstracts on electromagnetics and a publication on parametric geologic modelling {cite:p}`HeagySEG2014, Kang2015, Heagy2015, Cockett2016`.

Throughout the course of my graduate research, I have started and contributed to several open source software projects to support, test, and validate the geophysical simulation and inversion framework that is the main focus of my thesis. My main focus with this software was on inheritance, composition, terminology, and the interfaces between simulation and inversion components -- the elements that define the framework. This is demonstrated by my personal contribution of 267,614 lines of code over the last five years, which have been reduced over 4.5 fold to 59,111 lines of code while increasing possibilities and geophysical applications. For an up-to-date, detailed analysis on code contribution and attribution over time, please see: <https://www.openhub.net/p/simpeg-geophysics>. This is perhaps the most salient distinguishment between the focus on framework development as opposed to a script or executable that is aimed at a specific type of geophysical inversion. The major software packages that have been created are:

1. `SimPEG`, a framework for simulation and parameter estimation in geophysics (<https://github.com/simpeg/simpeg>);
2. `discretize`, a finite volume package for simulation in the context of inverse problems (<https://github.com/simpeg/discretize>); and
3. `pymatsolver`, a common interface to several matrix solvers and packages (<https://github.com/rowanc1/pymatsolver>).

These projects have seen significant investment from my colleagues in testing, applying, and expanding the capabilities of the framework to other geophysical applications. This open, collaborative work has involved colleagues across industry, government, and six universities. Currently {sc}`SimPEG` includes methods for: vadose zone flow {cite:p}`Cockett2017`; direct current resistivity and induced polarization {cite:p}`Kang2016`; time-domain and frequency-domain electromagnetics {cite:p}`Heagy2016, HeagyEM62017`; magnetotellurics {cite:p}`Rosenkjaer2016`; magnetics and gravity {cite:p}`Miller2017`; and several examples of other inverse problems {cite:p}`mttutorial`. All software has been released under the permissive MIT license, to encourage reuse, adaptation, and sustained contribution to these ideas.
