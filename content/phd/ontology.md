---
title: Frameworks & Ontologies
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
subject: Appendix A
journal: PhD Thesis
open_access: true
license: CC-BY-4.0
doi: 10.14288/1.0362383
numbering:
  enumerator: A.%s
  heading_1: true
  heading_2: true
  heading_3: false
  heading_4: false
  heading_5: false
  heading_6: false
---

+++

A computational science framework provides a set of standards such that individual scientists can contribute software components, in this case components used in simulation or inversion routines, with the confidence that those components will work with other components in that framework. As such, the standards of the framework define the responsibilities of each component (or class of component) and the required interfaces of all components. The term framework is commonly used in the software community, however, the formal term for the component organization and their properties is an ontology. The use of ontologies in the sciences to formally describe domain knowledge has exploded in recent decades, especially in domains of artificial intelligence, chemistry, and biology, but also more recently in the geosciences {cite:p}`Sharman2004, Ma2011`. A computational ontology (rather than the underlying discipline of philosophical ontology) is a "formal explicit specification of a shared conceptualization" {cite:p}`Gruber1993`. {cite:t}`Noy2002` summarizes the purpose of computational ontologies as enabling a shared understanding of the structure of information and systematically enabling knowledge and information reuse. Practically ontologies can (a) provide access and discoverability to heterogeneous information; (b) act as a common language to lower the barrier to transfer of ideas; and (c) act as a specification for interoperability, for example, as a communication protocol or application programming interface {cite:p}`Sharman2004`. The techniques for building ontologies amount to capturing, synthesizing, organizing, and digitizing the relationships between concepts, conceptual inheritance patterns, and behaviour. Ontologies are most commonly used for storing and organizing data, for example, connecting genetic data with phenotypic data in bioinformatics. However, ontologies are also used in defining tasks, workflows, and problem solving methods {cite:p}`Fensel1997, Bard2004`. In more mature interdisciplinary fields, this research is becoming core to scientists' day to day research; for example, {cite:t}`Stein2008` notes that "all current biomedical cyberinfrastructure efforts use ontologies." As a result of successes in other fields, geoscience integration is currently the target of major funding initiatives across the world (e.g. EarthCube - 11 year NSF project, \$35M in 2015; CIMIC Footprints Project - NSERC Project, 24 Universities, 30 Industry, \$13M). Many of the current efforts are focused on computational science frameworks, formally describing geoscientific data (using ontologies), and formally describing methods of integrating disciplines. For example, a Common Component Architecture for high performance scientific computing {cite:p}`Armstrong1999` has been used as the basis for coupled forward integration of a number of geoscience simulation tools written by different authors {cite:p}`Peckham2013`. The research into these domain specific standards for interoperability is critical for sustainable interdisciplinary research.

The growth in complexity of geophysical data and analysis and the necessity for cross-disciplinary integrations is also coincident with the revolution of open source software communities, largely enabled through web-based interactions. Other research communities, for example `Astropy` in astronomy and `SciPy` in numerical computing, have embraced the open source approach for collaboration and research {cite:p}`astropy, scipy`. These pioneering efforts are now complemented by easy-to-use, ubiquitous web-based repositories and version-control systems (e.g. GitHub), that have removed many of the barriers associated with management and collaboration. The growth of such systems, coupled with the maturity of individual geophysical subdisciplines (e.g. potential fields, electromagnetics), presents an opportunity to develop a computational framework and associated ontology for geophysical simulation and inversion. An ontology is an embodiment of concepts, relationships, and behaviours in a specific scientific domain and can be (a) captured in special purpose languages (e.g. Web Ontology Language, Resource Description Framework), or (b) captured in general purpose computer programming languages (e.g. Python, Java, C++) {cite:p}`Sharman2004`. To research a geophysical simulation and inversion framework, I have chosen the latter approach for the purposes of utility, testing, and creating a framework/ontology that can be openly used and evolved by the geoscience community.
