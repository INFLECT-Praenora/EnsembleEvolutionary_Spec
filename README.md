# EnsembleEvolutionary_Spec

Specification documents for EnsGen, an evolutionary ensemble generator library
for patient ODE models, and its supporting process library builder.

## Documents

All documents are in the [`docs`](docs) folder.

### EnsGen library

- [EnsGen: Evolutionary Ensemble Generator Library - Specification](docs/EnsGen%20Evolutionary%20Ensemble%20Generator%20Library%20%E2%80%94%20Specification.md)
  - C/C++ shared library that evolves an ensemble of ODE models fitted to a
    patient's longitudinal time-series data, and keeps the ensemble current as
    new data arrive.
- [ensgen-architecture.png](docs/ensgen-architecture.png) - architecture diagram.

### Process Library Builder (plbuild)

- [EnsGen Process Library Builder - Specification](docs/plbuild-spec.md)
  - Offline tool that turns a chosen set of BioModels models into the curated
    process-library file EnsGen loads.
- [plbuild-pipeline.png](docs/plbuild-pipeline.png) - pipeline diagram.

### Presentation

- [Ensemble Generation for the CIRCLE Immune Digital Twin (PDF)](docs/Ensemble%20Generation%20for%20the%20CIRCLE%20Immune%20Digital%20Twin.pdf)
  - viewable in the browser.
- [Ensemble Generation for the CIRCLE Immune Digital Twin (PowerPoint)](docs/Ensemble%20Generation%20for%20the%20CIRCLE%20Immune%20Digital%20Twin.pptx)
  - editable source; download to view.
