# OpenPipeline Spatial

Extensible spatial single cell analysis pipelines for reproducible and large-scale spatial single cell processing using Viash and Nextflow.

[![ViashHub](https://img.shields.io/badge/ViashHub-openpipeline-7a4baa.svg)](https://www.viash-hub.com/packages/openpipeline_spatial)
[![GitHub](https://img.shields.io/badge/GitHub-viash--hub%2Fopenpipeline-blue.svg)](https://github.com/openpipelines-bio/openpipeline_spatial)
[![GitHub
License](https://img.shields.io/github/license/openpipelines-bio/openpipeline.svg)](https://github.com/openpipelines-bio/openpipeline_spatial/blob/main/LICENSE)
[![GitHub
Issues](https://img.shields.io/github/issues/openpipelines-bio/openpipeline.svg)](https://github.com/openpipelines-bio/openpipeline_spatial/issues)
[![Viash
version](https://img.shields.io/badge/Viash-v0.9.3-blue.svg)](https://viash.io)

## Overview

OpenPipeline Spatial extends the OpenPipeline ecosystem with specialized workflows and components for spatial transcriptomics analysis. It provides standardized, reproducible pipelines that are technology-agnostic and can be used for processing spatial omics data from various technologies and platforms, including 10x Genomics Xenium, NanoString CosMx and AtoMx, and Element Biosciences Aviti.

## Functionality

OpenPipeline Spatial executes a list of predefined tasks specifically designed for spatial omics data. These discrete steps are also provided as standalone components that can be executed individually with a standardized interface.

The following spatial-specific workflows are provided:

- [Ingestion](https://www.viash-hub.com/packages/openpipeline_spatial/latest/components?search=mapping): Whereas many technologies generate count matrices on-instrument, functionality is provided for the mapping & quantification of 10X Visum data.
- [Interoperability](https://www.viash-hub.com/packages/openpipeline_spatial/latest/components?search=convert): To make sure all spatial workflows are technology-agnostic, functionality is provided to convert count matrices from different technologies into a common format (H5Mu). In addition, functionality is provided to convert between various Spatial data formats (e.g. Seurat, SpatialExperiment, MuData, SpatialData).
- [QC](https://www.viash-hub.com/packages/openpipeline_spatial/latest/components?search=spatial_qc): Calculation of comprehensive quality control metrics.
- [Sample Processing](https://www.viash-hub.com/packages/openpipeline_spatial/latest/components?search=spatial_process_samples): Batch processing of multiple spatial samples, including count-based filtering, normalisation and dimensionality reduction.

## Extended functionality

Whereas this package only provides spatial-specific functionality, it is designed to work seamlessly with the core OpenPipeline package. This means that all core OpenPipeline workflows and components can be used in conjunction with the spatial-specific ones. For example, the **integration** and **cell type annotation** workflows can be applied to spatial data after it has been processed using the spatial-specific workflows.

``` mermaid lang="mermaid"
flowchart LR
  demultiplexing["Step 1: Ingestion"]
  ingestion["Step 2: QC"]
  process_samples["Step 3: Process Samples"]
  integration["Step 4: Integration"]
  downstream["Step 5: Downstream Analysis"]
  demultiplexing-->ingestion-->process_samples-->integration-->downstream
```

## Execution via CLI or Seqera Cloud

The openpipeline_spatial package is available via [Viash
Hub](https://www.viash-hub.com/packages/openpipeline_spatial/latest/), where
you can receive instructions on how to run the end-to-end workflow as
well as individual subworkflows or components.

It’s possible to run the workflow directly from Seqera Cloud. The necessary Nextflow schema files have been [built and provided with the workflows](https://packages.viash-hub.com/vsh/openpipeline_spatial/-/tree/build/main/target/nextflow?ref_type=heads) in order to use the form-based input. However, Seqera Cloud can not deal with multiple-value parameters for batch processing of multiple samples. Therefore, it’s better to use Viash Hub also here for launching the workflow on Seqera Cloud.

* Navigate to the [Viash Hub package page](https://www.viash-hub.com/packages/openpipeline_spatial/latest/), select the workflow you want to launch and click the `launch` button.
* Select the execution environment of choice (e.g. `Seqera Cloud`, `CLI` or `Executable`)
* Fill in the form with the required parameters and launch the workflow.

## Support
For issues specific to spatial analysis, please use the GitHub issues tracker. For general OpenPipeline questions, refer to the main OpenPipeline documentation.