# COMAVI documentation

Use this page to find the shortest path for your task. If you are new to COMAVI,
start with the browser notebook or the worked local example.

## Run COMAVI

| Goal | Start here |
|---|---|
| Score variants in a browser | [Open the guided Colab notebook](https://colab.research.google.com/github/la424/comavi/blob/main/notebooks/COMAVI_colab.ipynb) |
| Prepare or reuse a structural system | [System-setup wizard](COMAVI_SYSTEM_SETUP_WIZARD.md) |
| Check a local installation without FoldX | [Worked HBB–HBA1 example](../examples/README.md) |
| Configure genes from the command line | [Adding your own genes](adding_your_own_genes.md) |
| Verify the published benchmark | [Verification entry point](../verification/README.md) |

The Colab workflow is recommended for a first run. The command-line workflow is
best for automation, larger batches, and already-curated structure sets.

## Interpret the outputs

- [Results synthesis](COMAVI_results_synthesis.md) — plain-language walkthrough
  of the benchmark and result fields.
- [Canonical benchmark ledger](COMAVI_v7_canonical_benchmark_ledger.md) —
  authoritative metrics, provenance, and recomputation record.
- [Design decisions](design_decisions.md) — why the operating pipeline makes its
  current modeling and reporting choices.
- [`reference_outputs/`](../reference_outputs/) — committed CSV, workbook, JSON,
  and report artifacts.
- [`figures/`](../figures/README.md) — index of manuscript figures and their data
  sources.

COMAVI reports modeled structural disruption and a proposed mechanism. It does
not issue a pathogenicity verdict.

## Setup-wizard validation

- [Supported behavior and limits](COMAVI_SYSTEM_SETUP_WIZARD.md)
- [Validation report](COMAVI_SYSTEM_SETUP_WIZARD_VALIDATION.md)
- [Machine-readable validation record](COMAVI_SYSTEM_SETUP_WIZARD_VALIDATION.json)
- [Public Colab workflow validation](COLAB_PUBLIC_WORKFLOW_VALIDATION.md)

## Study reproduction and historical records

- [Methods and metrics sketch](methods_metrics_sketch.md)
- [Pre-publication checkpoint](CHECKPOINT_pre_publication.md)
- [Phase 4 handoff](PHASE4_HANDOFF.md)
- [`archive/`](../archive/) — superseded derivation and grading material; not part
  of the operating pipeline.
- [`maintainers/`](maintainers/) — release and repository-administration records;
  not required to run COMAVI.

## Contributing or reporting a problem

See [`CONTRIBUTING.md`](../CONTRIBUTING.md) for supported environments, checks,
and the information to include in a reproducible issue.
