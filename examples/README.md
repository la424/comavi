# Worked examples

Each subdirectory is a complete, self-contained MAVIS input set with the frozen
benchmark values it should reproduce. Use these to verify an installation, to see
what a valid config and variant file look like before writing your own, or as a
template to copy.

Everything except the final scoring step runs **without a FoldX licence**:
structures come from RCSB, and `run.py --dry-run` validates the config and expands
the variant list without invoking FoldX.

| Example | System | Structure | What it demonstrates |
|---|---|---|---|
| `kras_craf/` | KRAS - RAF1 RBD-CRD | 6XI7, x-ray 1.95 A | Experimental complex; three oncogenic variants that are **structurally silent** on all three axes |

## Why the KRAS example is the one to start with

The three variants (G12D, G12V, Q61H) are unambiguously pathogenic and among the
best-characterised oncogenic mutations in human biology — and MAVIS calls all
three structurally silent, correctly. They act by impairing GTP hydrolysis, a
catalytic mechanism, not by destabilising the fold or the RAF1 interface, and the
benchmark grades their binding axis as `neutral` on that basis.

That makes the example a useful check on what MAVIS is and is not. A structural
score near zero is *not* a benign call: mechanism and pathogenicity are separate
axes, and this example is the clearest demonstration of the separation in the
whole benchmark. Anyone who reads the output as "MAVIS thinks G12D is fine" has
misread the tool.

## Running it

```bash
# 1. structures from RCSB (no licence needed)
python examples/kras_craf/fetch_structures.py --out examples/kras_craf/structures

# 2. validate config + variants, stop before FoldX
python run.py --config examples/kras_craf/systems.yaml \
              --variants examples/kras_craf/variants.csv \
              --structures examples/kras_craf/structures \
              --out /tmp/mavis_kras --dry-run

# 3. full scoring (needs your own FoldX 5.x Linux binary)
python run.py --config examples/kras_craf/systems.yaml \
              --variants examples/kras_craf/variants.csv \
              --structures examples/kras_craf/structures \
              --out /tmp/mavis_kras --foldx /path/to/foldx

# 4. compare against the frozen benchmark values
python examples/kras_craf/compare_to_reference.py /tmp/mavis_kras/structural_results.csv
```

Step 2 is the useful smoke test: it exercises config parsing, structure
resolution, and variant expansion — the three things that actually break when
someone sets up their own system.

## Reproducing headline metrics instead

To check the published numbers rather than a single system, use the self-test,
which runs off cached intermediates and needs no structures and no FoldX:

```bash
python verification/verify_stage6.py \
  --intermediate inputs/intermediate/mavis_v7_results_with_nbhd.csv \
  --am inputs/AM_variants_mavis_mechanism_test.xlsx \
  --scripts-dir scripts
```

See `verification/README.md` for what each check covers.

## Tolerance when comparing

FoldX is stochastic; repeated runs of the same mutation differ slightly.
`compare_to_reference.py` sets its tolerance to 3x the largest per-axis replicate
SD recorded in the reference run, so it flags real disagreement rather than
ordinary run-to-run noise. Mechanism strings are compared exactly — a mechanism
flip is a genuine difference even when the underlying dG moved only a little,
which is exactly the borderline case worth surfacing.
