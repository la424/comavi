# Worked examples

Each subdirectory is a complete, self-contained COMAVI input set with the frozen
benchmark values it should reproduce. Use these to verify an installation, to see
what a valid config and variant file look like before writing your own, or as a
template to copy.

Everything except the final scoring step runs **without a FoldX licence**:
structures come from RCSB and the AlphaFold DB, and `run.py --dry-run` validates
the config and expands the variant list without invoking FoldX.

| Example | System | Structures | What it demonstrates |
|---|---|---|---|
| `hemoglobin_dimer/` | HBB - HBA1 (alpha1-beta1) | 2HHB x-ray 1.74 A + AlphaFold monomers | HBB E6V (sickle cell): a famous pathogenic variant that is **structurally silent** on all three axes |

## Why hemoglobin E6V is the example to start with

HBB Glu6Val is the sickle-cell variant — one of the best-characterised pathogenic
mutations in human medicine. COMAVI scores it as structurally silent on all three
axes, and the benchmark grades all three expectations as `neutral`. That is the
correct answer: E6V does not destabilise the globin fold or the alpha1-beta1
interface. It creates a *new* hydrophobic surface patch that drives polymerisation
of deoxy-HbS into fibres — an intermolecular gain-of-function that a single
alpha-beta dimer cannot express, and that no ddG on these three axes is designed
to detect.

So the example is a check on what COMAVI is and is not. **A structural score near
zero is not a benign call.** Mechanism and pathogenicity are separate axes, and
E6V is the cleanest demonstration in the benchmark: maximally pathogenic,
structurally quiet. Anyone reading the output as "COMAVI thinks sickle cell is
fine" has misread the tool.

## Running it

```bash
# 1. structures: 2HHB chains A+B from RCSB, HBB/HBA1 monomers from AlphaFold DB
python examples/hemoglobin_dimer/fetch_structures.py \
       --out examples/hemoglobin_dimer/structures

# 2. validate config + variants, stop before FoldX
python run.py --config examples/hemoglobin_dimer/systems.yaml \
              --variants examples/hemoglobin_dimer/variants.csv \
              --structures examples/hemoglobin_dimer/structures \
              --out /tmp/comavi_hb --dry-run

# 3. full scoring (needs your own FoldX 5.x binary)
python run.py --config examples/hemoglobin_dimer/systems.yaml \
              --variants examples/hemoglobin_dimer/variants.csv \
              --structures examples/hemoglobin_dimer/structures \
              --out /tmp/comavi_hb --foldx /path/to/foldx

# 4. compare against the frozen benchmark values
python examples/hemoglobin_dimer/compare_to_reference.py \
       /tmp/comavi_hb/structural_results.csv
```

Step 2 is the useful smoke test: it exercises config parsing, structure
resolution, and variant expansion — the three things that actually break when
someone sets up their own system.

## The two structure sources are not interchangeable

The complex axis runs on **crystal chains** (2HHB A+B, a measured interface, so no
pLDDT gating applies to the binding axis). The monomer axis runs on **AlphaFold
models**, which is what `monomer_pdb` in the config points at.

The config's `monomer_offset: -1` is calibrated to AlphaFold numbering, which keeps
the initiator Met that mature globin numbering drops — so CSV position 6 maps to
model residue 7 = Glu. Crystal chain B residue 7 is *also* a glutamate, so
substituting crystal chains for the monomer models would score **Glu7 instead of
Glu6** and return a plausible-looking but wrong result. `fetch_structures.py` gets
both sources right and asserts the residue identity at the mutated position before
you spend FoldX time on it.

## Reproducing headline metrics instead

To check the published numbers rather than a single system, use the self-test,
which runs off cached intermediates and needs no structures and no FoldX:

```bash
python verification/verify_stage6.py \
  --intermediate inputs/intermediate/comavi_v7_results_with_nbhd.csv \
  --am inputs/AM_variants_comavi_mechanism_test.xlsx \
  --scripts-dir scripts
```

See `verification/README.md` for what each check covers.

## Tolerance when comparing

FoldX is stochastic; repeated runs of the same mutation differ slightly.
`compare_to_reference.py` sets its tolerance to 3x the largest per-axis replicate
SD in the reference run, with a 0.15 kcal/mol floor (one axis here has SD 0.0 from
a single run, and a zero tolerance would flag ordinary noise). Mechanism strings
are compared exactly — a mechanism flip is a genuine difference even when the
underlying ddG moved only a little, which is the borderline case worth surfacing.
