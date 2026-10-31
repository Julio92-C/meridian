# Project memory snapshot — drafting context

Relevant entries from the parent project's auto-memory, copied here so
this isolated session has the context without crossing project
boundaries.

## Pipeline state (project memory)

The pipeline is at **Phase 1 / alpha**. The R/00–R/12 scaffold has been
smoke-tested end-to-end against real data but each stage is still being
validated against ground-truth outputs. Reverse-engineering per stage
from the hand-edited reference scripts is the current focus.

**Application to drafting**: be honest in the Conclusion — frame as
"we present" not "we provide a production-ready". Use language like
"in active development" or "alpha release". Do not promise stability
guarantees.

## Validation dataset (reference memory)

Real reference data and ground-truth scripts live at:
`PC_JC_2024-11-29_Julio_Gallus/` (chicken gut microbiome, first batch).
The dataset includes ~46 hand-edited R scripts and known-good output
figures/tables. Each pipeline stage is validated by comparing its
output against the corresponding hand-edited script.

**Application to drafting**: Section 3 (Application / Validation) should
describe this benchmark setup. Cite the specific stages where validation
has been completed. Use `{{TODO: insert benchmark table}}` for the
quantitative comparison table — the user will fill numbers from the
actual validation runs.

## Known divergence from ground truth (project memory)

The pipeline's R/02 fixes a taxid-parse bug present in the ground-truth
`cleanData.R` (variant-suffix taxids were lost in the original). This
means the pipeline produces ~+110 (sample, GENE) rows compared to the
ground-truth output. This is a **correctness improvement, not a bug**.

**Application to drafting**: if the validation table includes resistome
counts, flag this as a deliberate fix in the discussion. Do not present
the row-count difference as a failure-to-reproduce.

## Drafting principle (feedback memory, adapted)

Per-project flexibility is a load-bearing design choice — never
hardcode sample IDs, control names, or column slicing in any module
code. The manuscript's Implementation section should highlight this
config-driven flexibility as a feature, since it directly contrasts
with how the reference / ground-truth scripts had to be hand-edited
for every new study.
