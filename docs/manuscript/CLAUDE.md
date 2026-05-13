# Bioinformatics Applications Note — drafting workspace

This directory bootstraps a Claude Code session whose sole job is to draft
a short Applications Note manuscript about the Metagenomics Pipeline
Automation tool. The pipeline source lives one directory up at `..`.

## Read these first

1. `brief.md` — full task spec (structure, length, journal style)
2. `pipeline_facts.md` — pre-extracted facts about the pipeline so you do
   not need to re-read 12 R modules to write the Methods section
3. `references.md` — target journal + template papers with DOIs
4. `memory_snapshot.md` — relevant project memory (alpha state, ground-truth
   validation context)

## Output

Single file: `applications_note.md` in this directory.

## What this session IS for

- Drafting the manuscript
- Fetching the L-ARRAP and rANOMALY template papers if WebFetch is allowed
- Writing factual claims grounded in the facts file (no fabricated numbers)

## What this session is NOT for

- Editing pipeline R scripts under `../R/`
- Running the pipeline
- Modifying any config or settings

Those tasks belong to the parent project's main Claude session.

## Permissions note

The previous attempt was blocked because the `Write` tool was denied.
If that happens again, ask the user to grant Write permission for this
directory rather than working around it.
