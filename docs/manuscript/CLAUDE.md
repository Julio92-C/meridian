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

- `applications_note.md` — drafting source under version control. Edit
  freely.
- `applications_note.docx` — **hand-formatted** Word deliverable. See
  the rule below.

## `.docx` is source-of-truth, NOT a render target

`applications_note.docx` carries manual Word formatting (track-changes
acceptance state, custom paragraph styles, figure sizing, layout tweaks)
that pandoc cannot round-trip. **Do NOT regenerate `.docx` from `.md`
via `pandoc -o applications_note.docx`** — every prior re-render
silently clobbered Julio's styling.

When Julio reviews the manuscript in Word, sync changes the OTHER way
(.docx → .md):

1. `pandoc applications_note.docx -t gfm --wrap=none -o _tmp.md`
2. `diff -u applications_note.md _tmp.md` — identify SEMANTIC changes;
   ignore pandoc-style noise (escape backslashes, `<>` autolinks, HTML
   `<figure>` wrappers, table column-spec width, list-item indentation)
3. Apply each semantic change to `applications_note.md` via targeted
   `Edit` calls, preserving the `.md`'s own styling (backticks around
   code identifiers, multi-line affiliations, image-alt-text captions)
4. Delete `_tmp.md` and commit

If a downstream step genuinely needs a fresh `.docx`, use the
`tracked-diff` skill (writes `<w:ins>`/`<w:del>` OOXML directly into the
docx zip) rather than a full pandoc render — it preserves far more of
the underlying formatting.

## Local baselines: `_backup/`

Drop a Word-formatted snapshot into `docs/manuscript/_backup/` before
any risky editing pass. The directory is gitignored — never committed
to history, always available locally for restore.

A repo-root pre-commit hook (`.githooks/pre-commit`) aborts commits
where a staged `.docx` shrinks by more than 8% vs `HEAD`, catching the
exact pandoc-clobber failure mode. Override the threshold with
`DOCX_SHRINK_PCT=20 git commit`, or bypass entirely for a deliberate
large-figure removal with `git commit --no-verify`.

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
