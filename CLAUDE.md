# MERIDIAN — repo-level rules

This file is read at session start by any Claude Code instance launched
from the repo root. Subdirectory CLAUDE.md files (e.g.
`docs/manuscript/CLAUDE.md`) add scoped rules on top of these.

## `applications_note.docx` is hand-formatted — do NOT re-render

`docs/manuscript/applications_note.docx` carries manual Word formatting
(track-changes state, custom paragraph styles, figure sizing, layout
tweaks) that pandoc cannot round-trip. **Never run `pandoc ... -o
applications_note.docx`** — every prior re-render silently stripped
Julio's styling, costing ~70 kB of embedded format per pass.

Propagation rules:

- `.md` → `.docx`: **don't**. The `.docx` is the deliverable, not a
  render target. If you really need a downstream `.docx` from the
  current `.md`, use the `tracked-diff` skill (writes
  `<w:ins>`/`<w:del>` OOXML directly into a docx zip and preserves
  most underlying formatting).
- `.docx` → `.md`: yes, on demand. Workflow:
  1. `pandoc docs/manuscript/applications_note.docx -t gfm --wrap=none -o _tmp.md`
  2. `diff -u docs/manuscript/applications_note.md _tmp.md` — extract
     SEMANTIC changes only; ignore pandoc style noise (escape
     backslashes, `<>` autolinks, HTML `<figure>` wrappers, table
     column-spec padding, list-item indentation).
  3. Apply each semantic change to `applications_note.md` via targeted
     `Edit` calls, preserving the `.md`'s own styling (backticks
     around code identifiers, multi-line affiliations, image-alt-text
     captions).
  4. Delete `_tmp.md` and commit.

## Safeguards

- A pre-commit hook at `.githooks/pre-commit` aborts any commit where a
  staged `.docx` shrinks by more than 8% vs `HEAD` (override the
  threshold with `DOCX_SHRINK_PCT=20 git commit`). Bypass entirely for
  a deliberate large-figure removal with `git commit --no-verify`.
- The hook is enabled via `git config core.hooksPath .githooks` (set
  once per clone; re-run after a fresh clone).
- Local Word-formatted snapshots of the manuscript belong in
  `docs/manuscript/_backup/` (gitignored). Drop a copy there before any
  risky editing pass.
