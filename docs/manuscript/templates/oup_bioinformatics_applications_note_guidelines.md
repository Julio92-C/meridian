# OUP Bioinformatics — Applications Note submission guidelines

Compiled 2026-05-29 from
- https://academic.oup.com/bioinformatics/pages/author-guidelines
- https://academic.oup.com/bioinformatics/pages/instructions_for_authors
- https://academic.oup.com/bioinformatics/pages/submission_online
- https://academic.oup.com/bioinformatics/pages/scope_guidelines

## Scope and form factor

An Applications Note is a short description of:

- novel software,
- a new algorithm implementation,
- a database, or
- a web service / network interface.

Notes cannot describe trivial utilities or tools that require
significant installation effort. The application name should appear
in the title.

## Length

> "Up to 4 pages; this is approx. 2,600 words or 2,000 words plus one
> figure."

So: ~2,600 words text-only, or ~2,000 words + 1 figure.

Our current `applications_note.md` runs ~900–1,000 words in the body
(excluding refs and front matter) plus one figure — well under
the limit, with room to add the per-domain PCoA / Sankey /
utility-engine details flagged in `manuscript_edit_plan.md`.

## Abstract structure (de-facto standard in this journal)

The OUP author-guidelines page itself does NOT explicitly enumerate
abstract sub-headings for Applications Notes, but every modern
Applications Note in this journal uses the same five labels:

- **Motivation:** problem statement / gap.
- **Results:** what the tool does and what it achieves.
- **Availability and Implementation:** where the software lives, how
  to install / invoke it, license.
- **Contact:** corresponding author's email.
- **Supplementary information:** (optional) pointer to supplementary
  PDF / data.

The OUP guidelines explicitly state that:

> "Additional supplementary data can be published online-only by the
> journal. This supplementary material should be referred to in the
> abstract of the Application Note."

Target total abstract length: ~150 words (this is what L-ARRAP,
rANOMALY, and other similar notes use).

## Software availability — hard requirements

The author guidelines spell out:

- Software must be "freely available to non-commercial users".
- "Availability and Implementation must be clearly stated in the
  article."
- Software must remain available for "a full two years following
  publication."
- Web services must require no mandatory user registration.
- Implementation details must be clearly documented.

Our note already meets this: MIT licence, GitHub URL with
`Rscript run_pipeline.R config.yaml` invocation line, dependency
pinning via `renv`.

## Reference style

The OUP guidelines page surveyed today does not enumerate the citation
format for Applications Notes specifically. The journal uses a
**numbered (Vancouver-style) citation system in the body, with a
numbered reference list**, matching what `applications_note.md`
already does ([1], [2], …).

No explicit maximum number of references is documented on the public
guidelines page, but Applications Notes in this journal typically
have 10–20 references. We have 13 — appropriate.

## Templates

The OUP guidelines state:

> "Initial submissions are not required to follow the template. However,
> because manuscripts are published online following acceptance, it is
> important that any revised version of the manuscript supplied by the
> author follow the provided templates."

Word and LaTeX templates are linked from the Author guidelines page
(behind the publisher's submission portal — exact URLs were not
surfaced in the public author-guidelines text we fetched). Initial
submission of `applications_note.md` (or its `.docx` / `.pdf` export)
is acceptable without templating; the template should be applied
after acceptance.

## Backup journals (in order of fit)

Already enumerated in `references.md`:

1. **Bioinformatics Advances** — OUP sister journal, more flexible format.
2. **F1000Research** — rolling open review (rANOMALY's venue).
3. **GigaScience** — rewards reproducibility, FAIR data, containers.
4. **BMC Bioinformatics** — Software section.
5. **Briefings in Bioinformatics** — if scope expands (L-ARRAP's venue).

## What is still ambiguous

- The author-guidelines page we fetched does not enumerate the abstract
  word limit or reference cap explicitly. The 150-word abstract and
  10–20 reference range are inferred from peer Applications Notes in
  this journal.
- We were not able to download the LaTeX / Word template directly
  via WebFetch (Silverchair watermarked links require a session
  token). Pull the templates from the submission portal after
  account creation, or from the OUP Bioinformatics author-guidelines
  page in a browser.
