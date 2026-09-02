# What the template's static gates cost, and four measurements that lied

Measured 2026-09-02, while adopting `vale`, `markdownlint` and `conformance` from
`liulab-repo-template` and raising pyright to `standard`.

**Four findings. (1) The gates were cheap and the estimate was not: vale found 5 errors where a
`wc -w` estimate had predicted six files over their word cap, because vale counts parsed prose and
`wc` counts table markup and link syntax — `AGENTS.md` measures 891 words against its cap of 1000,
not the 1297 `wc` reports. (2) Raising pyright to `standard` does not enforce an annotation rule:
it adds five checks over `basic` and leaves both parameter-annotation rules at `none`, so a wholly
unannotated function is silent under both modes. (3) Where a readability metric has a passing
exemplar in the same file, matching the exemplar converges and optimising the metric does not —
comparable effort moved a changelog 12.46 → 12.41 one way and 12.46 → 10.82 the other, while
cutting 3% of the words and *raising* the entry count. (4) Four separate measurements taken during
this work returned plausible numbers that answered a different question than the one asked; every
one was a tool that accepts a narrowing argument and silently widens it.**

This is a measurement, not a decision. What it decided lives in
[ADR-0024](../adr/0024-the-template-is-this-repos-authority.md) and
[ADR-0025](../adr/0025-python-is-3-13-only.md).

## Method

Run on the maintainer's box against the repository at `chore/template-sync`. Tool versions are
those the lockfile pins: vale 3.19.0, markdownlint-cli2 0.23.2 (markdownlint 0.41.1), pyright
1.1.411. Every count below is from the tool's own summary line, not from a wrapper. Where a
measurement is a comparison, the variants were run back to back on one tree.

## 1. What the gates actually found

| Gate | Errors on first run | After |
| --- | --- | --- |
| `vale` | 5 in 67 files | 0 |
| `markdownlint` | 966 in 40 files | 0 |
| `conformance` | passes, 1 waiver | unchanged |

The vale errors were four `Lab.Jargon` and one `Lab.Readability`. **No length cap was breached at
any point**, which contradicted the estimate this work was planned against.

The estimate came from `wc -w`, which reported `AGENTS.md` at 1297 words against a cap of 1000 and
five glossary files between 1355 and 2661. Vale's own count is 891 for `AGENTS.md`. The difference
is that vale parses the markdown and counts prose, while `wc` counts table pipes, code spans, link
targets and emphasis markers as words. A repository whose agent-facing documents are table-heavy —
this one — is measured far more harshly by `wc` than by the rule that actually gates it.

The 966 markdownlint errors were mostly one convention: 649 `MD060` table-column-style, 142 `MD049`
emphasis-style, 130 `MD010` hard tabs. The linter's own fixer resolved 937 of them. The residue was
16 bare code fences, one delimiter row matching no style, two identical headings in a research
note, and 4 `MD024` duplicate-heading errors discussed below.

Two consequences of the mass reformat worth recording, because neither is visible in the counts:
five `console` blocks lost their `$` prompts to `MD014` and were retyped as `bash`, and the
tab-to-space conversion reached pasted `chrom.sizes` and `.tsv` output in four documents. Names and
numbers survived; those blocks are no longer literally tab-delimited.

## 2. pyright's `standard` mode does not check annotations

The plan for raising pyright assumed `standard` was the mechanism behind this repository's
full-annotation rule. It is not. Measured by extracting the rule tables from the pyright bundle and
confirmed with a probe file:

| Mode | `def f(x): return x` |
| --- | --- |
| `basic` | 0 errors |
| `standard` | 0 errors |
| `strict` | 3 errors, and 556 across this tree |

`standard` adds exactly five checks over `basic`: possibly-unbound variable, incompatible method
override, incompatible variable override, overlapping overload, and function member access. Both
`reportMissingParameterType` and `reportUnknownParameterType` stay at `none`.

Naming those two rules beside the mode is what enforces the rule. Their cost here, and the cost of
the third member of the same family, measured on this tree with the config isolated in an otherwise
empty directory:

| Configuration | Errors |
| --- | --- |
| `standard` alone | 0 |
| `+ reportMissingParameterType` | 0 |
| `+ reportUnknownParameterType` | 0 |
| `+ reportUnknownVariableType` | 0 |

All four are free here, but the last is not free everywhere: the same file, `scripts/conformance.py`
with its untyped `yaml` import, costs 20 to 25 errors in the template repository. It is free here
only because `scripts/` sits in `extraPaths` rather than `include`, and **`extraPaths` resolves
without checking** — a directory listed there contributes zero diagnostics and is
indistinguishable from a directory that passes.

The two rules also differ in kind, which matters for anyone adopting them elsewhere.
`reportMissingParameterType` is syntactic and local: `Any` satisfies it and nothing outside the
repository can turn it red. `reportUnknownParameterType` is semantic and transitive: it fails code
that *is* annotated, such as a bare `list` or `Callable`, and its verdict depends on whether every
dependency ships complete types. This repository holds the transitive rule because `pixi.lock` is
committed and CI resolves locked, so an upstream release cannot change the verdict without a
deliberate re-lock that somebody reviews.

## 3. Matching an exemplar beats optimising the metric

`Lab.Readability` is a whole-file Flesch-Kincaid average capped at grade 11. `CHANGELOG.md` failed
at 12.46. Splitting the file at its first released-version heading:

| Section | Words | Grade |
| --- | --- | --- |
| Whole file | 21,290 | 12.46 |
| Unreleased only | 16,592 | 13.81 |
| Released history only | 4,730 | passes |

The shipped history passes on its own, so nothing was being dragged down by accumulation — the
unreleased half was the entire cause, and it was 3.5 times the size of every shipped release
combined. That ruled out the theory that the cap was mis-scoped for a document of this type: a cap
that only becomes inapplicable once the writing degrades is not a scoping argument.

Two approaches were then tried, with comparable effort:

| Approach | Grade |
| --- | --- |
| Shortening sentences to move the number | 12.46 → 12.41 |
| Rewriting to match the released entries' shape | 12.46 → 10.82 |

The second cut 663 words, three percent of the file, and *raised* the entry count from 158 to 161.
So the text was barely shortened; it was resegmented. The metric is a ratio over structure, and
shaving words inside sentences already committed to changes numerator and denominator together.
An exemplar encodes where sentences end, which is the variable being measured.

The margin after the rewrite is 0.18 grades, over a file that only grows. This cap is therefore a
rule about how the next entry is written, not a state the tree can be put into once.

## 4. Four measurements that returned plausible wrong numbers

Every one of these produced a number in the range a person would have guessed, which is why none
announced itself.

| What was run | What it silently did instead |
| --- | --- |
| `pixi run vale <file>` on three split files | Ran from the manifest root and graded the repository's real file three times — three identical grades across three different word counts |
| `pixi run pyright -p <config>` with the config in a scratch directory | Took that directory as the project root and swept up unrelated probe files, attributing 8 phantom errors to the repository |
| The same, fixed by using absolute include paths | Still scanned the config's own directory; the fix is a directory holding nothing but the config |
| `pixi run vale <file>` in the template repository | Appended the filename to a task that already ended in `.`, grading all 18 files while reporting one |

The shared shape is a tool that accepts a narrowing argument and silently widens it. The defence
that would have caught all four is the same: predict the number before running, and treat
"different inputs, identical outputs" or a suspiciously round zero as a tooling failure until shown
otherwise.

A fifth artefact of the same family is not a measurement at all. The comment above
`conformance.py`'s loader states that it exits 2 when it cannot list tracked files. It exits 1,
which is also what a rule violation exits, so "the checker never ran" and "this repository violates
a rule" are indistinguishable to anything reading the code:

| Condition | Exit |
| --- | --- |
| No `.git` to list tracked files from | 1 |
| Rules failing | 1 |
| Clean repository | 0 |

That comment is worse than an unenforced claim, and the reason generalises. An unenforced claim is
visibly a claim. A wrong one arrives in the form of an answer, and an answer is what stops you
looking — the same mechanism as a plausible number. Both were read as reassurance here before
either was tested.

## What this suggests about where these defects live

Three defects in the shared configuration surfaced during this work: pyright's `standard` mode not
covering the annotation rule, the docs generator not validating the navigation, and `MD024` firing
on a correctly formatted changelog. All three are age-dependent. One needs code somebody stopped
annotating, one needs a page to have moved, and one needs a second release.

The template's own conformance job proves that a *fresh* repository passes, which is exactly the
class of repository that cannot exercise any of them. The blind spot is structural rather than a
run of three coincidences, and the remedy is not machinery — an artificially aged fixture would be
a fake, and no rule can check for defects nobody has thought of. The remedy is this audit
direction, an older repository read against the template, repeated. Filed as an issue it would be
closed by whoever did one audit; the value is entirely in the repetition.

The `MD024` and pyright findings were reproduced independently in `liulab-repo-template` by the
session working there, which also found the navigation and exit-code defects. Half of what is
recorded above is theirs.
