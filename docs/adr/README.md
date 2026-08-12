# Architecture decision records

One decision per file, one paragraph: `# N. <the decision as a claim>` and one to three sentences on
what was there, what was decided and why the obvious reading lost. Twelve lines is the ceiling, and a
`**Status.**` line only where a record was amended or superseded. The bar for writing a new one at all
is stated in [`AGENTS.md`](../../AGENTS.md).

Every record here is system-wide. This repository has no per-context ADR directories, so this
directory is the whole tree.

**Numbers are permanent and gaps stay.** Take the next number not present anywhere in the tree, never
reuse one. A number that was taken stays taken — including by a record that was deleted, superseded,
or reserved on a branch that never landed — so a gap is expected and is never filled.

Cite a record by number (`ADR-0003`), never by path: records move and numbers do not.
