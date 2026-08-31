"""The ``genome xref`` sub-app — identifiers to and from **Gene id stem**s, and symbols.

A thin Typer wrapper over :mod:`genome.xref`: it translates arguments, makes one call on
an :class:`~genome.xref.xref.XrefSet` and renders. It ships from this package so that what
a hop prints and what the answer type holds change in one place. No **Assembly** is named
and no **Genome** is opened by either command — an identifier is a name and not a place.

Examples
--------
>>> from genome.xref.cli import app
>>> [command.name for command in app.registered_commands]
['ids', 'symbols']
"""

from __future__ import annotations

import json as _json

import typer
from typer._click import Command as _ClickCommand
from typer._click import Context as _ClickContext
from typer.core import TyperGroup as _TyperGroup

from genome.xref import NAMESPACES as _NAMESPACES
from genome.xref import SYMBOL as _SYMBOL
from genome.xref import ResolvedStems as _ResolvedStems
from genome.xref import ResolvedSymbols as _ResolvedSymbols
from genome.xref import ResolvedXrefIds as _ResolvedXrefIds
from genome.xref import XrefSet as _XrefSet

#: What a failed xref hop raises, and every one of them already names its next action: a
#: species, a source or a **Namespace** that resolves to nothing is a ``LookupError``; a set
#: that could not be fetched, and a directory an interrupted download left unfinished, are
#: ``RuntimeError``s; a publisher's file that does not match its pin and a stored slice this
#: package did not write are ``ValueError``s; and a file that went away under the read is an
#: ``OSError``. Its own name rather than the annotation sub-app's list reused, because an
#: **Xref set** belongs to no assembly and no annotation — the two lists are alike by
#: coincidence rather than by construction.
_XREF_ERRORS = (ValueError, OSError, RuntimeError, LookupError)

#: What to do when neither direction is named, or both are. Spelled once because the two
#: mistakes are one mistake, and it names both flags because which one is wanted is exactly
#: what the caller has not said.
_XREF_DIRECTION_HELP = (
    "error: name exactly one direction. `--to-stems NAMESPACE` reads the ids as that "
    "namespace and answers in gene id stems; `--from-stems NAMESPACE` reads them as stems "
    "and answers in that namespace. An identifier does not say which system it belongs to, "
    "so nothing here infers the direction from the strings you passed."
)

#: The **Namespace**s ``--to-stems`` offers: every one this package knows except the symbol
#: one, which is not a conversion this command makes at all. Derived rather than written out
#: again, so a namespace added to the package reaches the help without being spelled twice —
#: and so the help cannot go on advertising a conversion the command refuses.
_TO_STEMS_NAMESPACES: tuple[str, ...] = tuple(name for name in _NAMESPACES if name != _SYMBOL)

#: What to do when the symbol namespace is asked *toward* the hub. The reason is the API's
#: own — answering here would match this release's approved spellings and nothing else — but
#: the next action is this surface's to give: the exception names ``match_symbols(symbols)``,
#: which is a Python call and no next step for someone in a shell, so the command that
#: answers it is named instead. Checked before the set is prepared, so the refusal costs no
#: download; the same reason ``_XREF_DIRECTION_HELP`` is checked where the flags are read.
_XREF_SYMBOL_HELP = (
    "error: a gene symbol is not read toward gene id stems with --to-stems. It would match "
    "this release's approved spellings and nothing else, so a table spelling a gene the way "
    "the authority used to would come back unresolved rather than matched — which is what "
    "happens to 31 of EpiFactors' 801 rows. Run `genome xref symbols SPECIES SYMBOL...` "
    "instead, which matches approved, previous and alias spellings and says on each match "
    "which kind it was. `--from-stems symbol` is the labelling direction and is answered "
    "here."
)

#: The columns one **Symbol match** is printed as, which are the keys
#: :meth:`~genome.xref.xref.SymbolMatch.as_json` writes and in its order — so the text
#: rendering and ``--json`` cannot drift apart, and every cell printed is a value the API put
#: in the answer rather than one assembled here. ``kind`` is what makes a spelling the
#: authority retired distinguishable from its current one, which is the whole point of
#: matching a symbol rather than converting it, so it is a column and never a footnote.
_SYMBOL_MATCH_COLUMNS: tuple[str, ...] = ("symbol", "gene_id_stem", "kind")

#: Those three behind the spelling that was asked about, which is the mapping key and so
#: belongs to no match. It earns a column of its own rather than being assumed equal to the
#: match's: folded, ``brca1`` is what was asked and ``BRCA1`` is what the authority spells.
_SYMBOL_COLUMNS: tuple[str, ...] = ("asked", *_SYMBOL_MATCH_COLUMNS)

#: What the flat spelling this group's own name shadows is told, on stderr. ``genome xref
#: SPECIES ID...`` converted identifiers before the tree grew, and the ids it was handed
#: are the very tokens Click resolves a group's subcommand from — so the notice says which
#: command actually ran rather than only that something is deprecated.
_FLAT_SPELLING_NOTICE = (
    "DeprecationWarning: `genome xref SPECIES ID...` is deprecated; it is `genome xref ids "
    "SPECIES ID...` now, and that is what ran. `genome xref symbols SPECIES SYMBOL...` "
    "matches gene symbols."
)


class _XrefGroup(_TyperGroup):
    """The ``xref`` group, whose own name is also a command's old flat spelling.

    Alone among the six sub-apps this one collides: ``genome xref`` named a command before
    it named a group, so the species a caller passes arrives exactly where Click looks for
    a subcommand and the invocation would die on ``No such command 'Homo sapiens'`` — the
    bare refusal the deprecated aliases exist to avoid. A first token that is no subcommand
    of this group is therefore handed to ``ids`` unconsumed, with
    :data:`_FLAT_SPELLING_NOTICE` on stderr so ``--json`` on stdout still parses. Deleted
    with the alias table in :mod:`genome.cli`, in the release after this one.

    Examples
    --------
    >>> from genome.xref.cli import app
    >>> app.info.cls.__name__
    '_XrefGroup'
    """

    def resolve_command(
        self, ctx: _ClickContext, args: list[str]
    ) -> tuple[str | None, _ClickCommand | None, list[str]]:
        """Resolve a subcommand, falling back to ``ids`` on the flat spelling it replaced."""
        if args and self.get_command(ctx, args[0]) is None:
            typer.echo(_FLAT_SPELLING_NOTICE, err=True)
            return "ids", self.get_command(ctx, "ids"), args
        return super().resolve_command(ctx, args)


app = typer.Typer(
    help="Convert identifiers against one published xref set, and match gene symbols.",
    no_args_is_help=True,
    cls=_XrefGroup,
)


@app.command("ids")
def xref_ids(
    species: str = typer.Argument(
        ...,
        help="Species an xref set exists for, e.g. 'Homo sapiens' — the slug "
        "'homo_sapiens' names the same one. A species none exists for names the ones "
        "that do rather than answering nothing.",
    ),
    ids: list[str] = typer.Argument(
        ...,
        help="The identifiers to convert. Each comes back in the order you passed it, "
        "with its version suffix and its namespace's CURIE prefix accepted either way.",
    ),
    to_stems: str | None = typer.Option(
        None,
        "--to-stems",
        metavar="NAMESPACE",
        help=f"Read the ids as this namespace and answer in gene id stems: "
        f"{', '.join(_TO_STEMS_NAMESPACES)}, whichever of them this set carries. A gene "
        f"symbol is not among them — it matches spellings the authority has retired and "
        f"each match carries which kind it was, so `genome xref symbols` answers it. "
        f"Exactly one of this and --from-stems is named.",
    ),
    from_stems: str | None = typer.Option(
        None,
        "--from-stems",
        metavar="NAMESPACE",
        help="Read the ids as gene id stems and answer in this namespace, `symbol` "
        "included — which gives the authority's one current approved spelling, the one a "
        "figure axis wants. A versioned gene id is accepted and reduced to its stem, so an "
        "annotation's own ids go straight in.",
    ),
    source: str | None = typer.Option(
        None,
        "--source",
        help="Answer from this xref source rather than the species' default one. Which "
        "publisher answers is a scientific choice and not a detail: NCBI and Ensembl "
        "agree on 57.6% of human gene-level (GeneID, ENSG) pairs.",
    ),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Convert identifiers to and from gene id stems against one published xref set.

    The way a column of Entrez GeneIDs from a GEO series, UniProt accessions from a
    mass-spec run or HGNC ids from a curated resource reaches this package's answers,
    without writing Python and without the hand-built join everyone in the lab writes
    slightly differently. No assembly is named and no genome is opened: an identifier is a
    name and not a place.

    **The direction is named, never inferred.** `--to-stems NAMESPACE` reads the ids as
    that namespace and answers in gene id stems; `--from-stems NAMESPACE` reads them as
    stems and answers in that namespace. A string does not say which system it belongs to,
    so `HGNC:11998` asked the wrong way answers *nothing found* rather than quietly turning
    around. There is no third direction: Entrez to HGNC is two calls and the join is yours,
    which keeps the hop visible in your pipeline rather than invisible in ours.

    **A gene symbol is the one namespace these two directions do not mirror.**
    `--from-stems symbol` is answered here and gives the authority's single current approved
    spelling. The other way round is `genome xref symbols`, because a symbol also matches
    spellings the authority has retired and each match carries which kind it was — so
    `--to-stems symbol` exits 2 naming that command rather than matching approved spellings
    alone, which is what drops 31 of EpiFactors' 801 rows.

    **The pairs go to stdout, tab-separated, so the output pipes** — `cut -f2` is the
    answer, `cut -f1` says what asked for it — and the heading, the publisher's URL and the
    counts go to stderr. An id naming two genes prints two rows rather than whichever came
    first, and **an id that resolved to nothing gets a row too, with an empty second
    column**: what your list holds and this release does not is the one thing a hand-rolled
    join drops silently. `--json` carries the same answer, keyed by what was asked about,
    with those ids under `unresolved`.

    Omitting `--source` answers from the species' default xref source, so everyone in the
    lab reaches for the same one without discussing it. It is a default and not a
    recommendation: naming a source is how the scientific choice gets made deliberately,
    and every answer names the source and the release that produced it either way. **A
    default is per species and per question**, so `--from-stems symbol` — the one question
    here that is about symbols — is answered by the source that carries them, `hgnc` for
    human and `alliance_bgi` for mouse and worm, rather than by the identifier default.

    Naming a species prepares its set, which the first time is a download. **The lab's CPU
    cluster compute nodes have no internet**, so a set must be constructed once from a login
    node — by running this there, or from Python — before a job that needs it is submitted;
    after that it is read from the **Data dir** and shared by every project on the machine.

    Exits with code 2 when no direction is named or both are, and when the symbol namespace
    is asked toward the hub; and with code 1 when no set exists for the species — the message
    names the ones that do — when the source is not one this package prepares, when the set
    is not here and cannot be fetched, when the namespace is not one the set carries, and
    when a directory holds a set left unfinished.
    """
    # Which way the hop goes and which namespace it is in arrive on one flag, so neither can
    # be given without the other. Naming both, or neither, is the one thing left to check.
    if to_stems is not None and from_stems is None:
        to_hub, namespace = True, to_stems
    elif from_stems is not None and to_stems is None:
        to_hub, namespace = False, from_stems
    else:
        typer.echo(_XREF_DIRECTION_HELP, err=True)
        raise typer.Exit(code=2)

    # The one direction-and-namespace pair this command does not answer, refused where the
    # flags are read rather than where the set is prepared: it is another command's question
    # however the set turns out, so nothing is downloaded to find that out.
    if to_hub and namespace.strip().lower() == _SYMBOL:
        typer.echo(_XREF_SYMBOL_HELP, err=True)
        raise typer.Exit(code=2)

    # Which default an unnamed source is filled in with is the question's to decide, and
    # the question is on the flags — so the namespace goes to the constructor that decides,
    # and nothing here chooses between two of them.
    try:
        xrefs = _XrefSet.for_namespace(species, namespace, source, progressbar=not json)
        answer = xrefs.to_stems(ids, namespace) if to_hub else xrefs.from_stems(ids, namespace)
    except _XREF_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(answer.as_json()))
        return
    _report_xref(answer, source_url=xrefs.source_url, to_hub=to_hub)


def _report_xref(
    answer: _ResolvedStems | _ResolvedXrefIds, *, source_url: str, to_hub: bool
) -> None:
    """Print one hop as the pairs it found, with what answered and what it missed beside them.

    The pairs to stdout and everything else to stderr, for the reason `annotation gene-list`
    splits them: two tab-separated columns are what a shell pipeline wants, and which
    publisher asserted them is what a reader wants. **Every id asked about gets at least one
    row**, the ones that resolved to nothing getting one with an empty second column, so
    nothing leaves this command shorter than it arrived.

    Which noun each column holds is the surface's to say and not the answer's: the answer
    carries the namespace and the direction is what the caller named, so it is passed in
    rather than guessed back out of which type came through.
    """
    asked_as = f"{answer.namespace} ids" if to_hub else "gene id stems"
    answered_as = "gene id stems" if to_hub else f"{answer.namespace} ids"
    rows = [(asked, found) for asked, values in answer.resolved.items() for found in values]
    typer.echo(
        f"{asked_as} -> {answered_as} for {answer.species} ({answer.source} {answer.release})",
        err=True,
    )
    typer.echo(f"  source  {source_url}", err=True)
    typer.echo(
        f"  {len(answer.resolved)} resolved, {len(rows)} {answered_as}, "
        f"{len(answer.unresolved)} this release names none for",
        err=True,
    )
    for asked, found in rows:
        typer.echo(f"{asked}\t{found}")
    for asked in answer.unresolved:
        typer.echo(f"{asked}\t")


@app.command("symbols")
def match_symbols(
    species: str = typer.Argument(
        ...,
        help="Species an xref set exists for, e.g. 'Homo sapiens' — the slug "
        "'homo_sapiens' names the same one. The species fixes the authority, so a symbol "
        "is matched against that authority's spellings and no other's.",
    ),
    symbols: list[str] = typer.Argument(
        ...,
        help="The gene symbols to match, answered in the order you passed them. "
        "Surrounding whitespace goes; case does not, unless --case-insensitive is named.",
    ),
    source: str | None = typer.Option(
        None,
        "--source",
        help="Answer from this xref source rather than the species' default one for "
        "symbols, which is `hgnc` for human and `alliance_bgi` for mouse and worm. Naming "
        "one is deliberate and is never overridden: a source that carries no symbol says "
        "so and names the one that does, rather than matching nothing.",
    ),
    case_insensitive: bool = typer.Option(
        False,
        "--case-insensitive",
        help="Fold case on both sides — your spelling and the authority's — and still "
        "answer with every gene matched. Off by default: the species is fixed by the set, "
        "so 'Brca1' asked of a human set is a mouse spelling asked of the wrong authority.",
    ),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Print the genes each gene symbol names, and which kind of spelling matched.

    The way a gene list copied out of a paper becomes usable without first finding its ids —
    and without the join that silently drops every row spelling its gene the way the
    authority used to. No assembly is named and no genome is opened: a symbol is a name and
    not a place.

    **A symbol is matched, never converted.** Approved, previous *and* alias spellings are
    matched, every **Gene id stem** any of them names comes back, and each match says which
    kind of spelling it was — so ambiguity is what you are handed rather than something
    resolved on your behalf. `ADCY3` is HGNC's approved symbol for one gene and a symbol it
    retired from another, and both are printed. This is why it is a command of its own and
    not a third direction of `genome xref ids`: matching approved spellings alone would drop
    exactly the rows this exists for — 31 of EpiFactors' 801 human rows spell their gene the
    way HGNC spelled it years ago. The opposite hop, a stem to the authority's one current
    approved spelling, is `genome xref ids --from-stems symbol`.

    **The matches go to stdout, tab-separated, so the output pipes** — `cut -f3` is the
    answer, `cut -f1` says what asked for it and `cut -f4` says which kind of spelling
    matched — and the heading, the publisher's URL, the counts and what this source could
    not have matched go to stderr. A symbol naming two genes prints two rows rather than
    whichever came first, and **a symbol this release matched nothing for gets a row too,
    with every other column empty**. Column 2 is the authority's own spelling, which is not
    always the one asked about: folded, `brca1` asked comes back as `BRCA1` matched.

    **Matching is exact by default**, because the species is fixed by the set:
    `--case-insensitive` folds both sides and still answers with every gene matched rather
    than picking one.

    **What this source could not have matched is printed too.** Only HGNC publishes previous
    and alias spellings typed; mouse and worm match current approved symbols alone, their
    authorities' typed spellings belonging to publishers that cannot be pinned or cannot be
    fetched. So the answer says which kinds it could match and why the others are
    missing — without which *this gene is not in the release* and *this source does not
    publish the spelling you used* would both be silence.

    **Omitting `--source` answers from the species' default source for symbols**, which is
    not the same row as its default for identifiers: human's identifiers come from
    `alliance`, whose cross-reference file publishes no human symbol at all, and its symbols
    from `hgnc`; mouse's and worm's from `alliance_bgi`. A default is per species and per
    question for that reason, and every answer names the source and release that produced it
    either way. Naming a source is still how the scientific choice gets made deliberately,
    and a named one is never swapped — so `--source alliance` here exits 1 saying that set
    carries no symbol, rather than quietly answering from somebody else's file.

    Naming a species prepares its set, which the first time is a download. **The lab's CPU
    cluster compute nodes have no internet**, so a set must be constructed once from a login
    node — by running this there, or from Python — before a job that needs it is submitted.

    Exits with code 1 when no set exists for the species — the message names the ones that
    do — when the source is not one this package prepares, when a named source carries no
    symbols at all — the message names the one that does — when the set is not here and
    cannot be fetched, and when a directory holds a set left unfinished.
    """
    try:
        xrefs = _XrefSet.for_symbols(species, source, progressbar=not json)
        answer = xrefs.match_symbols(symbols, case_insensitive=case_insensitive)
    except _XREF_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(answer.as_json()))
        return
    _report_symbols(answer, source_url=xrefs.source_url)


def _report_symbols(answer: _ResolvedSymbols, *, source_url: str) -> None:
    """Print one set's symbol matches, with what asserted them and what limits them beside.

    The matches to stdout and everything else to stderr, for the reason `ids` splits them:
    tab-separated columns are what a shell pipeline wants, and which authority asserted them
    is what a reader wants. **Every symbol asked about gets at least one row**, the ones that
    matched nothing getting one with every other column empty, so nothing leaves this command
    shorter than it arrived.

    The limits line prints only when there is one: the line above it already names the kinds
    that *were* matched on, so a set that carries all three has nothing left to explain.
    """
    matching = "case-insensitive" if answer.case_insensitive else "exact"
    rows = [(asked, match) for asked, matches in answer.resolved.items() for match in matches]
    typer.echo(
        f"gene symbols -> gene id stems for {answer.species} ({answer.source} {answer.release})",
        err=True,
    )
    typer.echo(f"  source   {source_url}", err=True)
    typer.echo(f"  columns  {', '.join(_SYMBOL_COLUMNS)}", err=True)
    typer.echo(f"  matching {matching}, on {', '.join(answer.kinds)} spellings", err=True)
    typer.echo(
        f"  {len(answer.resolved)} resolved, {len(rows)} matches, "
        f"{len(answer.unresolved)} this release matched nothing for",
        err=True,
    )
    if answer.limits is not None:
        typer.echo(f"  limits   {answer.limits}", err=True)
    for asked, match in rows:
        written = match.as_json()
        typer.echo("\t".join([asked, *(str(written[column]) for column in _SYMBOL_MATCH_COLUMNS)]))
    for asked in answer.unresolved:
        typer.echo(asked + "\t" * len(_SYMBOL_MATCH_COLUMNS))
