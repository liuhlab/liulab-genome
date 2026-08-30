"""Tests for the ``genome tf`` sub-app."""

from __future__ import annotations

import json as _json
from pathlib import Path

from genome.annotation import AnnotationRegistry, GtfAnnotation
from genome.assembly import metadata
from genome.cli import app
from genome.tf.cofactor import CofactorTable, cofactor_table
from genome.tf.gene import TFGeneTable, tf_gene_table

from .._cli import output, runner

#: The two assemblies a census ships for, with the annotation each one defaults to. Named
#: rather than derived: which species has a census is exactly what the TF command is about,
#: so the pairing is written down where the test controls it. Two of them because the
#: censuses are two publishers' and their columns differ — mouse's has none beyond the
#: uniform four — and a command that only ever met Lambert's would not know that.
_TF_ASSEMBLY, _TF_ANNOTATION = "hg38", "gencode_v50"
_MOUSE_ASSEMBLY, _MOUSE_ANNOTATION = "mm39", "gencode_vM39"


def _register(assembly: str, assembly_dir: Path, gtf: Path, name: str) -> GtfAnnotation:
    """Register ``gtf`` under ``assembly_dir``, so a command has something to report on."""
    return AnnotationRegistry.locate(assembly, assembly_dir).register_path(gtf, name)


def _gtf_declaring(*gene_ids: str) -> str:
    """A GTF declaring one gene, transcript and exon per id, in the shape GENCODE has.

    The ids are the only thing that varies: nothing reading it looks at a coordinate, so
    every gene sits on the same interval rather than inviting anyone to.
    """
    return "".join(
        f'chrI\ttest\t{feature}\t1\t100\t.\t+\t.\tgene_id "{gene_id}"; transcript_id "{gene_id}_t";\n'
        for gene_id in gene_ids
        for feature in ("gene", "transcript", "exon")
    )


class TestTFGeneListCommand:
    """``genome tf gene-list`` — an assembly's TF genes, shaped like the annotation one.

    The shipped censuses answer, as the shipped curated lists answer for gene categories:
    which genes are transcription factors is the census's judgement and no fixture stands
    in for it, so every expectation below is read off the shipped file. What is asserted
    here is the command and not the crossing — ``tests/annotation/`` owns that — so: the
    stdout/stderr split that makes the output pipe, the record ``--json`` emits, and a
    non-zero exit naming the next action for each of the three ways it can fail.
    """

    def _census(self, assembly: str) -> TFGeneTable:
        """The census shipped for ``assembly``'s species, whichever publisher wrote it."""
        species = metadata.assembly_metadata(assembly).species
        assert species is not None, f"{assembly} has no species in the assembly table"
        census = tf_gene_table(species)
        assert census is not None, f"no census ships for {species}"
        return census

    def _registered(
        self,
        liulab_data: Path,
        *gene_ids: str,
        assembly: str = _TF_ASSEMBLY,
        name: str = _TF_ANNOTATION,
    ) -> None:
        """Register a GTF declaring ``gene_ids`` as ``assembly``'s ``name``, where it lives."""
        source = liulab_data / f"{assembly}.{name}.gtf"
        source.write_text(_gtf_declaring(*gene_ids))
        _register(assembly, liulab_data / "genome" / assembly, source, name)

    def _positive(self, assembly: str, count: int) -> list[str]:
        """``count`` gene ids, one per assessed-positive stem, versioned as GENCODE spells them."""
        return [f"{stem}.1" for stem in self._census(assembly).assessed_positive[:count]]

    def test_ids_naming_and_json_work_for_human_and_a_mouse_census_answers_empty(
        self, liulab_data: Path
    ) -> None:
        gene_ids = self._positive(_TF_ASSEMBLY, 2)
        self._registered(liulab_data, *gene_ids)

        result = runner.invoke(app, ["tf", "gene-list", _TF_ASSEMBLY])
        assert result.exit_code == 0
        assert result.stdout == "".join(f"{gene_id}\n" for gene_id in gene_ids)
        # Whose judgement it is must be printed and must not cost the pipe, so it goes
        # beside the ids: the heading, the census's own attribution, and what the crossing
        # cost — the stems the census holds that this annotation carries no gene for.
        assert f"{_TF_ASSEMBLY} / {_TF_ANNOTATION}" in result.stderr
        assert "Homo sapiens" in result.stderr
        assert self._census(_TF_ASSEMBLY).provenance.attribution() in result.stderr
        unresolved = len(self._census(_TF_ASSEMBLY).assessed_positive) - len(gene_ids)
        assert f"2 genes, 2 gene ids, {unresolved} stems" in result.stderr

        named_ids = self._positive(_TF_ASSEMBLY, 1)
        self._registered(liulab_data, *named_ids, name="mine")
        named = runner.invoke(app, ["tf", "gene-list", _TF_ASSEMBLY, "--annotation", "mine"])
        assert named.exit_code == 0
        assert named.stdout.splitlines() == named_ids
        assert f"{_TF_ASSEMBLY} / mine" in named.stderr

        # The record `--json` emits carries the genes, the provenance and the stems this
        # annotation carries no gene for, not dropped.
        census = self._census(_TF_ASSEMBLY)
        stem = census.assessed_positive[0]
        self._registered(liulab_data, f"{stem}.1", name="json")
        json_result = runner.invoke(
            app, ["tf", "gene-list", _TF_ASSEMBLY, "--annotation", "json", "--json"]
        )
        assert json_result.exit_code == 0
        payload = _json.loads(json_result.stdout)
        assert list(payload) == [
            "assembly",
            "annotation",
            "species",
            "provenance",
            "genes",
            "gene_ids",
            "unresolved",
        ]
        assert (payload["assembly"], payload["annotation"]) == (_TF_ASSEMBLY, "json")
        assert payload["gene_ids"] == [f"{stem}.1"]
        assert payload["provenance"]["pubmed_id"] == census.provenance.pubmed_id
        assert payload["unresolved"]
        cells = dict(
            zip(census.columns, census.rows[census.gene_id_stems.index(stem)], strict=True)
        )
        gene = payload["genes"][0]
        assert gene["dbd_family"] == cells["dbd_family"]
        assert gene["judgements"]["tf_assessment"] == cells["tf_assessment"]

        # The verdict travels with the census that reached it, so which one spoke is a
        # fact about the assembly's species and never about which one came first.
        mouse_ids = self._positive(_MOUSE_ASSEMBLY, 1)
        self._registered(liulab_data, *mouse_ids, assembly=_MOUSE_ASSEMBLY, name=_MOUSE_ANNOTATION)
        mouse = runner.invoke(app, ["tf", "gene-list", _MOUSE_ASSEMBLY])
        assert mouse.exit_code == 0
        assert mouse.stdout.splitlines() == mouse_ids
        assert "Mus musculus" in mouse.stderr
        assert self._census(_MOUSE_ASSEMBLY).provenance.attribution() in mouse.stderr
        assert self._census(_TF_ASSEMBLY).provenance.publisher not in mouse.stderr

        # AnimalTFDB ships the four uniform columns and no more, so every mouse gene's
        # judgements are empty. A surface reaching for a **TF assessment** that is not
        # there would raise here rather than print, which is why nothing does.
        mouse_json = runner.invoke(app, ["tf", "gene-list", _MOUSE_ASSEMBLY, "--json"])
        assert mouse_json.exit_code == 0
        mouse_gene = _json.loads(mouse_json.stdout)["genes"][0]
        assert mouse_gene["judgements"] == {}
        assert mouse_gene["dbd_family"]  # the uniform four are there all the same

    def test_an_unregistered_annotation_an_unsupported_species_or_an_unknown_assembly_exits_one(
        self, liulab_data: Path
    ) -> None:
        unregistered = runner.invoke(app, ["tf", "gene-list", _TF_ASSEMBLY])
        assert unregistered.exit_code == 1
        assert unregistered.stdout == ""
        assert f"genome annotation register {_TF_ASSEMBLY} {_TF_ANNOTATION}" in output(unregistered)

        # Human gene ids registered for a worm assembly: the species is the assembly's
        # own and never what the GTF happens to hold, so this is refused rather than
        # answered.
        self._registered(
            liulab_data, *self._positive(_TF_ASSEMBLY, 1), assembly="ce11", name="wormbase_ws298"
        )
        no_census = runner.invoke(app, ["tf", "gene-list", "ce11"])
        assert no_census.exit_code == 1
        assert no_census.stdout == ""
        assert "no TF census ships" in output(no_census)
        assert "Caenorhabditis elegans" in output(no_census)
        assert "Homo sapiens" in output(no_census)  # …and what may be asked about instead

        # Not the same fact as no census ships: the question was which species this is,
        # and no row answered it. An unlisted local key is the ordinary way in.
        self._registered(
            liulab_data, *self._positive(_TF_ASSEMBLY, 1), assembly="tiny", name="mine"
        )
        unknown_species = runner.invoke(app, ["tf", "gene-list", "tiny", "--annotation", "mine"])
        assert unknown_species.exit_code == 1
        assert unknown_species.stdout == ""
        assert "nothing says what species 'tiny' is" in output(unknown_species)
        assert "Homo sapiens" in output(unknown_species)


class TestTFCofactorListCommand:
    """``genome tf cofactor-list`` — an assembly's cofactors, shaped like ``tf gene-list``.

    The shipped tables answer, as the shipped censuses answer for TF genes: which genes a
    publisher lists as cofactors is that publisher's judgement and no fixture stands in for
    it, so every expectation below is read off the shipped file. What is asserted here is
    the command and not the crossing — ``tests/annotation/`` owns that — so: the
    stdout/stderr split that makes the output pipe, the record ``--json`` emits, and a
    non-zero exit naming the next action for each of the three ways it can fail.

    One test is about neither: a worm assembly answers here while ``tf gene-list`` refuses
    the same one, because a publisher assessed worm cofactors and none has released a worm
    TF census. That asymmetry is the publishers' shape rather than a defect, and it is
    pinned so that nobody smooths it away.
    """

    def _table(self, assembly: str) -> CofactorTable:
        """The cofactor table shipped for ``assembly``'s species, whichever publisher wrote it."""
        species = metadata.assembly_metadata(assembly).species
        assert species is not None, f"{assembly} has no species in the assembly table"
        table = cofactor_table(species)
        assert table is not None, f"no cofactor table ships for {species}"
        return table

    def _registered(
        self,
        liulab_data: Path,
        *gene_ids: str,
        assembly: str = _MOUSE_ASSEMBLY,
        name: str = _MOUSE_ANNOTATION,
    ) -> None:
        """Register a GTF declaring ``gene_ids`` as ``assembly``'s ``name``, where it lives."""
        source = liulab_data / f"{assembly}.{name}.gtf"
        source.write_text(_gtf_declaring(*gene_ids))
        _register(assembly, liulab_data / "genome" / assembly, source, name)

    def _listed(self, assembly: str, count: int) -> list[str]:
        """``count`` gene ids, one per listed stem, versioned as GENCODE spells them."""
        return [f"{stem}.1" for stem in self._table(assembly).cofactor_stems[:count]]

    def test_only_the_gene_ids_reach_stdout_the_annotation_may_be_named_and_json_carries_it_all(
        self, liulab_data: Path
    ) -> None:
        gene_ids = self._listed(_MOUSE_ASSEMBLY, 2)
        self._registered(liulab_data, *gene_ids)

        result = runner.invoke(app, ["tf", "cofactor-list", _MOUSE_ASSEMBLY])
        assert result.exit_code == 0
        assert result.stdout == "".join(f"{gene_id}\n" for gene_id in gene_ids)
        # Whose list it is must be printed and must not cost the pipe, so it goes beside
        # the ids: the heading, the publishers' own attribution, and what the crossing
        # cost — the stems the table holds that this annotation carries no gene for.
        assert f"{_MOUSE_ASSEMBLY} / {_MOUSE_ANNOTATION}" in result.stderr
        assert "Mus musculus" in result.stderr
        assert self._table(_MOUSE_ASSEMBLY).provenance.attribution() in result.stderr
        unresolved = len(self._table(_MOUSE_ASSEMBLY).cofactor_stems) - len(gene_ids)
        assert f"2 cofactors, 2 gene ids, {unresolved} stems" in result.stderr

        named_ids = self._listed(_MOUSE_ASSEMBLY, 1)
        self._registered(liulab_data, *named_ids, name="mine")
        named = runner.invoke(app, ["tf", "cofactor-list", _MOUSE_ASSEMBLY, "--annotation", "mine"])
        assert named.exit_code == 0
        assert named.stdout.splitlines() == named_ids
        assert f"{_MOUSE_ASSEMBLY} / mine" in named.stderr

        # The record `--json` emits carries the cofactors, the provenance and the stems
        # this annotation carries no gene for, not dropped.
        table = self._table(_MOUSE_ASSEMBLY)
        stem = table.cofactor_stems[0]
        self._registered(liulab_data, f"{stem}.1", name="json")
        json_result = runner.invoke(
            app, ["tf", "cofactor-list", _MOUSE_ASSEMBLY, "--annotation", "json", "--json"]
        )
        assert json_result.exit_code == 0
        payload = _json.loads(json_result.stdout)
        assert list(payload) == [
            "assembly",
            "annotation",
            "species",
            "provenance",
            "cofactors",
            "gene_ids",
            "unresolved",
        ]
        assert (payload["assembly"], payload["annotation"]) == (_MOUSE_ASSEMBLY, "json")
        assert payload["gene_ids"] == [f"{stem}.1"]
        # One provenance record per publisher that contributed, never one flattened row.
        assert [source["pubmed_id"] for source in payload["provenance"]["sources"]] == [
            source.pubmed_id for source in table.provenance.sources
        ]
        assert payload["unresolved"]
        cells = dict(zip(table.columns, table.rows[table.gene_id_stems.index(stem)], strict=True))
        cofactor = payload["cofactors"][0]
        assert (cofactor["symbol"], cofactor["source"]) == (cells["symbol"], cells["source"])
        assert cofactor["is_cofactor"] is True
        assert cofactor["classifications"]["animaltfdb_category"] == cells["animaltfdb_category"]

    def test_a_worm_assembly_answers_although_tf_gene_list_refuses_the_same_one(
        self, liulab_data: Path
    ) -> None:
        # AnimalTFDB assessed worm cofactors and no publisher has censused worm
        # transcription factors, so one command answers and the other refuses for one
        # registered annotation. WormBase's ids carry no version, and a stem carrying
        # none resolves to itself, so they are registered exactly as the table spells them.
        gene_ids = list(self._table("ce11").cofactor_stems[:2])
        self._registered(liulab_data, *gene_ids, assembly="ce11", name="wormbase_ws298")

        answered = runner.invoke(app, ["tf", "cofactor-list", "ce11"])
        refused = runner.invoke(app, ["tf", "gene-list", "ce11"])

        assert answered.exit_code == 0
        assert answered.stdout.splitlines() == gene_ids
        assert "Caenorhabditis elegans" in answered.stderr
        assert refused.exit_code == 1
        assert "no TF census ships" in output(refused)

    def test_an_unregistered_annotation_an_unsupported_species_or_an_unknown_assembly_exits_one(
        self, liulab_data: Path
    ) -> None:
        unregistered = runner.invoke(app, ["tf", "cofactor-list", _MOUSE_ASSEMBLY])
        assert unregistered.exit_code == 1
        assert unregistered.stdout == ""
        next_action = f"genome annotation register {_MOUSE_ASSEMBLY} {_MOUSE_ANNOTATION}"
        assert next_action in output(unregistered)

        # Mouse gene ids registered for a yeast assembly: the species is the assembly's
        # own and never what the GTF happens to hold, so this is refused rather than
        # answered.
        self._registered(
            liulab_data, *self._listed(_MOUSE_ASSEMBLY, 1), assembly="sacCer3", name="ensgene_v101"
        )
        no_table = runner.invoke(app, ["tf", "cofactor-list", "sacCer3"])
        assert no_table.exit_code == 1
        assert no_table.stdout == ""
        assert "no cofactor table ships" in output(no_table)
        assert "Saccharomyces cerevisiae" in output(no_table)
        assert "Mus musculus" in output(no_table)  # …and what may be asked about instead

        # Not the same fact as no table ships: the question was which species this is, and
        # no row answered it. The message says which shipped table could not be chosen, so
        # a cofactor question is never refused with a sentence about a census.
        self._registered(
            liulab_data, *self._listed(_MOUSE_ASSEMBLY, 1), assembly="tiny", name="mine"
        )
        unknown_species = runner.invoke(
            app, ["tf", "cofactor-list", "tiny", "--annotation", "mine"]
        )
        assert unknown_species.exit_code == 1
        assert unknown_species.stdout == ""
        assert "nothing says what species 'tiny' is, so no cofactor table" in output(
            unknown_species
        )
        assert "Mus musculus" in output(unknown_species)
