"""Tests for genome.annotation.registry — the four states, and what may be asked of them.

``AnnotationRegistry`` and the three scans it settles at construction, the **Default
annotation** rule, the status report, the gene categories an annotation declares, and the
same questions addressed by assembly name rather than opened.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pytest

from genome.annotation import GtfAnnotation
from genome.annotation.curated import (
    GeneCategoryNotDeclaredError,
    GeneListAssemblyMismatchError,
    NoGeneCategoriesError,
    curated_gene_list,
)
from genome.annotation.registration import (
    ChromosomeMismatchError,
    MergeSource,
    annotation_dir,
    annotation_register_command,
    register_merged_gtf,
)
from genome.annotation.registry import (
    AnnotationNotRegisteredError,
    AnnotationRegistry,
    AnnotationStatus,
    AnnotationStatusRow,
    annotation_status,
    default_annotation,
    gene_list,
    gene_lists,
    list_annotations,
    list_broken_annotations,
)
from genome.assembly.registration import AssemblyDir
from genome.store.completion import UnfinishedRegistrationError, record_path, work_dir

from ..conftest import FakeFetch
from .conftest import (
    _CHIMERA,
    _CURATED,
    _CURATED_ASSEMBLY,
    _FOOD,
    _FOOD_COMPONENT,
    _GTF,
    _NAME,
    _WORM,
    _WORM_COMPONENT,
    _register_by_path,
    _row,
    _write_chrom_sizes,
)


class TestListAnnotations:
    """Registered means a record that agrees with disk — never a database file."""

    def test_registered_means_a_record_agreeing_with_disk_never_just_a_database_file(
        self, tmp_path: Path
    ) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        annotation = _register_by_path(tmp_path, src, "finished")
        assert list(list_annotations(tmp_path)) == ["finished"]

        # A database file with no record beside it (a build killed part-way) is not
        # registered...
        halfway = annotation_dir(tmp_path, "halfway")
        halfway.mkdir(parents=True)
        (halfway / "halfway.db").write_bytes(b"half a database")
        assert list(list_annotations(tmp_path)) == ["finished"]

        # ...and neither is a record that no longer agrees with what is on disk.
        annotation.db.write_bytes(b"truncated")
        assert list_annotations(tmp_path) == {}


class TestListBrokenAnnotations:
    """The complement of ``list_annotations``: what is on disk and cannot be trusted.

    Every directory under ``gtf/`` is registered, broken, or not begun; this reports
    the middle one, which is otherwise invisible to anything that lists.
    """

    def test_broken_means_untrusted_files_never_a_fresh_or_finished_directory(
        self, tmp_path: Path
    ) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)

        _register_by_path(tmp_path, src, "no-record")
        record_path(annotation_dir(tmp_path, "no-record")).unlink()
        broken = list_broken_annotations(tmp_path, "tiny")
        assert list(broken) == ["no-record"]
        assert broken["no-record"].directory == annotation_dir(tmp_path, "no-record")
        assert "holds files but no .completion.json" in broken["no-record"].problem

        disagreeing = _register_by_path(tmp_path, src, "disagreeing")
        disagreeing.db.write_bytes(b"truncated")
        broken = list_broken_annotations(tmp_path, "tiny")
        assert "disagreeing.db" in broken["disagreeing"].problem

        # A finished annotation is never reported broken...
        _register_by_path(tmp_path / "finished", src, "mine")
        assert list_broken_annotations(tmp_path / "finished", "tiny") == {}

        # ...and neither is a fresh or empty one: ADR-0007 says an absent or empty
        # directory is a fresh registration, not a broken one. A run interrupted before
        # it downloaded anything must not be reported.
        annotation_dir(tmp_path, "fresh").mkdir(parents=True)
        work_dir(annotation_dir(tmp_path, "fresh")).mkdir()
        assert "fresh" not in list_broken_annotations(tmp_path, "tiny")

        # And an assembly with no annotations at all answers empty rather than raising.
        assert list_broken_annotations(tmp_path / "no-such-directory", "sacCer3") == {}

    def test_a_broken_one_leaves_the_others_listed(self, tmp_path: Path) -> None:
        # The invariant: one broken annotation must not stop the rest being found, and
        # nothing here raises — reporting is the whole point.
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        _register_by_path(tmp_path, src, "healthy")
        _register_by_path(tmp_path, src, "damaged")
        record_path(annotation_dir(tmp_path, "damaged")).unlink()

        assert list(list_annotations(tmp_path)) == ["healthy"]
        assert list(list_broken_annotations(tmp_path, "tiny")) == ["damaged"]

    def test_the_repair_it_names_depends_on_whether_the_table_or_a_path_can_be_shown(
        self, tmp_path: Path
    ) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)

        # A name the table offers is repaired by name.
        _register_by_path(tmp_path, src, "ensgene_v101")
        record_path(annotation_dir(tmp_path, "ensgene_v101")).unlink()
        offered = list_broken_annotations(tmp_path, "sacCer3")
        assert offered["ensgene_v101"].repair == (
            "genome annotation register sacCer3 ensgene_v101 --force"
        )
        assert offered["ensgene_v101"].repair in offered["ensgene_v101"].problem

        # An unlisted one is repaired from the path its own record remembers.
        unlisted = _register_by_path(tmp_path, src, "mine")
        unlisted.db.write_bytes(b"truncated")
        broken = list_broken_annotations(tmp_path, "tiny")
        assert broken["mine"].repair == f"genome annotation register-gtf tiny {src} mine --force"

        # No record survives to say which GTF it was built from, so there is no path to
        # print: the command is named with the one thing it still needs filled in,
        # rather than a path that would not run.
        _register_by_path(tmp_path, src, "unknowable")
        record_path(annotation_dir(tmp_path, "unknowable")).unlink()
        broken = list_broken_annotations(tmp_path, "tiny")
        assert (
            broken["unknowable"].repair
            == "genome annotation register-gtf tiny <path> unknowable --force"
        )

        # A record survives, but the path it names is gone — same placeholder, for the
        # same reason: the path it remembers would not run either.
        gone = _register_by_path(tmp_path, src, "gone")
        gone.db.write_bytes(b"truncated")
        src.unlink()
        broken = list_broken_annotations(tmp_path, "tiny")
        assert str(src) not in broken["gone"].repair
        assert broken["gone"].repair == "genome annotation register-gtf tiny <path> gone --force"


class TestDefaultAnnotation:
    """The one rule that decides a default, wherever the question is asked from."""

    def test_the_flag_the_sole_registration_or_an_explicit_choice_decides_it(self) -> None:
        assert default_annotation([_row()], []) == _NAME
        # Everyone in the lab reaches for the same one, whatever this machine happens
        # to hold...
        assert default_annotation([_row()], ["something_else"]) == _NAME
        # ...unless a caller explicitly chooses another, which wins over the flag.
        assert default_annotation([_row()], [_NAME], explicit="something_else") == "something_else"

        # With nothing flagged, the sole registered annotation wins, or there is none.
        assert default_annotation([_row(default=False)], ["only_one"]) == "only_one"
        assert default_annotation([], []) is None
        assert default_annotation([], ["one", "two"]) is None


class TestAnnotationStatus:
    """What an assembly's table offers, set against what is registered on this machine."""

    def test_a_fresh_machine_reports_the_shipped_table_serialized_and_untouched(
        self, fake_fetch: FakeFetch, liulab_data: Path
    ) -> None:
        # The case it most needs to serve: a fresh machine, where the answer is
        # entirely the shipped table's.
        payload = annotation_status("sacCer3")

        assert payload.assembly == "sacCer3"
        assert payload.directory == liulab_data / "genome" / "sacCer3"
        assert payload.default_annotation == "ensgene_v101"
        rows = payload.annotations
        assert [(r.name, r.offered, r.registered) for r in rows] == [("ensgene_v101", True, False)]
        assert rows[0].provider == "UCSC"
        assert rows[0].path is None

        # `--json` is this report rendered, so a row's fields and the payload's keys are
        # one spelling: a surface reads attributes and never names a key of its own.
        assert payload.as_json() == {
            "assembly": "sacCer3",
            "directory": str(liulab_data / "genome" / "sacCer3"),
            "default_annotation": "ensgene_v101",
            "annotations": [asdict(row) for row in payload.annotations],
        }

        # What the closing line of `genome annotation list` needs: the default's own state,
        # so "not registered here" and "broken here" are told apart by the report itself.
        default = payload.default_row
        assert default is payload.annotations[0]
        assert default is not None
        assert not default.registered

        # No default decided, so no row is about one — the state a fresh unlisted
        # assembly is in, and not an error.
        nothing = annotation_status("tiny")
        assert nothing.default_annotation is None
        assert nothing.default_row is None

        # And asking creates nothing and fetches nothing.
        annotation_status("hg38")
        assert fake_fetch.calls == []
        assert not (liulab_data / "genome" / "hg38").exists()

    def test_a_registered_or_broken_annotation_is_reported_whether_offered_or_not(
        self, tmp_path: Path
    ) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        offered_dir, unlisted_dir = tmp_path / "offered", tmp_path / "unlisted"
        annotation = _register_by_path(offered_dir, src, "ensgene_v101")
        _register_by_path(unlisted_dir, src, "mine")

        offered = annotation_status("sacCer3", cache_dir=offered_dir)
        unlisted = annotation_status("tiny", cache_dir=unlisted_dir)

        offered_rows = offered.annotations
        assert [(r.name, r.offered, r.registered) for r in offered_rows] == [
            ("ensgene_v101", True, True)
        ]
        assert offered_rows[0].path == str(annotation.gtf)

        unlisted_rows = unlisted.annotations
        assert [(r.name, r.offered, r.registered) for r in unlisted_rows] == [("mine", False, True)]
        assert unlisted_rows[0].provider is None
        assert unlisted_rows[0].broken is False
        assert unlisted_rows[0].problem is None
        assert unlisted.default_annotation == "mine"  # nothing flagged, and it is alone

        # The bug this closes: half-registered and never-fetched looked identical here.
        # No row lists the second one and no record vouches for it either, so nothing
        # used to mention it at all — a broken annotation is reported as broken, and
        # never as simply absent, whether it is offered or not.
        record_path(annotation_dir(offered_dir, "ensgene_v101")).unlink()
        broken_unlisted = _register_by_path(unlisted_dir, src, "gone")
        broken_unlisted.db.write_bytes(b"truncated")

        broken_offered = annotation_status("sacCer3", cache_dir=offered_dir)
        broken_unlisted_payload = annotation_status("tiny", cache_dir=unlisted_dir)

        broken_offered_rows = broken_offered.annotations
        assert [(r.name, r.offered, r.registered, r.broken) for r in broken_offered_rows] == [
            ("ensgene_v101", True, False, True)
        ]
        assert broken_offered_rows[0].repair == (
            "genome annotation register sacCer3 ensgene_v101 --force"
        )
        assert "holds files but no .completion.json" in str(broken_offered_rows[0].problem)
        assert broken_offered_rows[0].path is None

        unlisted_rows = [row for row in broken_unlisted_payload.annotations if row.name == "gone"]
        assert [(r.name, r.offered, r.registered, r.broken) for r in unlisted_rows] == [
            ("gone", False, False, True)
        ]
        assert unlisted_rows[0].repair == f"genome annotation register-gtf tiny {src} gone --force"

    def test_one_broken_annotation_does_not_hide_the_others(self, tmp_path: Path) -> None:
        src = tmp_path / "ann.gtf"
        src.write_text(_GTF)
        assembly_dir = tmp_path / "asm"
        _register_by_path(assembly_dir, src, "healthy")
        _register_by_path(assembly_dir, src, "damaged")
        record_path(annotation_dir(assembly_dir, "damaged")).unlink()

        payload = annotation_status("tiny", cache_dir=assembly_dir)

        rows = payload.annotations
        assert [(r.name, r.registered, r.broken) for r in rows] == [
            ("damaged", False, True),
            ("healthy", True, False),
        ]
        # A broken one is not a registered one, so it never becomes the sole-registered
        # default either.
        assert payload.default_annotation == "healthy"


class TestAnnotationRegistry:
    """One assembly's annotations, and the four states each of them can be in.

    The registry is the one place *registered / broken / offered / not begun* is
    assembled: everything that used to rebuild that four-way state — a genome opening,
    the status report, the error a name nobody registered earns — asks this instead.
    """

    def _registered(self, tmp_path: Path, name: str) -> Path:
        """Register the fixture GTF under ``name`` and return the source it came from."""
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(tmp_path, source, name)
        return source

    def test_the_four_states_are_settled_in_one_construction_bound_to_its_own_directory(
        self, tmp_path: Path
    ) -> None:
        self._registered(tmp_path, "ensgene_v101")
        self._registered(tmp_path, "damaged")
        record_path(annotation_dir(tmp_path, "damaged")).unlink()

        registry = AnnotationRegistry.locate("sacCer3", tmp_path)

        assert registry.registered == ["ensgene_v101"]
        assert [entry.name for entry in registry.broken] == ["damaged"]
        assert [record.name for record in registry.offered] == ["ensgene_v101"]
        assert registry.default == "ensgene_v101"

        # The assembly dir travels with the registry rather than being re-derived from
        # the data root at each question.
        elsewhere = tmp_path / "elsewhere"
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(elsewhere, source, "mine")

        assert AnnotationRegistry.locate("tiny", elsewhere).registered == ["mine"]
        assert AnnotationRegistry.locate("tiny", tmp_path / "unused").registered == []

        # Its status is exactly what the status report answers, and nothing is
        # created by asking.
        assert registry.status() == annotation_status("sacCer3", cache_dir=tmp_path)
        nothing = AnnotationRegistry.locate("sacCer3", tmp_path / "unused" / "genome" / "sacCer3")
        assert nothing.registered == []
        assert nothing.broken == []
        assert nothing.default == "ensgene_v101"
        assert not (tmp_path / "unused" / "genome").exists()

    def test_path_resolves_registered_names_and_explains_unregistered_or_broken_ones(
        self, tmp_path: Path
    ) -> None:
        self._registered(tmp_path, "mine")
        registry = AnnotationRegistry.locate("tiny", tmp_path)
        assert registry.path("mine") == annotation_dir(tmp_path, "mine") / "mine.gtf"

        offering = AnnotationRegistry.locate("sacCer3", tmp_path)
        with pytest.raises(AnnotationNotRegisteredError) as unregistered:
            offering.path("no_such_annotation")
        message = str(unregistered.value)
        assert "no_such_annotation" in message
        assert "ensgene_v101" in message  # what the table does offer

        # Not the command that would itself raise and demand --force: the one that works.
        self._registered(tmp_path, "ensgene_v101")
        record_path(annotation_dir(tmp_path, "ensgene_v101")).unlink()
        with pytest.raises(AnnotationNotRegisteredError) as broken:
            AnnotationRegistry.locate("sacCer3", tmp_path).path("ensgene_v101")
        assert "genome annotation register sacCer3 ensgene_v101 --force" in str(broken.value)

    def test_registering_adopts_the_result_and_the_default_is_the_flag_unless_overridden(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        by_path = AnnotationRegistry.locate("tiny", tmp_path / "path")
        assert by_path.registered == []

        annotation = by_path.register_path(source, "mine")

        assert by_path.registered == ["mine"]
        assert by_path.path("mine") == annotation.gtf
        assert by_path.default == "mine"  # nothing flagged, and it is alone

        fake_fetch.serve("tiny.gtf.gz")
        by_name = AnnotationRegistry.locate("tiny", tmp_path / "name")

        fetched = by_name.register(_NAME, progressbar=False, metadata=_row())

        assert by_name.registered == [_NAME]
        assert by_name.path(_NAME) == fetched.gtf

        flagged = AnnotationRegistry.locate("sacCer3", tmp_path / "flagged")
        assert flagged.default == "ensgene_v101"  # the table's flag

        flagged.register_path(source, "mine")
        assert flagged.default == "ensgene_v101"  # never displaced by a registration

        explicit = AnnotationRegistry.locate("sacCer3", tmp_path / "explicit", default="mine")
        assert explicit.default == "mine"
        assert explicit.registered == []  # need not be registered to win

    def test_a_broken_directory_either_names_a_repair_command_or_is_repaired_by_force(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        source = self._registered(tmp_path, "mine")
        record_path(annotation_dir(tmp_path, "mine")).unlink()
        registry = AnnotationRegistry.locate("tiny", tmp_path)
        assert [entry.name for entry in registry.broken] == ["mine"]

        registry.register_path(source, "mine", force=True)
        assert registry.broken == []
        assert registry.registered == ["mine"]

        # A registry always knows its assembly, so the repair a half-built directory
        # names is a command a shell can run rather than a Python call with the
        # assembly left to guess at.
        directory = annotation_dir(tmp_path, "WS298")
        directory.mkdir(parents=True)
        (directory / "WS298.db").write_bytes(b"half a database")
        with pytest.raises(UnfinishedRegistrationError) as excinfo:
            registry.register_path(source, "WS298")
        assert f"genome annotation register-gtf tiny {source} WS298 --force" in str(excinfo.value)

        # Chrom.sizes defaults to the assembly's own...
        chrom_dir = tmp_path / "chrom"
        _write_chrom_sizes(chrom_dir, "chrI", "chrII", "chrIII")
        with pytest.raises(ChromosomeMismatchError):
            AnnotationRegistry.locate("tiny", chrom_dir).register_path(
                data_dir / "ensembl_style.gtf", _NAME
            )

        # ...but what a Genome opened somewhere of its own passes — the chrom.sizes it
        # actually prepared, rather than the one the layout would name for its
        # assembly — is the one checked against.
        elsewhere = tmp_path / "elsewhere"
        sizes = _write_chrom_sizes(elsewhere, "chrI", "chrII", "chrIII", assembly="elsewhere")
        with pytest.raises(ChromosomeMismatchError):
            AnnotationRegistry(
                AssemblyDir.locate("tiny", elsewhere), chrom_sizes=sizes
            ).register_path(data_dir / "ensembl_style.gtf", _NAME)


# ---------------------------------------------------------------------------------------
# Which genes are in a category, and what an annotation that cannot say answers
# ---------------------------------------------------------------------------------------


def _declared(annotation: str) -> tuple[str, ...]:
    """The categories a shipped curated list declares, in file order — never a fixed set.

    Which categories exist is the curated list's to say and differs per annotation, so
    every assertion below reads them off the shipped file rather than naming them.
    """
    listed = curated_gene_list(annotation)
    assert listed is not None, f"no curated gene list ships for {annotation}"
    return tuple(listed.categories)


def _curated_ids(annotation: str, category: str) -> tuple[str, ...]:
    """The gene ids a shipped curated list puts in ``category``."""
    listed = curated_gene_list(annotation)
    assert listed is not None
    return listed.categories[category].gene_ids


def _register_merged(
    assembly_dir: Path,
    name: str,
    source: Path,
    *,
    assembly: str = _CHIMERA,
    food: str = _FOOD,
) -> GtfAnnotation:
    """Register a merged annotation of the worm/food pair under ``name``, from one GTF.

    A real :func:`register_merged_gtf`, not a hand-written record: what a merge writes
    into ``details`` is what the gene-category path reads back, so standing that in would
    test the stand-in. The fixture GTF is merged twice under two component names, which is
    all the record needs to carry — nothing here reads the features.

    ``food`` names the second contributor's annotation, so a caller may merge in one no
    curated list ships for and test what a contributor that cannot contribute does.
    """
    chrom_sizes = _write_chrom_sizes(
        assembly_dir,
        *(f"chrI__{component}" for component in (_WORM_COMPONENT, _FOOD_COMPONENT)),
        assembly=assembly,
    )
    return register_merged_gtf(
        assembly_dir,
        name,
        [
            MergeSource(_WORM_COMPONENT, _WORM, source),
            MergeSource(_FOOD_COMPONENT, food, source),
        ],
        separator="__",
        chrom_sizes=chrom_sizes,
    )


class TestGeneList:
    """``AnnotationRegistry.gene_list`` — the genes one annotation puts in one category.

    The shipped curated lists answer here, since which categories exist is data and no
    fixture may pretend otherwise. What is asserted is structure: who contributed, what
    order the ids arrive in, and — the point of the whole surface — that an annotation
    which cannot answer raises rather than answering with nothing.
    """

    def _registry(self, tmp_path: Path, *, assembly: str = _CURATED_ASSEMBLY) -> AnnotationRegistry:
        """Register the fixture GTF under a curated annotation's name and open the registry."""
        tmp_path.mkdir(parents=True, exist_ok=True)
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(tmp_path, source, _CURATED, assembly=assembly)
        return AnnotationRegistry.locate(assembly, tmp_path)

    def test_a_plain_annotation_answers_with_one_unattributed_source_and_the_curated_ids(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path)
        category = _declared(_CURATED)[0]

        answer = registry.gene_list(category, _CURATED)

        assert (answer.assembly, answer.annotation, answer.category) == (
            _CURATED_ASSEMBLY,
            _CURATED,
            category,
        )
        # One contributor, and `component` is None rather than the assembly's own name:
        # attribution is what a merge needs, and there is nothing here to attribute.
        assert [source.component for source in answer.sources] == [None]
        assert [source.annotation for source in answer.sources] == [_CURATED]
        # The two sentences travel with the ids, because they are what says whether
        # these ids mean what the caller's metric needs.
        assert answer.sources[0].description.strip()
        assert answer.sources[0].source.strip()
        assert answer.gene_ids == list(_curated_ids(_CURATED, category))
        assert answer.gene_ids  # never an empty answer: a declared category has genes

    def test_the_two_kinds_of_absence_are_told_apart_and_each_names_whats_available(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path)

        with pytest.raises(GeneCategoryNotDeclaredError) as declared:
            registry.gene_list("no_such_category", _CURATED)
        message = str(declared.value)
        assert "no_such_category" in message
        for category in _declared(_CURATED):
            assert category in message

        # The fact #111 exists for: *no categories are declared* and *this category is
        # not declared* are different, and neither is an empty answer.
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(tmp_path, source, "mine")
        nothing_declared_registry = AnnotationRegistry.locate("tiny", tmp_path)
        with pytest.raises(NoGeneCategoriesError) as nothing_declared:
            nothing_declared_registry.gene_list("rRNA", "mine")
        nothing_message = str(nothing_declared.value)
        assert "mine" in nothing_message
        assert _CURATED in nothing_message  # …and which annotations do declare categories

        assert isinstance(declared.value, GeneCategoryNotDeclaredError)
        assert isinstance(nothing_declared.value, NoGeneCategoriesError)
        assert not isinstance(declared.value, NoGeneCategoriesError)
        assert not isinstance(nothing_declared.value, GeneCategoryNotDeclaredError)

    def test_naming_no_annotation_asks_the_default_and_an_unregistered_or_missing_one_names_it(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path / "default")
        category = _declared(_CURATED)[0]

        assert registry.default == _CURATED
        assert registry.gene_list(category).annotation == _CURATED

        # Callers pass names, never paths, and a name nothing registered is resolved by
        # the same `path` every other question goes through — so the message is the
        # same one.
        unregistered = AnnotationRegistry.locate(_CURATED_ASSEMBLY, tmp_path / "unregistered")
        with pytest.raises(AnnotationNotRegisteredError) as excinfo:
            unregistered.gene_list("rRNA", _CURATED)
        assert f"genome annotation register {_CURATED_ASSEMBLY} {_CURATED}" in str(excinfo.value)

        no_default = AnnotationRegistry.locate("tiny", tmp_path / "no-default")
        with pytest.raises(ValueError, match="annotation") as no_default_excinfo:
            no_default.gene_list("rRNA")
        assert "default_gtf" in str(no_default_excinfo.value)

        # A name is unique only within its assembly, so a list found by name alone is not
        # yet known to be about this reference — and answering would hand back another
        # species' genes under this one's name.
        wrong_assembly = self._registry(tmp_path / "wrong", assembly="tiny")

        with pytest.raises(GeneListAssemblyMismatchError) as excinfo:
            wrong_assembly.gene_list(_declared(_CURATED)[0], _CURATED)

        message = str(excinfo.value)
        assert "tiny" in message
        assert _CURATED_ASSEMBLY in message

    def test_gene_lists_returns_every_declared_category_or_raises_if_none_are(
        self, tmp_path: Path
    ) -> None:
        answers = self._registry(tmp_path).gene_lists(_CURATED)
        assert [answer.category for answer in answers] == list(_declared(_CURATED))
        assert all(answer.gene_ids for answer in answers)

        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(tmp_path, source, "mine")
        with pytest.raises(NoGeneCategoriesError):
            AnnotationRegistry.locate("tiny", tmp_path).gene_lists("mine")


class TestMergedGeneList:
    """A **Merged annotation**'s genes, attributed to the component each came from.

    The case #111 was opened over: a chimera of a worm and the bacterium it eats, whose
    ribosomal RNA is both species' and must not arrive as one number.
    """

    def _registry(self, tmp_path: Path, name: str = f"{_WORM}+{_FOOD}") -> AnnotationRegistry:
        """Register a merged annotation under ``name`` and open the chimera's registry."""
        tmp_path.mkdir(parents=True, exist_ok=True)
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_merged(tmp_path, name, source)
        return AnnotationRegistry.locate(_CHIMERA, tmp_path)

    def test_each_contributor_answers_for_its_own_component_and_ids_are_never_deduped(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path)
        shared = next(category for category in _declared(_WORM) if category in _declared(_FOOD))

        answer = registry.gene_list(shared, f"{_WORM}+{_FOOD}")

        assert [(source.component, source.annotation) for source in answer.sources] == [
            (_WORM_COMPONENT, _WORM),
            (_FOOD_COMPONENT, _FOOD),
        ]
        assert answer.annotation == f"{_WORM}+{_FOOD}"
        assert answer.assembly == _CHIMERA
        assert answer.gene_ids == [
            *_curated_ids(_WORM, shared),
            *_curated_ids(_FOOD, shared),
        ]

    def test_a_contributor_that_cannot_answer_is_left_out_rather_than_raised_over(
        self, tmp_path: Path
    ) -> None:
        # One component's annotation is curated and the other's is not, so only one of the
        # two can contribute. That is an omission and not a failure: the answer carries the
        # contributor that can answer, and stays attributed to its component.
        gtf = tmp_path / "ann.gtf"
        gtf.write_text(_GTF)
        _register_merged(tmp_path, "merged", gtf, food="nobody_curated_this")
        registry = AnnotationRegistry.locate(_CHIMERA, tmp_path)
        category = _declared(_WORM)[0]

        answer = registry.gene_list(category, "merged")

        assert [contributed.component for contributed in answer.sources] == [_WORM_COMPONENT]
        assert answer.gene_ids == list(_curated_ids(_WORM, category))

    def test_a_category_no_contributor_declares_raises_and_gene_lists_is_their_union(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path)

        with pytest.raises(GeneCategoryNotDeclaredError):
            registry.gene_list("no_such_category", f"{_WORM}+{_FOOD}")

        union = list(dict.fromkeys([*_declared(_WORM), *_declared(_FOOD)]))
        assert [answer.category for answer in registry.gene_lists()] == union

        # The name of a merge is the +-join of what went in, but it cannot say which
        # component each half came from — and that is exactly what attribution needs, so
        # who contributed comes from the record and not from splitting the name.
        renamed = self._registry(tmp_path / "renamed", name="merged")
        shared = next(category for category in _declared(_WORM) if category in _declared(_FOOD))

        answer = renamed.gene_list(shared, "merged")

        assert [source.component for source in answer.sources] == [
            _WORM_COMPONENT,
            _FOOD_COMPONENT,
        ]


class TestAddressedByAssembly:
    """``gene_list`` and ``gene_lists`` addressed by assembly name rather than opened.

    A registry for the length of the call, exactly as ``annotation_status`` is, so there
    is no second code path to keep in step.
    """

    def test_it_answers_for_an_assembly_named_rather_than_opened_or_raises_if_nothing_is(
        self, tmp_path: Path
    ) -> None:
        source = tmp_path / "ann.gtf"
        source.write_text(_GTF)
        _register_by_path(tmp_path, source, _CURATED, assembly=_CURATED_ASSEMBLY)
        category = _declared(_CURATED)[0]

        answer = gene_list(_CURATED_ASSEMBLY, category, cache_dir=tmp_path)

        assert answer.annotation == _CURATED
        assert [entry.category for entry in gene_lists(_CURATED_ASSEMBLY, cache_dir=tmp_path)] == (
            list(_declared(_CURATED))
        )

        with pytest.raises(AnnotationNotRegisteredError):
            gene_list(_CURATED_ASSEMBLY, "rRNA", annotation=_CURATED, cache_dir=tmp_path / "empty")


def _state_row(**overrides: object) -> AnnotationStatusRow:
    """An offered-but-not-registered status row, with any field overridden by keyword."""
    fields: dict[str, object] = {
        "name": "gencode_v50",
        "offered": True,
        "registered": False,
        "broken": False,
        "default": True,
        "provider": "GENCODE",
        "version": "v50",
        "url": "https://example.org/gencode_v50.gtf.gz",
        "sha256": None,
        "path": None,
        "problem": None,
        "repair": None,
    }
    fields.update(overrides)
    return AnnotationStatusRow(**fields)  # type: ignore[arg-type]


def _state(default: str | None, *rows: AnnotationStatusRow) -> AnnotationStatus:
    """A status report over ``rows``, naming ``default`` as the default annotation."""
    return AnnotationStatus(
        assembly="hg38",
        directory=Path("/data/genome/hg38"),
        default_annotation=default,
        annotations=rows,
    )


class TestTheStateARowIsIn:
    """``AnnotationStatusRow.state`` — the four-way state, in the words a surface prints.

    It was the CLI's to derive, which meant the precedence between ``broken`` and
    ``registered`` was written twice: once as the row's invariant, once as an ``if``
    ordering in a renderer. It is the row's, and the renderer chooses column widths.
    """

    def test_the_four_states_and_brokens_precedence_over_both(self) -> None:
        assert _state_row(registered=True).state == "registered"
        assert _state_row().state == "offered, not registered"
        assert _state_row(offered=False, registered=True, default=False).state == (
            "registered, not offered"
        )
        # Broken beats registered rather than reading as neither: a directory nothing
        # vouches for is not registered, so reporting the absence of a registration would
        # be true and useless. What needs acting on is that it is broken.
        assert _state_row(broken=True).state == "broken"
        assert _state_row(offered=False, broken=True).state == "broken"  # and beats unlisted


class TestTheDefaultAnnotationLine:
    """``AnnotationStatus.default_summary`` — the closing line, and what to do about it.

    Four answers, and the two that name a command take it off an interface: the broken
    one from the row's own ``repair``, the absent one from
    :func:`annotation_register_command`. Neither is concatenated here or in a surface.
    """

    def test_the_four_summary_lines(self) -> None:
        assert _state(None).default_summary == "default: (none)"
        assert (
            _state("gencode_v50", _state_row(registered=True)).default_summary
            == "default: gencode_v50"
        )

        absent = _state("gencode_v50", _state_row())
        assert absent.default_summary == (
            "default: gencode_v50 — not registered here; register it with "
            "`genome annotation register hg38 gencode_v50`"
        )
        # The command it names is the one the package spells once, not a copy of it.
        assert annotation_register_command("hg38", "gencode_v50") in absent.default_summary

        repair = "genome annotation register hg38 gencode_v50 --force"
        broken = _state("gencode_v50", _state_row(broken=True, repair=repair))
        assert broken.default_summary == (
            f"default: gencode_v50 — broken here; repair it with `{repair}`"
        )

        # Named by the table, and no row here is about it: there is no row to read
        # `registered` or `broken` off, and the answer is the one that registers it.
        fresh = _state("gencode_v50")
        assert fresh.default_row is None
        assert "not registered here" in fresh.default_summary


class TestTheJsonKeysAndTheirOrder:
    """``as_json`` — what the status report and its rows serialize as, key for key.

    ``--json`` is what a script parses, so a key renamed, dropped or reordered is a break
    whether or not anything in this suite notices. These assert the whole list rather than
    a key inside it, which is the only form that fails on an addition.
    """

    def test_a_status_row_and_the_status_around_it_serialize_with_no_derived_keys(self) -> None:
        assert list(_state_row().as_json()) == [
            "name",
            "offered",
            "registered",
            "broken",
            "default",
            "provider",
            "version",
            "url",
            "sha256",
            "path",
            "problem",
            "repair",
        ]
        # Writing `state` out would be a second spelling of the precedence for a parser to
        # disagree with, and the three fields it comes from are all here already.
        assert "state" not in _state_row(broken=True).as_json()

        status = _state(
            "gencode_v50", _state_row(), _state_row(name="mine", offered=False, registered=True)
        )
        assert list(status.as_json()) == [
            "assembly",
            "directory",
            "default_annotation",
            "annotations",
        ]
        assert [row["name"] for row in status.as_json()["annotations"]] == ["gencode_v50", "mine"]
        assert status.as_json()["directory"] == "/data/genome/hg38"

    def test_the_whole_report_survives_json_dumps_unchanged(self) -> None:
        # The last step the CLI takes, taken here: nothing in a report is a type json
        # cannot write, which is what makes the `--json` path total.
        status = _state("gencode_v50", _state_row(broken=True, repair="x", problem="y"))

        assert json.loads(json.dumps(status.as_json())) == status.as_json()
