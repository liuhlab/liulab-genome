"""Tests for genome.aligner — the Aligner abstraction, STAR, and the mixin.

The pure-logic tests stub out the STAR binary (its resolution, version, and the
subprocess call), so they run anywhere. A couple of integration tests build a
real index from a toy FASTA + GTF and are skipped when STAR is not on ``PATH``.
"""

from __future__ import annotations

import json
import shutil
import types
from pathlib import Path
from typing import TYPE_CHECKING, cast

import pandas as pd
import pytest

import genome.aligner.aligner as aligner_mod
from genome.aligner.chromap import Chromap
from genome.aligner.chromap import _kwargs_to_flags as _chromap_kwargs_to_flags
from genome.aligner.mixin import AlignerMixin, _resolve_aligner
from genome.aligner.star import STAR, _kwargs_to_flags
from genome.external import ToolNotFoundError

if TYPE_CHECKING:
    from genome.genome import Genome

_STAR_PRESENT = shutil.which("STAR") is not None
_CHROMAP_PRESENT = shutil.which("chromap") is not None

# A 10 kb single-chromosome genome — large enough for STAR to index, small
# enough to build in a second or two.
_SEQ = "ACGTACGTAGGCATCGATCG" * 500

# One gene, two exons (GTF is 1-based inclusive): exon1 [101, 300], exon2
# [601, 800]. STAR should derive the single intron chr1:301-600 (+).
_TOY_GTF = (
    'chr1\ttoy\texon\t101\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\ttoy\texon\t601\t800\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
)


def _make_genome(tmp_path: Path, gtfs: dict[str, str] | None = None) -> Genome:
    """Return a minimal Genome-like stub backed by a real FASTA + GTF(s) on disk.

    ``gtfs`` maps annotation name -> GTF text; each is written to ``tmp_path`` and
    exposed through a ``get_gtf_path`` matching :meth:`Genome.get_gtf_path`.
    """
    fasta = tmp_path / "tiny.fa"
    fasta.write_text(">chr1\n" + _SEQ + "\n")

    gtf_text = gtfs or {"toy": _TOY_GTF}
    gtf_paths: dict[str, Path] = {}
    for name, content in gtf_text.items():
        path = tmp_path / f"{name}.gtf"
        path.write_text(content)
        gtf_paths[name] = path

    stub = types.SimpleNamespace(
        assembly="tiny",
        files=types.SimpleNamespace(fasta=fasta),
        chrom_sizes=pd.Series({"chr1": len(_SEQ)}),
        get_gtf_path=lambda name: gtf_paths[name],
    )
    return cast("Genome", stub)


@pytest.fixture
def stub_star(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> STAR:
    """A STAR (bound to the ``toy`` annotation) with faked binary + version."""
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    monkeypatch.setattr(aligner_mod, "_resolve", lambda name: f"/fake/{name}")
    monkeypatch.setattr(STAR, "_detect_version", lambda _self: "0.0-test")
    return STAR(_make_genome(tmp_path), gtf="toy")


@pytest.fixture
def stub_chromap(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Chromap:
    """A Chromap with faked binary + version (no annotation — chromap needs none)."""
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    monkeypatch.setattr(aligner_mod, "_resolve", lambda name: f"/fake/{name}")
    monkeypatch.setattr(Chromap, "_detect_version", lambda _self: "0.0-test")
    return Chromap(_make_genome(tmp_path))


@pytest.fixture
def captured_run(monkeypatch: pytest.MonkeyPatch) -> list[list[str]]:
    """Record (and suppress) every ``Aligner._run`` invocation's argument list."""
    calls: list[list[str]] = []
    monkeypatch.setattr(aligner_mod.Aligner, "_run", lambda self, args: calls.append(list(args)))
    return calls


# -- pure logic -------------------------------------------------------------


def test_kwargs_to_flags_scalar_and_list() -> None:
    flags = _kwargs_to_flags({"genomeSAindexNbases": 11, "genomeFastaFiles": ["a.fa", "b.fa"]})
    assert flags == ["--genomeSAindexNbases", "11", "--genomeFastaFiles", "a.fa", "b.fa"]


def test_missing_tool_prints_instructions_and_raises(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    def _missing(name: str) -> str:
        raise ToolNotFoundError("nope")

    monkeypatch.setattr(aligner_mod, "_resolve", _missing)
    with pytest.raises(ToolNotFoundError, match="required to build a star index"):
        STAR(_make_genome(tmp_path), gtf="toy")
    err = capsys.readouterr().err
    assert "STAR is not installed" in err
    assert "bioconda" in err


def test_index_dir_is_per_annotation(stub_star: STAR) -> None:
    # The annotation key names the genomeDir, so different GTFs never collide.
    assert stub_star.index_dir.parts[-4:] == ("genome", "tiny", "index", "star_toy")


def test_distinct_gtf_keys_use_distinct_index_dirs(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    monkeypatch.setattr(aligner_mod, "_resolve", lambda name: f"/fake/{name}")
    monkeypatch.setattr(STAR, "_detect_version", lambda _self: "0.0-test")

    genome = _make_genome(tmp_path, {"a": _TOY_GTF, "b": _TOY_GTF})
    star_a = STAR(genome, gtf="a")
    star_b = STAR(genome, gtf="b")

    assert star_a.index_dir != star_b.index_dir
    assert star_a.index_dir.name == "star_a"
    assert star_b.index_dir.name == "star_b"


def test_index_path_raises_before_build(stub_star: STAR) -> None:
    with pytest.raises(RuntimeError, match="No successful star index"):
        _ = stub_star.index_path


def test_index_writes_metadata_flag_and_returns_dir(
    stub_star: STAR, captured_run: list[list[str]]
) -> None:
    out = stub_star.index(threads=3)

    assert out == stub_star.index_dir == stub_star.index_path
    assert (out / ".success").is_file()
    assert len(captured_run) == 1

    meta = json.loads((out / "star.index.json").read_text())
    assert meta["aligner"] == "star"
    assert meta["version"] == "0.0-test"
    assert meta["assembly"] == "tiny"
    assert meta["parameters"]["threads"] == 3
    assert meta["parameters"]["gtf"] == "toy"
    # Small genome -> a reduced suffix-array size is auto-added.
    assert "genomeSAindexNbases" in meta["parameters"]

    cmd = meta["command"]
    assert cmd[0] == "STAR"
    assert "genomeGenerate" in cmd
    assert "--genomeFastaFiles" in cmd


def test_index_is_cached_and_overwrite_rebuilds(
    stub_star: STAR, captured_run: list[list[str]]
) -> None:
    stub_star.index()
    stub_star.index()  # success flag present -> reused, no rebuild
    assert len(captured_run) == 1

    stub_star.index(overwrite=True)  # forced
    assert len(captured_run) == 2


def test_index_emits_sjdb_flags_from_bound_gtf(
    stub_star: STAR, captured_run: list[list[str]], tmp_path: Path
) -> None:
    stub_star.index(sjdb_overhang=49)

    gtf_path = str((tmp_path / "toy.gtf").resolve())
    args = captured_run[0]
    assert "--sjdbGTFfile" in args
    assert gtf_path in args
    assert "--sjdbOverhang" in args
    assert "49" in args

    meta = json.loads((stub_star.index_dir / "star.index.json").read_text())
    assert meta["parameters"]["gtf"] == "toy"
    assert meta["parameters"]["sjdb_gtf_file"] == gtf_path
    assert meta["parameters"]["sjdb_overhang"] == 49


# -- mixin: build/get index entry points ------------------------------------


@pytest.fixture
def mixin_genome(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Genome:
    """A Genome-like object carrying :class:`AlignerMixin`, with STAR stubbed."""
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    monkeypatch.setattr(aligner_mod, "_resolve", lambda name: f"/fake/{name}")
    monkeypatch.setattr(STAR, "_detect_version", lambda _self: "0.0-test")
    monkeypatch.setattr(Chromap, "_detect_version", lambda _self: "0.0-test")

    stub = _make_genome(tmp_path)

    class _MixinGenome(AlignerMixin):
        pass

    genome = _MixinGenome()
    for attr in ("assembly", "files", "chrom_sizes", "get_gtf_path"):
        setattr(genome, attr, getattr(stub, attr))
    return cast("Genome", genome)


def test_resolve_aligner_is_case_insensitive() -> None:
    assert _resolve_aligner("star") is STAR
    assert _resolve_aligner("STAR") is STAR


def test_resolve_aligner_unknown_raises() -> None:
    with pytest.raises(ValueError, match="Unknown aligner 'bowtie'"):
        _resolve_aligner("bowtie")


def test_get_index_returns_built_index_path(
    mixin_genome: Genome, captured_run: list[list[str]]
) -> None:
    built = mixin_genome.build_star_index("toy")
    assert mixin_genome.get_index("star", gtf="toy") == built


def test_get_star_index_matches_generic_get_index(
    mixin_genome: Genome, captured_run: list[list[str]]
) -> None:
    mixin_genome.build_star_index("toy")
    assert mixin_genome.get_star_index("toy") == mixin_genome.get_index("star", gtf="toy")


def test_get_index_raises_before_build(mixin_genome: Genome) -> None:
    with pytest.raises(RuntimeError, match="No successful star index"):
        mixin_genome.get_index("star", gtf="toy")


def test_get_index_unknown_aligner_raises(mixin_genome: Genome) -> None:
    with pytest.raises(ValueError, match="Unknown aligner 'bowtie'"):
        mixin_genome.get_index("bowtie", gtf="toy")


# -- integration (require a real STAR) --------------------------------------


@pytest.mark.skipif(not _STAR_PRESENT, reason="STAR not on PATH")
def test_real_star_index_builds(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    star = STAR(_make_genome(tmp_path), gtf="toy")

    out = star.index(threads=2)

    assert out == star.index_path
    assert out.name == "star_toy"
    names = {p.name for p in out.iterdir()}
    assert {"SA", "SAindex", "Genome", "genomeParameters.txt"} <= names
    assert (out / ".success").is_file()


@pytest.mark.skipif(not _STAR_PRESENT, reason="STAR not on PATH")
def test_real_star_index_with_gtf(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    star = STAR(_make_genome(tmp_path), gtf="toy")

    out = star.index(sjdb_overhang=49, threads=2)

    # GTF-derived annotation files are produced.
    assert (out / "geneInfo.tab").is_file()
    sjdb = (out / "sjdbList.out.tab").read_text()
    # The single intron between the two exons: chr1:301-600.
    assert "301" in sjdb
    assert "600" in sjdb

    meta = json.loads((out / "star.index.json").read_text())
    assert meta["parameters"]["gtf"] == "toy"
    assert meta["parameters"]["sjdb_gtf_file"] == str((tmp_path / "toy.gtf").resolve())


# ===========================================================================
# Chromap — a second aligner: one index per assembly, no annotation, a file
# (not a directory) as the artifact, and hyphenated long flags.
# ===========================================================================


# -- pure logic -------------------------------------------------------------


def test_chromap_kwargs_to_flags_hyphenates_and_expands_lists() -> None:
    # chromap's long options are hyphenated, so underscores become hyphens.
    flags = _chromap_kwargs_to_flags({"min_frag_length": 30, "read_format": ["r1", "bc"]})
    assert flags == ["--min-frag-length", "30", "--read-format", "r1", "bc"]


def test_missing_chromap_prints_instructions_and_raises(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    def _missing(name: str) -> str:
        raise ToolNotFoundError("nope")

    monkeypatch.setattr(aligner_mod, "_resolve", _missing)
    with pytest.raises(ToolNotFoundError, match="required to build a chromap index"):
        Chromap(_make_genome(tmp_path))
    err = capsys.readouterr().err
    assert "chromap is not installed" in err
    assert "bioconda" in err


def test_chromap_index_dir_is_per_assembly(stub_chromap: Chromap) -> None:
    # No annotation -> no per-GTF suffix; one chromap index serves the assembly.
    assert stub_chromap.index_dir.parts[-3:] == ("tiny", "index", "chromap")


def test_chromap_index_path_raises_before_build(stub_chromap: Chromap) -> None:
    with pytest.raises(RuntimeError, match="No successful chromap index"):
        _ = stub_chromap.index_path


def test_chromap_index_writes_metadata_flag_and_returns_file(
    stub_chromap: Chromap, captured_run: list[list[str]]
) -> None:
    out = stub_chromap.index()

    # The artifact is a single file inside the index dir, not the dir itself.
    assert out == stub_chromap.index_path
    assert out.name == "chromap.index"
    assert out.parent == stub_chromap.index_dir
    assert (stub_chromap.index_dir / ".success").is_file()
    assert len(captured_run) == 1

    meta = json.loads((stub_chromap.index_dir / "chromap.index.json").read_text())
    assert meta["aligner"] == "chromap"
    assert meta["version"] == "0.0-test"
    assert meta["assembly"] == "tiny"

    cmd = meta["command"]
    assert cmd[0] == "chromap"
    assert "--build-index" in cmd
    assert "--ref" in cmd
    assert "--output" in cmd


def test_chromap_default_index_command_is_minimal(
    stub_chromap: Chromap, captured_run: list[list[str]], tmp_path: Path
) -> None:
    # With no tuning kwargs, the command is exactly build-index over ref -> output.
    stub_chromap.index()
    fasta = str(tmp_path / "tiny.fa")
    artifact = str(stub_chromap.index_dir / "chromap.index")
    assert captured_run[0] == ["--build-index", "--ref", fasta, "--output", artifact]


def test_chromap_index_is_cached_and_overwrite_rebuilds(
    stub_chromap: Chromap, captured_run: list[list[str]]
) -> None:
    stub_chromap.index()
    stub_chromap.index()  # success flag present -> reused, no rebuild
    assert len(captured_run) == 1

    stub_chromap.index(overwrite=True)  # forced
    assert len(captured_run) == 2


def test_chromap_index_emits_kmer_window_flags(
    stub_chromap: Chromap, captured_run: list[list[str]]
) -> None:
    stub_chromap.index(kmer=20, window=10, min_frag_length=25)

    args = captured_run[0]
    assert args[args.index("--kmer") + 1] == "20"
    assert args[args.index("--window") + 1] == "10"
    assert args[args.index("--min-frag-length") + 1] == "25"

    meta = json.loads((stub_chromap.index_dir / "chromap.index.json").read_text())
    assert meta["parameters"]["kmer"] == 20
    assert meta["parameters"]["window"] == 10
    assert meta["parameters"]["min_frag_length"] == 25


# -- mixin: build/get index entry points ------------------------------------


def test_resolve_aligner_includes_chromap() -> None:
    assert _resolve_aligner("chromap") is Chromap
    assert _resolve_aligner("CHROMAP") is Chromap


def test_build_and_get_chromap_index(mixin_genome: Genome, captured_run: list[list[str]]) -> None:
    built = mixin_genome.build_chromap_index()
    assert mixin_genome.get_index("chromap") == built
    assert built.name == "chromap.index"


def test_get_chromap_index_matches_generic_get_index(
    mixin_genome: Genome, captured_run: list[list[str]]
) -> None:
    mixin_genome.build_chromap_index()
    assert mixin_genome.get_chromap_index() == mixin_genome.get_index("chromap")


def test_get_chromap_index_raises_before_build(mixin_genome: Genome) -> None:
    with pytest.raises(RuntimeError, match="No successful chromap index"):
        mixin_genome.get_chromap_index()


# -- integration (require a real chromap) -----------------------------------


@pytest.mark.skipif(not _CHROMAP_PRESENT, reason="chromap not on PATH")
def test_real_chromap_index_builds(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    chromap = Chromap(_make_genome(tmp_path))

    out = chromap.index()

    assert out == chromap.index_path
    assert out.name == "chromap.index"
    assert out.is_file()  # chromap's index is a single file, not a directory
    assert out.stat().st_size > 0
    assert (chromap.index_dir / ".success").is_file()

    meta = json.loads((chromap.index_dir / "chromap.index.json").read_text())
    assert meta["aligner"] == "chromap"
    assert meta["version"]  # a real version string was detected
