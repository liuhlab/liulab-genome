"""GTF annotation registration — place a GTF and build its gffutils database.

A reference assembly may carry several gene annotations (e.g. GENCODE,
RefSeq, Ensembl). Each is registered under a unique ``name`` and lives in its
own directory beside the assembly's sequence files::

    <LIULAB_DATA>/genome/<assembly>/gtf/<name>/
        <name>.gtf      # the (unzipped) annotation, copied in
        <name>.db       # the gffutils SQLite database built from it

:func:`register_gtf` does the placement + database build; :func:`list_annotations`
discovers what has already been registered. Richer GTF query types build on the
:class:`GtfAnnotation` records these return — that work lives elsewhere; this
module is only the I/O boundary.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

import gffutils

#: Subdirectory under an assembly's data dir holding all its GTF annotations.
_GTF_SUBDIR = "gtf"


@dataclass(frozen=True)
class GtfAnnotation:
    """A registered GTF annotation: its name and the on-disk GTF + database paths."""

    name: str
    gtf: Path
    db: Path


def _annotations_root(assembly_dir: Path) -> Path:
    """Return ``<assembly_dir>/gtf``, the parent of every annotation directory."""
    return assembly_dir / _GTF_SUBDIR


def annotation_dir(assembly_dir: Path, name: str) -> Path:
    """Return the directory holding the annotation registered as ``name``."""
    return _annotations_root(assembly_dir) / name


def _annotation_files(assembly_dir: Path, name: str) -> GtfAnnotation:
    """Resolve the GTF + database paths for ``name`` (without checking existence)."""
    directory = annotation_dir(assembly_dir, name)
    return GtfAnnotation(name=name, gtf=directory / f"{name}.gtf", db=directory / f"{name}.db")


def list_annotations(assembly_dir: Path) -> dict[str, GtfAnnotation]:
    """Return registered annotations (those with a built database), keyed by name."""
    root = _annotations_root(assembly_dir)
    if not root.is_dir():
        return {}
    found: dict[str, GtfAnnotation] = {}
    for directory in sorted(root.iterdir()):
        if not directory.is_dir():
            continue
        annotation = _annotation_files(assembly_dir, directory.name)
        if annotation.db.exists():
            found[annotation.name] = annotation
    return found


def register_gtf(
    assembly_dir: Path,
    gtf: str | Path,
    name: str,
    *,
    force: bool = False,
    disable_infer_genes: bool = True,
    disable_infer_transcripts: bool = True,
) -> GtfAnnotation:
    """Register an unzipped GTF under ``name`` and build its gffutils database.

    Gene/transcript inference is disabled by default — standard annotation GTFs
    (GENCODE, Ensembl, RefSeq) already declare ``gene``/``transcript`` features,
    and inferring them is the classic gffutils slow path. Enable it only for a
    bare exon-level GTF.
    """
    source = Path(gtf)
    if not source.is_file():
        raise FileNotFoundError(f"GTF file not found: {source}")
    if source.suffix == ".gz":
        raise ValueError(
            f"GTF must be unzipped, got {source.name!r}; decompress it first (e.g. `gunzip`)."
        )

    annotation = _annotation_files(assembly_dir, name)
    if not force and (annotation.db.exists() or annotation_dir(assembly_dir, name).exists()):
        raise FileExistsError(
            f"annotation {name!r} is already registered for this assembly at "
            f"{annotation_dir(assembly_dir, name)}; pass force=True to overwrite."
        )

    annotation.db.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() != annotation.gtf.resolve():
        shutil.copy2(source, annotation.gtf)

    gffutils.create_db(
        str(annotation.gtf),
        str(annotation.db),
        force=force,
        keep_order=True,
        merge_strategy="create_unique",
        sort_attribute_values=True,
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )
    return annotation
