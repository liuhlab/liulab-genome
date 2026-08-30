"""What registering an annotation takes, shared by the four modules that need it.

Registering a GTF is the setup half of nearly every test here — whether the subject is
placement, the registry's four states or a **Gene id stem** — so the two ways in live
beside the four test modules rather than being spelled again in each.

Databases are built **for real** from the committed ``tiny.gtf`` fixtures: the build is
most of what registering an annotation is, so substituting gffutils would test nothing.
That needs only gffutils (pure Python + SQLite), never the native bioinformatics
binaries, so nothing here is gated on a tool skip.

Nothing touches the network either: the shared ``fake_fetch`` fixture replaces the
package's one fetch step with a copy out of ``tests/data``, and the annotation table is
injected as an in-memory :class:`AnnotationMetadata` record rather than faked as a TSV —
the shipped table is tested in test_metadata beside it.
"""

from __future__ import annotations

from pathlib import Path

from genome.annotation import AnnotationRegistry, GtfAnnotation
from genome.annotation.metadata import AnnotationMetadata
from genome.assembly.registration import AssemblyDir

# A minimal but valid GTF: one gene with a transcript and an exon. Standard
# gene/transcript features are declared, so the default no-inference path applies.
_GTF = (
    "\n".join(
        [
            'chrI\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
            'chrI\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
            'chrI\ttest\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
        ]
    )
    + "\n"
)


# A bare exon-level GTF: exon lines and nothing else, which is what gene/transcript
# inference exists for. Built with inference off it yields a database of exons alone.
_BARE_GTF = (
    "\n".join(
        [
            'chrI\ttest\texon\t1\t50\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
            'chrI\ttest\texon\t60\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
        ]
    )
    + "\n"
)


#: sha256 of the committed ``tiny.gtf`` — the *unpacked* bytes ``tiny.gtf.gz`` yields.
_TINY_GTF_SHA256 = "255f43bd9abef76424d1c2d89a40cccc1a36215409bbc8f32dcead49ca3baf5e"


#: A URL that is nothing like any provider's, so using it can only come from a row.
_PINNED_URL = "https://mirror.example.invalid/annotations/tiny.gtf.gz"


#: The name the fixture annotation is registered under throughout.
_NAME = "ensgene_v101"


#: An annotation whose curated gene list ships today, and the assembly the table files it
#: under. Named rather than derived: registering it under *another* assembly is what the
#: guard below is about, so the pairing has to be written down somewhere the test controls.
_CURATED = "gencode_v50"


_CURATED_ASSEMBLY = "hg38"


#: The two contributors a **Merged annotation** of worm and its food is made of, and the
#: components whose sequences their features sit on.
_WORM, _WORM_COMPONENT = "wormbase_ws298", "ce11"


_FOOD, _FOOD_COMPONENT = "refseq_rs_2025_06_26", "ecHT115"


_CHIMERA = f"{_WORM_COMPONENT}_{_FOOD_COMPONENT}"


def _row(
    *,
    name: str = _NAME,
    url: str = _PINNED_URL,
    sha256: str | None = None,
    default: bool = True,
) -> AnnotationMetadata:
    """An in-memory annotation row for the ``tiny`` assembly."""
    return AnnotationMetadata(
        assembly="tiny",
        name=name,
        provider="UCSC",
        version="ensGene.v101",
        url=url,
        sha256=sha256,
        default=default,
    )


def _register_by_name(
    assembly_dir: Path,
    assembly: str,
    name: str,
    *,
    force: bool = False,
    progressbar: bool = True,
    metadata: AnnotationMetadata | None = None,
    check_chromosomes: bool = True,
) -> GtfAnnotation:
    """Register ``name`` through a registry bound to ``assembly_dir``.

    A registry addressed by directory rather than opened, which is what these tests want:
    the assembly's ``chrom.sizes`` is found under the directory exactly as an opened
    assembly's is.
    """
    return AnnotationRegistry(AssemblyDir(assembly=assembly, path=assembly_dir)).register(
        name,
        force=force,
        progressbar=progressbar,
        metadata=metadata,
        check_chromosomes=check_chromosomes,
    )


def _register_by_path(
    assembly_dir: Path,
    gtf: str | Path,
    name: str,
    *,
    assembly: str = "tiny",
    chrom_sizes: str | Path | None = None,
    force: bool = False,
    check_chromosomes: bool = True,
    disable_infer_genes: bool = True,
    disable_infer_transcripts: bool = True,
) -> GtfAnnotation:
    """Register the GTF at ``gtf`` under ``assembly_dir`` through a registry.

    The setup half of most of this module: a registry addressed by directory rather than
    opened, answering with the paths, which is what a test asserting on files wants. Left
    to itself it checks against ``<assembly_dir>/<assembly>.chrom.sizes`` — absent in most
    of these directories, so nothing is checked — and ``chrom_sizes`` names another file
    where a test has written one somewhere else.
    """
    return AnnotationRegistry(
        AssemblyDir(assembly=assembly, path=assembly_dir), chrom_sizes=chrom_sizes
    ).register_path(
        gtf,
        name,
        force=force,
        check_chromosomes=check_chromosomes,
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )


def _write_chrom_sizes(assembly_dir: Path, *names: str, assembly: str = "tiny") -> Path:
    """Write ``<assembly>.chrom.sizes`` naming ``names``, where an assembly's own sits."""
    assembly_dir.mkdir(parents=True, exist_ok=True)
    path = assembly_dir / f"{assembly}.chrom.sizes"
    path.write_text("".join(f"{name}\t10000\n" for name in names))
    return path
