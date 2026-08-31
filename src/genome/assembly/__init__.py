"""One **Assembly** on disk — which reference it is, where its files live, how a locus becomes bases.

The whole Assembly context, and its own I/O with it: nothing outside this package fetches
a FASTA, derives a companion file or decides where an assembly's directory is.

- :mod:`~genome.assembly.metadata` is the curated table — what the lab supports and how
  each reference is named across databases.
- :mod:`~genome.assembly.source` resolves a name into where its bytes come from, and
  :mod:`~genome.assembly.download` gets them; :mod:`~genome.assembly.registration` is the
  half of registering that has nothing to do with where they came from — the **Assembly
  dir** layout, the staging, the **Completion marker**.
- :mod:`~genome.assembly.fasta` derives the companions with **External tool**s and
  :mod:`~genome.assembly.twobit` reads bases back out of the ``.2bit``.
- :mod:`~genome.assembly.chimera` is the naming rules a **Chimera** obeys and
  :mod:`~genome.assembly.chimera_build` writes one; :mod:`~genome.assembly.components`
  is what such a build recorded, read back.
- :mod:`~genome.assembly.genome` is all of it opened: :class:`~genome.assembly.genome.Genome`.

**What this file re-exports stops short of two modules, and the omission is load-bearing.**
:class:`~genome.assembly.genome.Genome` and the chimera build both reach into
:mod:`genome.annotation`, which reaches back here for the **Assembly dir** — so importing
either of them from this file would make ``import genome.annotation`` run the whole
open-a-genome stack through a package it is halfway through importing. ``Genome`` is
exported from :mod:`genome` itself, which is where a caller holds it anyway, and so is
:class:`~genome.assembly.chimera_build.AmbiguousDefaultAnnotationError` — the one error
the chimera build raises. It is public at the root because of this import edge, not
because it belongs to this context any less than
:class:`~genome.assembly.chimera.ChimeraNamingError` beside it.

Examples
--------
>>> from genome.assembly import assembly_metadata
>>> assembly_metadata("hg38").ncbi_name
'GRCh38'
"""

from genome.assembly.chimera import (
    ChimeraNamingError,
    derive_name,
    split_name,
    split_suffixed,
    suffixed,
)
from genome.assembly.components import (
    COMPONENTS_UNCHANGED,
    COMPONENTS_UNKNOWN,
    ChimeraDetails,
    ComponentDetails,
    components_status,
    read_chimera_details,
)
from genome.assembly.download import (
    EXPECTED_FROM_RECORD,
    EXPECTED_FROM_TABLE,
    RegisteredAssembly,
    UCSCGenomeDownloader,
    VerifiedAssembly,
    assembly_table_row,
    register_assembly,
    verify_assembly,
)
from genome.assembly.fasta import GenomeFiles, prepare_fasta, read_chrom_sizes
from genome.assembly.metadata import (
    AssemblyMetadata,
    assembly_metadata,
    assembly_table,
    format_table_row,
    lookup_assembly,
)
from genome.assembly.registration import (
    ANNOTATIONS_SUBDIR,
    INDEXES_SUBDIR,
    AssemblyDir,
    AssemblyRegistration,
    assembly_data_dir,
    assembly_repair_command,
)
from genome.assembly.twobit import TwoBit

__all__ = [
    "ANNOTATIONS_SUBDIR",
    "COMPONENTS_UNCHANGED",
    "COMPONENTS_UNKNOWN",
    "EXPECTED_FROM_RECORD",
    "EXPECTED_FROM_TABLE",
    "INDEXES_SUBDIR",
    "AssemblyDir",
    "AssemblyMetadata",
    "AssemblyRegistration",
    "ChimeraDetails",
    "ChimeraNamingError",
    "ComponentDetails",
    "GenomeFiles",
    "RegisteredAssembly",
    "TwoBit",
    "UCSCGenomeDownloader",
    "VerifiedAssembly",
    "assembly_data_dir",
    "assembly_metadata",
    "assembly_repair_command",
    "assembly_table",
    "assembly_table_row",
    "components_status",
    "derive_name",
    "format_table_row",
    "lookup_assembly",
    "prepare_fasta",
    "read_chimera_details",
    "read_chrom_sizes",
    "register_assembly",
    "split_name",
    "split_suffixed",
    "suffixed",
    "verify_assembly",
]
