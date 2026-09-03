"""What belongs to no context and is reached by all of them.

Four modules that know about bytes, directories and hashes, and nothing about assemblies,
annotations, motifs or genes. Each is here because every context needs it and none owns it:

- :mod:`~genome.store.data_dir` reads ``$LIULAB_DATA`` — the **Data dir** the whole
  package files into. Which roots sit under it is each context's own to declare.
- :mod:`~genome.store.fetch` is the package's one fetch step, and imports nothing else
  from this package at all.
- :mod:`~genome.store.completion` is the **Completion marker**: the record a finished
  build writes and the working area it uses until it does.
- :mod:`~genome.store.checksum` digests a file, and refuses one that disagrees with what
  was expected of it.
- :mod:`~genome.store.prepared` is the pipeline the three **Prepared set**s share, built
  out of the four above and out of nothing else.

**Nothing here imports from** :mod:`genome.assembly` **or** :mod:`genome.annotation`, which
is what makes this a place both of them can reach without reaching each other.

**This file re-exports no callable on purpose, and the fetch step is why.** Every caller
inside and outside the package holds a *module* — ``from genome.store import fetch``,
``from genome.store.completion import CompletionRecord`` — and the fetch step's one patch
point depends on that: a callable re-exported here is a second reference that
``monkeypatch.setattr`` on the module would never reach, which is exactly the bug the
suite's offline guard exists to prevent.

**The exception classes below are the one exemption, because nothing patches one.** They
are re-exported so a caller can name in an ``except`` what this package hands them — a
**Assembly dir** that disagrees with its record, a **Prepared set** nothing has
prepared yet — rather than importing from a module the API reference declares free to move. The
exemption is theirs alone: a function or a non-exception class added to ``__all__``
re-opens the hole the paragraph above closes.

Examples
--------
>>> from genome.store import RegistrationMismatchError
>>> issubclass(RegistrationMismatchError, RuntimeError)
True
"""

from genome.store.checksum import ChecksumMismatchError
from genome.store.completion import (
    RegistrationError,
    RegistrationMismatchError,
    UnfinishedRegistrationError,
)
from genome.store.prepared import (
    PreparedChecksumError,
    PreparedDecodeError,
    PreparedSetNotDownloadedError,
)

__all__ = [
    "ChecksumMismatchError",
    "PreparedChecksumError",
    "PreparedDecodeError",
    "PreparedSetNotDownloadedError",
    "RegistrationError",
    "RegistrationMismatchError",
    "UnfinishedRegistrationError",
]
