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

This file re-exports nothing on purpose. Every caller inside and outside the package holds
a *module* — ``from genome.store import fetch``, ``from genome.store.completion import
CompletionRecord`` — and the fetch step's one patch point depends on that: a name
re-exported here is a second reference that ``monkeypatch.setattr`` on the module would
never reach, which is exactly the bug the suite's offline guard exists to prevent.
"""
