# 3. An assembly id is a local key; UCSC is the default source, not the namespace

The `assembly` argument is a free-form local key: it names the data directory and is the lookup key
into a curated cross-reference table carrying six identifiers across three naming authorities
(`hg38` / `GRCh38` / `GCF_000001405.40`). `lookup_assembly` returning `None` is legal and every field
is overridable at construction, so an assembly absent from the table is fine. Validation is a
property of the *source*, not of the identifier — `validate_assembly` runs only on the UCSC fetch
path, never on `fetch_genome_from`. The obvious reading, that an assembly *is* a UCSC id, loses
because it would foreclose the mixed sources the `path_or_url` seeding path already supports.
