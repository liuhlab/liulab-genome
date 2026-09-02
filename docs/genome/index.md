# Working with a genome

Four pages, in the order a project usually needs them: get the reference onto this machine,
read sequence out of it, attach an annotation, then build the index an aligner wants. Every
one of them starts from `Genome("<assembly>")`.

| Page | What it answers |
| --- | --- |
| [Assembly](assembly.md) | How a reference is prepared, and what files that leaves on disk. |
| [Sequences and regions](sequences.md) | How to name a locus, and what a `DNA` object does for you. |
| [Annotations](annotations.md) | How a GTF is registered, and how to ask it which genes it carries. |
| [Aligner indexes](aligner.md) | How a STAR or chromap index is built, and how to find it again. |

Transcription factors, motifs, gene identifiers and orthologs are keyed on a species rather
than on an assembly. They are under [Topics](../topics/index.md).
