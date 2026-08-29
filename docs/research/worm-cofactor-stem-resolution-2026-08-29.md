# Worm cofactor stem resolution against the annotation the lab registers on ce11

**Date:** 2026-08-29. **Method:** the WormBase annotation registered on GPU71FM
(`ce11/wormbase_ws298`), read straight out of the annotation database its registration built —
`SELECT id FROM features WHERE featuretype = 'gene'` — and joined against the **Gene id stem**s of
the shipped *C. elegans* **Cofactor table**. No package code was involved, so this measures the
annotation and the table rather than the resolver that reads them. It is the worm counterpart of the
human and mouse measurement taken the same day against the GENCODE registrations.

## The result

| | |
|---|---|
| gene features in `wormbase_ws298` | 46,926 |
| stems naming more than one gene id | **0** |
| cofactor stems | 317 |
| resolved | **317 / 317** |
| unresolved | 0 |

**Every worm cofactor resolves, and nothing collides.** This is the cleanest of the three species
measured: human loses three stems to Lambert's UniProt entry names and one to an id GENCODE 50 has
dropped, and mouse loses six, while worm loses none at all.

## Why it comes out this clean

WormBase gene ids are unversioned as published — `WBGene00000064`, never `WBGene00000064.1` — so a
stem and a gene id are the same string here and stem resolution is an identity rather than a
version-stripping join. That is the case the resolver already handles by answering with the stem
itself, so the path this exercises is the *absence* of an Ensembl-shaped assumption rather than the
version machinery. AnimalTFDB keys its worm cofactors on WormBase ids directly, so no crosswalk sits
between the publisher's file and the annotation — unlike the human table, whose EpiFactors half
reaches Ensembl only through HGNC.

## What it settles

The parent specification listed this as an open item: *whether worm `WBGene…` stems resolve against
the worm annotation the lab registers on ce11*, unverified, and explicitly not a reason to withhold
the table. It is now measured, and the answer is that they all do.

**This number cannot drift, which is unusual.** WS298 is not merely the latest WormBase release but
its **final** one: the project announced that WS298 "will be the final major WormBase release" and
that the site would afterwards be "anchored to the WS298 database release and entered into an
archival maintenance mode", with ongoing *C. elegans* curation moving to the Alliance of Genome
Resources —
<https://community.alliancegenome.org/t/announcing-the-final-release-of-wormbase/8461>, 2025-07-24.
The announcement is a separate post and is **not** in `letter.WS298`, which reads like an ordinary
release letter; the wormbase.org blog original returns 403 to automated clients, so cite the Alliance
copy. So 317 of 317 is a permanent fact about this annotation rather than a reading that a later
release could move, and pinning worm to WS298 is permanently equivalent to pinning it to the latest.
Anything built on worm after this should look to the Alliance, not to a WS299 that will not come.
