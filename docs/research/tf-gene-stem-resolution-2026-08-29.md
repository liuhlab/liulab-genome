# Gene id stem resolution against the annotations the lab actually registers

**Date:** 2026-08-29. **Method:** the three GENCODE annotations registered on GPU71FM
(`hg38/gencode_v50`, `hg19/gencode_v50lift37`, `mm39/gencode_vM39`), read straight out of the
annotation database each registration built — `SELECT id FROM features WHERE featuretype = 'gene'` —
and joined against the **Gene id stem**s of the two shipped censuses. No package code was involved,
so this measures the annotations and the censuses rather than the resolver that reads them.

## How many gene ids one stem names

| Annotation | gene features | stems naming more than one gene id |
| --- | --- | --- |
| `gencode_v50` | 78,733 | 0 |
| `gencode_v50lift37` | 80,315 | **9** |
| `gencode_vM39` | 78,289 | 0 |

The nine on `gencode_v50lift37`, in full:

```text
ENSG00000223274   ENSG00000223274.6_3    ENSG00000223274.1_PAR_Y
ENSG00000251823   ENSG00000251823.2_2    ENSG00000251823.1
ENSG00000263835   ENSG00000263835.1      ENSG00000263835.1_PAR_Y
ENSG00000263980   ENSG00000263980.1      ENSG00000263980.1_PAR_Y
ENSG00000264510   ENSG00000264510.1      ENSG00000264510.1_PAR_Y
ENSG00000264819   ENSG00000264819.1      ENSG00000264819.1_PAR_Y
ENSG00000265350   ENSG00000265350.1      ENSG00000265350.1_PAR_Y
ENSG00000265658   ENSG00000265658.1      ENSG00000265658.1_PAR_Y
ENSG00000266731   ENSG00000266731.1      ENSG00000266731.1_PAR_Y
```

Eight are pseudoautosomal-Y; the ninth, `ENSG00000251823`, is a lift-mapping duplicate. This is the
measurement behind the rule that resolution answers with **every** gene id a stem names and never
picks one — on this annotation, picking the first would be wrong nine times and silently.

**`gencode_v50` carries no `_PAR_Y` gene at all** (zero of 78,733). The pseudoautosomal duplicate is
a property of the lift, not of GENCODE 50, so a fixture asserting the two-ids case against an
`hg38`-shaped id is exercising the mechanism rather than reproducing a state that annotation is in.

## How much of a census a real annotation carries

| Annotation | stems resolved | gene ids | unresolved |
| --- | --- | --- | --- |
| `hg38/gencode_v50` | 1,636 / 1,639 | 1,636 | 3 |
| `hg19/gencode_v50lift37` | 1,633 / 1,639 | 1,633 | 6 |
| `mm39/gencode_vM39` | 1,605 / 1,611 | 1,605 | 6 |

Counted over the assessed-positive genes of each species' census.

The three unresolved on `gencode_v50` are `DUX1_HUMAN`, `DUX3_HUMAN` and `ENSG00000215271`. The
first two are the UniProt entry names Lambert et al. records in its Ensembl id column for a handful
of genes — they are not Ensembl gene ids and resolve against nothing, by construction. The third is
an Ensembl id GENCODE 50 no longer carries. All three ride back on the answer rather than being
dropped, which is the whole reason unresolved stems are reported.

## The finding worth keeping

**No TF gene collides in any annotation the lab currently registers.** Every one of the nine
colliding stems on `gencode_v50lift37` is a gene no census assessed, and the human census's
assessed-positive stems each name exactly one gene id on all three annotations — the `stems naming 2
ids` count is zero in every row of the second table.

So the "never pick one" rule is not idle — `gencode_v50lift37` really does collide nine times — but
nothing in the TF path exercises it against real data today. It is a guard against a state the
annotations are in and the censuses have not yet met, which is the right way round: the alternative
is discovering it when a census grows into one of those nine.
