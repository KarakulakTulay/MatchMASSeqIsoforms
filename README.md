# MatchMASSeqIsoforms

Iso-Seq pipeline gives different PBids for isoforms even though the isoforms are the same. This script takes two GFF files which are output of MAS-Seq Iso-Seq workflow, compares exon structures of isoforms and match them based on defined nucleotide length.

#### Example Usage:

`$ Rscript MatchMASSeqIsoforms.R sample1/scisoseq_transcripts.sorted.filtered_lite.gff sample2/scisoseq_transcripts.sorted.filtered_lite.gff sample1/scisoseq.seurat_info/isoforms_seurat/genes.tsv sample2/scisoseq.seurat_info/isoforms_seurat/genes.tsv IsoformMatch.csv max5diff max3diff &> IsoformMatch1.log &
`

#### Inputs: 
- scisoseq_transcripts.sorted.filtered_lite.gff: output of Iso-Seq workflow
- genes.tsv: output of Iso-Seq workflow 
- max5diff: max nucleotide length to match transcript at the 5' end (for exp. 50)
- max3diff: max nucleotide length to match transcript at the 3' end (for exp. 100)


