# rsomics-bam-to-bed

Convert BAM alignments to BED6 format.

## Usage

```
rsomics-bam-to-bed -i alignments.bam
rsomics-bam-to-bed -i alignments.bam --split   # split spliced reads at N in CIGAR
rsomics-bam-to-bed -i alignments.bam --ed      # use NM tag as score
rsomics-bam-to-bed -i alignments.bam --cigar   # add 7th CIGAR column
```

Output columns: `chrom  start(0-based)  end  name  score  strand`

## Origin

Port of `bedtools bamtobed` based on:
- Quinlan AR & Hall IM (2010). BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* 26(6): 841–842. DOI: 10.1093/bioinformatics/btq033
- The BEDTools source (MIT license) and black-box behavior testing.

License: MIT OR Apache-2.0  
Upstream credit: BEDTools <https://github.com/arq5x/bedtools2> (MIT)
