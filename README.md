# sanguine - a fwd + rev read merger for Sanger sequences with support for FASTQ output

A number of tools exist to merge pair-end NGS reads ([PEAR](https://sco.h-its.org/exelixis/web/software/pear/doc.html), [PANDASEQ](https://github.com/neufeld/pandaseq)) quickly and accurately, but fewer tools are optimized for Sanger sequences.  `sanguine` addresses this simpler use-case, and includes support for FASTQ output.

## Example
```
sanguine --fwd fwd.fastq --rev rev.fastq --out merged.fastq
```

## Algorithm

`sanguine` aligns the given reads using the [Rust-Bio library](https://github.com/rust-bio/rust-bio), trying both the supplied orientation and the reverse-complement.  Discrepencies between sequences are handled differently depending whether quality scores are available.

### FASTA input

If the supplied sequences lack quality scores, a sliding-window count of N-bases is used as a proxy.  In the example that follows, the presence of N's in the forward sequence suggests lower sequence quality, causing the reverse base to be favored where the sequences disagree (eg, position 5).

```
        012345678901

   fwd: AATANTGNCGAA
          ||  | ||||
   rev:   TACCGACGAATTAA

merged: AATACCGACGAATTAA
```

### FASTQ input

Given a discrepency between reads with explicit quality scores, `sanguine` will choose the base with the higher quality.  `sanguine` also combines input quality scores to produce a quality score for the merged sequence by treating the quality scores as measurements with uncertainty.

1. Q-score translated to probability: `P = 10<sup>-Q/10</sup>`
2. If corresponding bases disagree, the base with the higher score is chosen, and `1 - P` is used for the lower-quality base.
3. Uncertainty is computed as `(1-P)/3`

Base<sub>1</sub> | Qual<sub>1</sub> | Base<sub>2</sub> | Qual<sub>2</sub> | Merged base | Merged Q-score
---------------- | ---------------- | ---------------- | ---------------- | ----------- | --------------
A | 20 | A | 20 | A | 20
A | 20 | G | 5 | A | 20
A | 20 | G | 20 | G | 3

