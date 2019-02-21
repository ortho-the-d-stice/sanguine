// Copyright 2018
// Licensed under the GPL version 3
// This file may not be copied, modified, or distributed
// except according to those terms.
extern crate bio;

pub mod seqio;

use bio::alignment::pairwise::MIN_SCORE;
use bio::alignment::pairwise::*;
use bio::alignment::{Alignment, AlignmentOperation as Op};
use bio::alphabets::dna::revcomp;
use std::collections::VecDeque;
use std::fmt;
use std::str;

/// quality below which we use 'N'
const MIN_QUAL: u8 = 10;
pub const PHRED_OFFSET: u8 = 33;

// FIXME: use quick_error or something similar
#[derive(Debug, Clone)]
pub enum AlnError {
    NoAlignment,
    BadParams,
}

struct Seq<'a>(&'a [u8]);
impl<'a> fmt::Display for Seq<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", str::from_utf8(self.0).unwrap())
    }
}

/// given two measurements with uncertainties, return the "true" value that minimizes chi-squared
fn find_true(m1: f32, u1: f32, m2: f32, u2: f32) -> f32 {
    ((m1 / u1.powi(2)) + (m2 / u2.powi(2))) / (u1.powi(-2) + u2.powi(-2))
}

fn phred_to_p(q: f32) -> f32 {
    1.0 - 10_f32.powf(-q / 10.0)
}

fn p_to_phred(p: f32) -> f32 {
    -10.0 * (1.0 - p).log10()
}

/// switch low-quality bases to N
pub fn qual_n(base: u8, qual: u8, min_qual: u8) -> u8 {
    if qual >= min_qual {
        base
    } else {
        b'N'
    }
}

/// determine base identity and phred quality score consistent with two aligned bases.
/// currently, the sum of quality scores from a window around the base is also passed but unused.
pub fn merge_base_with_qual(
    b1: u8,
    q1: u8,
    _sum1: usize,
    b2: u8,
    q2: u8,
    _sum2: usize,
) -> (u8, u8) {
    let p1 = phred_to_p(q1 as f32);
    let p2 = phred_to_p(q2 as f32);

    if b1 == b2 {
        // bases match
        (
            b1,
            p_to_phred(find_true(p1, (1.0 - p1) / 3.0, p2, (1.0 - p2) / 3.0)).round() as u8,
        )
    } else if q1 >= q2 {
        // bases differ
        (
            b1,
            p_to_phred(find_true(p1, (1.0 - p1) / 3.0, 1.0 - p2, (1.0 - p2) / 3.0)).round() as u8,
        )
    } else {
        (
            b2,
            p_to_phred(find_true(p2, (1.0 - p2) / 3.0, 1.0 - p1, (1.0 - p1) / 3.0)).round() as u8,
        )
    }
}

/// pick base without explicit quality score.
/// if the bases differ, and neither is N, N-counts within a window are used as a tie-breaker.
/// if the bases are tied by this measure, N is returned.
fn merge_base(a: u8, a_window_n: u8, b: u8, b_window_n: u8) -> u8 {
    if a == b {
        a
    } else if a == b'N' {
        b
    } else if b == b'N' {
        a
    } else if a_window_n > b_window_n {
        b
    } else if b_window_n > a_window_n {
        a
    } else {
        b'N'
    }
}

fn merge_reads_with_qual(
    ops: &Vec<Op>,
    a: &[u8],
    aq: &[u8],
    b: &[u8],
    bq: &[u8],
) -> (Vec<u8>, Option<Vec<u8>>) {
    let mut sv: Vec<u8> = vec![];
    let mut qv: Vec<u8> = vec![];

    // circular buffer containing quality scores for @a
    // note that a HIGHER sum indicates higher quality
    let mut i_n: VecDeque<u8> = VecDeque::with_capacity(9);
    for _ in 0..4 {
        i_n.push_front(0);
    }
    for q in &aq[..5] {
        i_n.push_front(*q);
    }
    // circular buffer containing quality scores for @b
    let mut j_n: VecDeque<u8> = VecDeque::with_capacity(9);
    for _ in 0..4 {
        j_n.push_front(0);
    }
    for q in &bq[..5] {
        j_n.push_front(*q);
    }

    let mut i = 0;
    let mut j = 0;
    for op in ops.iter() {
        let i_n_ct = i_n.iter().map(|q| *q as usize).sum::<usize>();
        let j_n_ct = i_n.iter().map(|q| *q as usize).sum::<usize>();

        //println!("op={:?}, a[i]={}, i={}, i_n_ct={}, b[j]={}, j={}, j_n_ct={}", &op, a[i] as char, i, i_n_ct, b[j] as char, j, j_n_ct);

        match op {
            Op::Match | Op::Subst => {
                let (nb, nq) = merge_base_with_qual(a[i], aq[i], i_n_ct, b[j], bq[j], j_n_ct);
                sv.push(nb);
                qv.push(nq);

                i_n.pop_back();
                if i + 5 < a.len() {
                    i_n.push_front(aq[i]);
                }
                j_n.pop_back();
                if j + 5 < b.len() {
                    j_n.push_front(bq[j]);
                }
                i += 1;
                j += 1;
            }
            Op::Ins => {
                if i_n_ct <= j_n_ct {
                    sv.push(a[i]);
                    qv.push(aq[i]);
                }
                i_n.pop_back();
                i_n.push_front(aq[i]);
                i += 1;
            }
            Op::Del => {
                if j_n_ct <= i_n_ct {
                    sv.push(b[j]);
                    qv.push(bq[j]);
                }
                j_n.pop_back();
                i_n.push_front(bq[j]);
                j += 1;
            }
            _ => (),
        }
    }
    (sv, Some(qv))
}

/// Merge reads, substituting N where bases disagree.  Indels are handled based on a running tally of N-bases:
/// the read with fewer adjacent N's takes precident.
fn merge_reads(ops: &Vec<Op>, a: &[u8], b: &[u8]) -> Vec<u8> {
    let mut sv: Vec<u8> = vec![];

    let mut i = 0;
    // circular buffer containing quality proxy - 1 for N, 0 otherwise
    // note that a LOWER sum indicates higher quality
    let mut i_n: VecDeque<u8> = VecDeque::with_capacity(9);
    for _ in 0..4 {
        i_n.push_front(0);
    }
    for b in &a[..5] {
        i_n.push_front(if *b == b'N' { 1 } else { 0 })
    }
    let mut j = 0;
    // circular buffer containing quality proxy - 1 for N, 0 otherwise
    let mut j_n: VecDeque<u8> = VecDeque::with_capacity(9);
    for _ in 0..4 {
        j_n.push_front(0);
    }
    for b in &b[..5] {
        j_n.push_front(if *b == b'N' { 1 } else { 0 })
    }

    for op in ops.iter() {
        let i_n_ct = i_n.iter().cloned().sum::<u8>();
        let j_n_ct = j_n.iter().cloned().sum::<u8>();

        match op {
            Op::Match | Op::Subst => {
                sv.push(merge_base(a[i], i_n_ct, b[j], j_n_ct));

                i_n.pop_back();
                if i + 5 < a.len() {
                    i_n.push_front(if a[i + 5] == b'N' { 1 } else { 0 });
                }
                j_n.pop_back();
                if j + 5 < b.len() {
                    j_n.push_front(if b[i + 5] == b'N' { 1 } else { 0 });
                }
                i += 1;
                j += 1;
            }
            Op::Del => {
                if i_n_ct <= j_n_ct {
                    sv.push(a[i]);
                }
                i_n.pop_back();
                i_n.push_front(if i + 5 < a.len() && a[i + 5] == b'N' {
                    1
                } else {
                    0
                });
                i += 1;
            }
            Op::Ins => {
                if j_n_ct <= i_n_ct {
                    sv.push(b[j]);
                }
                j_n.pop_back();
                i_n.push_front(if j + 5 < b.len() && b[j + 5] == b'N' {
                    1
                } else {
                    0
                });
                j += 1;
            }
            _ => (),
        }
    }
    sv
}

/// Find the "best" alignment w/r/t to two issues:
/// 1. the "tgt" sequence might be a reverse-compliment
/// 2. the "semiglobal" mode can produce weird outputs like this:
///     AAAAAAAAAAAAAAAA-------------------AAAAAAAAAAAA
///            BBBBBBBBBBBBBBBBBBBBBBBBBBBB
///    which manifests as a double-clipping even, ie, starting with xclip & yclip,
///    or else ending w/ xclip & yclip.
///    we don't want this.
///
/// Returns tuple of: Alignment, Target (possibly reverse-compliment), and Target Quality
/// tgt and tgt_q are consumed and returned, possibly unaltered
fn best_aln(
    refr: &[u8],
    refr_q: Option<&[u8]>,
    tgt: Vec<u8>,
    tgt_q: Option<Vec<u8>>,
) -> (Alignment, Vec<u8>, Option<Vec<u8>>) {
    fn num_clips(ops: &[Op]) -> usize {
        ops.iter()
            .map(|op| match *op {
                Op::Xclip(_) | Op::Yclip(_) => 1,
                _ => 0,
            })
            .sum()
    }

    let xclip = -10; // * (refr.len() / 5) as i32;
    let yclip = -10; // * (_tgt.len() / 5) as i32;
    fn match_fn(a: u8, b: u8) -> i32 {
        if a == b'N' || b == b'N' {
            0
        } else {
            if a == b {
                1_i32
            } else {
                -5_i32
            }
        }
    };
    let mut clip_aligner = Aligner::with_scoring(Scoring {
        gap_open: -5,
        gap_extend: -1,
        match_fn: match_fn,
        match_scores: Some((1, -3)),
        xclip_prefix: xclip,
        xclip_suffix: xclip,
        yclip_prefix: yclip,
        yclip_suffix: yclip,
    });
    let mut reg_aligner = Aligner::with_scoring(Scoring {
        gap_open: -5,
        gap_extend: -1,
        match_fn: match_fn,
        match_scores: Some((1, -3)),
        xclip_prefix: MIN_SCORE,
        xclip_suffix: MIN_SCORE,
        yclip_prefix: MIN_SCORE,
        yclip_suffix: MIN_SCORE,
    });

    let (refr_n, tgt_n) = match (&refr_q, &tgt_q) {
        (&Some(ref rq), &Some(ref tq)) => (
            refr.iter()
                .zip(rq.iter())
                .map(|(ref base, ref q)| qual_n(**base, **q, 10))
                .collect(),
            tgt.iter()
                .zip(tq.iter())
                .map(|(ref base, ref q)| qual_n(**base, **q, 10))
                .collect(),
        ),
        _ => (refr.to_vec(), tgt.to_vec()),
    };

    let aln1 = {
        // default to "clip" aligner
        let a = clip_aligner.custom(refr_n.as_ref(), tgt_n.as_ref());
        if num_clips(&a.operations[..2]) > 1
            || num_clips(&a.operations[&a.operations.len() - 3..]) > 1
        {
            // only use "regular" aligner if clip-aligner is weird
            reg_aligner.custom(refr_n.as_ref(), tgt_n.as_ref())
        } else {
            a
        }
    };

    let rc: Vec<u8> = revcomp(tgt_n.as_slice());
    let aln2 = {
        let a = clip_aligner.custom(refr_n.as_slice(), rc.as_slice());
        if num_clips(&a.operations[..2]) > 1
            || num_clips(&a.operations[a.operations.len() - 3..]) > 1
        {
            // only use "regular" aligner if clip-aligner is weird
            reg_aligner.custom(refr_n.as_slice(), rc.as_ref())
        } else {
            a
        }
    };

    if aln1.score >= aln2.score {
        (aln1, tgt, tgt_q)
    } else {
        match tgt_q {
            Some(mut tgt_q) => {
                tgt_q.reverse();
                (aln2, rc, Some(tgt_q))
            }
            None => (aln2, rc, None),
        }
    }
}

/// Use Rust Bio's "clipped" aligner to align two overlapping reads and combine phred scores
/// where appropriate.  Reads can be in same direction, or reverse complements, as both will be
/// tried and the combination with the greater alignment score used.
///
/// Indels: if one read has an insertion relative to the other, the extra base will be included if
/// quality scores aren't provided.  If quality scores ARE present, the base will be included if
/// it is present in the read with the greater local quality (window size 5).
///
/// NOTE: quality scores should be adjusted, ie, 0 should mean 0
///
pub fn align(
    _tgt: &[u8],
    _tgt_q: Option<&[u8]>,
    refr: &[u8],
    refr_q: Option<&[u8]>,
    verbose: bool,
) -> Result<(Vec<u8>, Option<Vec<u8>>), AlnError> {
    // can't combine a sequence with quality and one without
    match (refr_q, _tgt_q) {
        (None, Some(_)) => return Err(AlnError::BadParams),
        (Some(_), None) => return Err(AlnError::BadParams),
        _ => (),
    }

    // try given direction and reverse-complement
    let (aln, tgt, mut tgt_q) = best_aln(
        refr,
        refr_q,
        _tgt.to_vec(),
        _tgt_q.and_then(|ref q| Some(q.to_vec())),
    );

    if verbose {
        println!("## score: {}", aln.score);
        println!("{}", aln.pretty(&refr, &tgt));
    }

    // we don't know which read is "first"
    let (pre, pre_q, post, post_q) = if aln.xstart > aln.ystart {
        (
            &refr[..aln.xstart],
            refr_q.and_then(|ref q| Some(&q[..aln.xstart])),
            &tgt[aln.yend..],
            match &tgt_q {
                &None => None,
                &Some(ref q) => Some(&q[aln.yend..])
            }
            //tgt_q.and_then(|ref q| Some(&q[aln.yend..])),
        )
    } else {
        (
            &tgt[..aln.ystart],
            match &tgt_q {
                &None => None,
                &Some(ref q) => Some((&q[..aln.ystart])),
            },
            //tgt_q.and_then(|ref q| Some(&q[..aln.ystart])),
            &refr[aln.xend..],
            refr_q.and_then(|ref q| Some(&q[aln.xend..])),
        )
    };

    if aln.score < 0 || aln.xlen <= 5 || aln.ylen < 5 {
        // effectively no overlap, often b/c one or both reads are low-quality
        return Err(AlnError::NoAlignment);
    }

    let (middle, mid_q): (Vec<u8>, Option<Vec<u8>>) = match (&refr_q, &tgt_q) {
        (&Some(ref rq), &Some(ref tq)) => merge_reads_with_qual(
            &aln.operations,
            &refr[aln.xstart..aln.xend],
            &rq[aln.xstart..aln.xend],
            &tgt[aln.ystart..aln.yend],
            &tq[aln.ystart..aln.yend],
        ),
        (&None, &None) => (
            merge_reads(
                &aln.operations,
                &tgt[aln.ystart..aln.yend],
                &refr[aln.xstart..aln.xend],
            ),
            None,
        ),
        _ => unreachable!(),
    };

    let mut merged = Vec::new();
    merged.extend(pre.iter());
    merged.extend(middle.iter());
    merged.extend(post.iter());

    let new_qual = refr_q.and_then(|_| {
        // FIXME: can we initialize this buffer initially and build as we go?
        let mut nq = vec![];
        nq.extend(pre_q.expect("").iter());
        nq.extend(mid_q.expect("").iter());
        nq.extend(post_q.expect("").iter());
        Some(nq)
    });

    Ok((merged, new_qual))
}

/// sequence quality score based on length and quality, or, barring q-scores,
/// the number of N-bases
pub fn score_seq(seq: &[u8], seq_q: &Option<Vec<u8>>) -> f32 {
    match seq_q {
        &Some(ref q) => q.iter().map(|q| phred_to_p(*q as f32)).sum(),
        &None => seq.iter().map(|b| if *b == b'N' { 0.0 } else { 1.0 }).sum(),
    }
}

#[cfg(test)]
mod test_aln {
    use super::*;

    #[test]
    fn single_base() {
        assert_eq!(merge_base_with_qual(b'A', 20, 0, b'A', 20, 0), (b'A', 20));

        // NOTE: b'A' is chosen arbitrarily b/c it's the first argument, and
        // quality score 3 corresponds w/ a 49.9% chance
        assert_eq!(merge_base_with_qual(b'A', 20, 0, b'T', 20, 0), (b'A', 3));
        // NOTE: 20 -> 0.99, 15 -> 0.968, 10 -> 0.9
        // this doesn't match my intuition, but then, the whole thing is hand-waving.
        // what's more, any base under 20 is often considered dubious.
        assert_eq!(merge_base_with_qual(b'T', 20, 0, b'A', 15, 0), (b'T', 10));
    }

    #[test]
    fn simple_merge() {
        // ATGGCCTCTTTATGAGCGATGGTTGCTCGGATAGTAGATA
        //      CTCTTTATGAGCGATGGTTGCTCGGATAGTAGATAGGGGCCC

        assert_eq!(
            align(
                b"ATGGCCTCTTTATGAGCGATGGTTGCTCGGATAGTAGATA",
                None,
                b"CTCTTTATGAGCGATGGTTGCTCGGATAGTAGATAGGGGCCC",
                None,
                false
            )
            .unwrap(),
            (
                b"ATGGCCTCTTTATGAGCGATGGTTGCTCGGATAGTAGATAGGGGCCC".to_vec(),
                None
            )
        );

        //     ATGCTTTATGAGCGATGGTTGCTCGGATAGTAGGCCTATA
        // GGCAATGCTTTATGAGCGATGGTTGCTCGGATAGTAGGCCT
        assert_eq!(
            align(
                b"ATGCTTTATGAGCGATGGTTGCTCGGATAGTAGGCCTATA",
                None,
                b"GGCAATGCTTTATGAGCGATGGTTGCTCGGATAGTAGGCCT",
                None,
                false
            )
            .unwrap(),
            (
                b"GGCAATGCTTTATGAGCGATGGTTGCTCGGATAGTAGGCCTATA".to_vec(),
                None
            )
        );
    }

    #[test]
    fn single_mismatch() {
        // TTTATGGCCTATAGGCTTTATGAGCGATGGTTGCTCGGATAGTAG
        //            x
        //       GCCTAAAGGCTTTATGAGCGATGGTTGCTCGGATAGTAGGGCCC

        let refr = b"TTTATGGCCTATAGGCTTTATGAGCGATGGTTGCTCGGATAGTAG".to_vec();
        let tgt = b"GCCTAAAGGCTTTATGAGCGATGGTTGCTCGGATAGTAGGGCCC".to_vec();

        // no quality scores
        assert_eq!(
            align(refr.as_slice(), None, tgt.as_slice(), None, true).unwrap(),
            (
                b"TTTATGGCCTANAGGCTTTATGAGCGATGGTTGCTCGGATAGTAGGGCCC".to_vec(),
                None
            )
        );

        // same but with scores - all bases same [high] quality
        let rq = vec![20; refr.len()];
        let mut tq = vec![20; tgt.len()];
        assert_eq!(
            align(refr.as_slice(), Some(&rq), tgt.as_slice(), Some(&tq), true)
                .unwrap()
                .0,
            b"TTTATGGCCTANAGGCTTTATGAGCGATGGTTGCTCGGATAGTAGGGCCC".to_vec()
        );

        // assign low quality score to inconsistent base - high quality base should take precidence
        tq[5] = 10;
        let (merged_seq, _merged_qual) =
            align(refr.as_slice(), Some(&rq), tgt.as_slice(), Some(&tq), false).unwrap();
        assert_eq!(
            merged_seq,
            b"TTTATGGCCTATAGGCTTTATGAGCGATGGTTGCTCGGATAGTAGGGCCC".to_vec()
        );
        let merged_qual = _merged_qual.unwrap();
        assert_eq!(merged_qual[10], 20);
        assert!(merged_qual[11] < 20);
        assert_eq!(merged_qual[12], 20);
    }

    #[test]
    fn simple_indel() {
        //                                                 ||||||
        let a = b"ATGTATCTGCTTTTTTGCACTCAAATTTAACCCATATGCTCGCATGGCGGTCT";
        let b = b"AAATTTAACCCATATGCTCGCCATGGCGGTCTACGTGGCACCTAATGGCAGGCCCCGAATCTAACGTTAGTCTGCCCGA";
        //                           ||-||||
        let merged = align(a, None, b, None, false).unwrap();
        assert_eq!(merged.0,
                   b"ATGTATCTGCTTTTTTGCACTCAAATTTAACCCATATGCTCGCCATGGCGGTCTACGTGGCACCTAATGGCAGGCCCCGAATCTAACGTTAGTCTGCCCGA"
                   .to_vec());
    }

    #[test]
    fn indel_with_n() {
        //                                                 ||||||
        let a = b"ATGTATCTGCTTTTTTGCACTCAAATTTAACCCATATGCTCGCATGGCGGTCT";
        let b = b"AAATTTAACCCATATGCTCGCCATNGCGGTCTACGTGGCACCTAATGGCAGGCCCCGAATCTAACGTTAGTCTGCCCGA";
        //                           ||-|| |
        let merged = align(a, None, b, None, false).unwrap();
        assert_eq!(merged.0,
                   b"ATGTATCTGCTTTTTTGCACTCAAATTTAACCCATATGCTCGCATGGCGGTCTACGTGGCACCTAATGGCAGGCCCCGAATCTAACGTTAGTCTGCCCGA"
                   .to_vec());
    }

    #[test]
    fn indel_with_qual() {
        let a: &[u8] = b"ATGTATCTGCTTTTTTGCACTCAAATTTAACCCATATGCTCGCATGGCGG";
        let b: &[u8] = b"ATGTATCTGCTTTTTTGCACTCAAATTTAACCCATATGCTGCATGGCGG";
        //               ATGTATCTGCTTTTTTGCACTCAAATTTAACCCATATGCTCGCATGGCGG
        //                                               *********
        let aq: &[u8] = &[40_u8; 50];
        let bq: &mut [u8] = &mut [40_u8; 49];
        for i in 40..49 {
            bq[i] = 10;
        }

        let mut reg_aligner = Aligner::with_scoring(Scoring {
            gap_open: -5,
            gap_extend: -1,
            match_fn: |a, b| {
                if a == b'N' || b == b'N' {
                    0
                } else {
                    if a == b {
                        1_i32
                    } else {
                        -5_i32
                    }
                }
            },
            match_scores: Some((1, -3)),
            xclip_prefix: MIN_SCORE,
            xclip_suffix: MIN_SCORE,
            yclip_prefix: MIN_SCORE,
            yclip_suffix: MIN_SCORE,
        });

        let aln = reg_aligner.custom(a, b);
        println!("{}", aln.pretty(a, b));

        // FIXME: what's up w/ the ordering here?
        let (merged, qual) = merge_reads_with_qual(&aln.operations, a, aq, b, bq);
        assert_eq!(a, merged.as_slice());
    }

    #[test]
    fn test_score_seq() {
        assert!(score_seq(b"ATGCCGT", &None) > score_seq(b"ATGNNNNNNNNNNN", &None));
        assert!(
            score_seq(b"ATGCCGT", &Some(vec![20; 7])) > score_seq(b"ATGCCGT", &Some(vec![10; 7]))
        );
    }
}
