// Copyright 2018 Distributed Bio
// Licensed under the GPL version 3
// This file may not be copied, modified, or distributed
// except according to those terms.

extern crate bio;
extern crate clap;
extern crate sanguine;
use bio::io::fasta::{Reader as FaReader, Record as FaRec, Writer as FaWriter};
use bio::io::fastq::{Reader as FqReader, Record as FqRec, Writer as FqWriter};
use clap::{App, Arg, ArgGroup};
use sanguine::seqio::*;
use sanguine::*;
use sanguine::{align, score_seq, AlnError};
use std::error::Error;
use std::fs::File;
use std::io::{stdout, Error as IOError, Write};

/// FIXME: handle more than Sanger quality offset
const LINKER_LEN: usize = 21;

#[derive(Debug)]
enum FailMode {
    BestQual,
    PolyN,
    Both,
}

fn main() -> Result<(), Box<Error>> {
    let matches = App::new("pair_merge")
        .version("0.1")
        .author("ortho.the.d.stice@gmail.com")
        .about("merge forward and reverse reads with or without quality scores")
        .arg(
            Arg::with_name("fwd")
                .short("f")
                .long("fwd")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("rev")
                .short("r")
                .long("rev")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("out")
                .short("o")
                .long("out")
                .required(true)
                .takes_value(true)
                .help("write joined reads to <out>; set to hyphen for STDOUT"),
        )
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .long("verbose")
                .required(false)
                .takes_value(false)
                .help("verbose mode -- print alignment to STDOUT"),
        )
        .arg(
            Arg::with_name("suppress")
                .short("s")
                .long("suppress-n")
                .required(false)
                .takes_value(false)
                .help("suppress low quality bases by converting to N"),
        )
        .arg(
            Arg::with_name("poly-n")
                .long("poly-n")
                .required(false)
                .takes_value(false)
                .help(
                    "how failing alignments are handled: by default, the read with the \
                     higher quality score is used, but this argument includes both reads\
                     joined with a poly-N linker",
                ),
        )
        .arg(
            Arg::with_name("both")
                .long("both")
                .required(false)
                .takes_value(false)
                .help(
                    "how failing alignments are handled: by default, the read with the \
                     higher quality score is used, but this argument includes both reads",
                ),
        )
        .group(
            ArgGroup::with_name("fail_handler")
                .args(&["poly-n", "both"])
                .required(false)
                .multiple(false), /*
                                  .help("how failing alignments are handled: by default, the read with the\
                                         higher quality score is used, but the user can choose to join the\
                                         reads with a poly-N or include both in the output")*/
        )
        .get_matches();

    let fwd_fname = matches.value_of("fwd").unwrap();
    let rev_fname = matches.value_of("rev").unwrap();
    let out_fname = matches.value_of("out").unwrap();
    //let out_id = matches.value_of("id").or(Some("merged")).unwrap();
    let suppress_n = matches.is_present("suppress");
    let verbose = matches.occurrences_of("verbose") > 0;
    let fail_mode = if matches.is_present("poly-n") {
        FailMode::PolyN
    } else if matches.is_present("both") {
        FailMode::Both
    } else {
        FailMode::BestQual
    };

    println!("fail_mode: {:?}", &fail_mode);

    let fwd_fh = File::open(fwd_fname)?;
    let rev_fh = File::open(rev_fname)?;

    let has_qual =
        fwd_fname.to_lowercase().ends_with("fq") || fwd_fname.to_lowercase().ends_with("fastq");

    let freader = open(
        fwd_fh,
        if has_qual {
            DataFmt::FASTQ
        } else {
            DataFmt::FASTA
        },
    );
    let rreader = open(
        rev_fh,
        if has_qual {
            DataFmt::FASTQ
        } else {
            DataFmt::FASTA
        },
    );

    let fwd_reads = freader.collect::<Result<Vec<Box<SeqRecord>>, IOError>>()?;
    let rev_reads = rreader.collect::<Result<Vec<Box<SeqRecord>>, IOError>>()?;

    if fwd_reads.len() != rev_reads.len() {
        return Err(From::from("equal number of fwd and rev reads required"));
    }

    let mut merged_seqs = vec![];
    let mut merged_q = vec![];

    for (_fwd, _rev) in fwd_reads.iter().zip(rev_reads.iter()) {
        let mut fwd = _fwd.seq().to_owned();
        let rev = _rev.seq().to_owned();

        let (fwd_q, rev_q) = if has_qual {
            let (a, b) = (_fwd.pos_annots(), _rev.pos_annots());
            match (a.as_slice(), b.as_slice()) {
                ([PosAnnot::PhredQual(f)], [PosAnnot::PhredQual(r)]) => {
                    (Some(f.to_owned()), Some(r.to_owned()))
                }
                _ => unreachable!(),
            }
        } else {
            (None, None)
        };

        let apply_suppression = |merged: Vec<u8>, qual: &Option<Vec<u8>>| -> Vec<u8> {
            if suppress_n && has_qual {
                if let Some(ref nq) = qual {
                    merged
                        .iter()
                        .zip(nq)
                        .map(|(ref base, ref q)| qual_n(**base, **q, 5))
                        .collect()
                } else {
                    unreachable!()
                }
            } else {
                merged
            }
        };

        match align(
            fwd.as_ref(),
            fwd_q.as_ref().map(|x| x.as_slice()),
            rev.as_ref(),
            rev_q.as_ref().map(|x| x.as_slice()),
            verbose,
        ) {
            Ok((merged, qual)) => {
                // apply filter to convert low-quality bases to N
                merged_seqs.push(apply_suppression(merged, &qual));
                merged_q.push(qual);
            }
            Err(e) => match e {
                AlnError::NoAlignment => {
                    match fail_mode {
                        FailMode::BestQual => {
                            if score_seq(fwd.as_slice(), &fwd_q)
                                >= score_seq(rev.as_slice(), &rev_q)
                            {
                                merged_seqs.push(apply_suppression(fwd, &fwd_q));
                                merged_q.push(fwd_q);
                            } else {
                                merged_seqs.push(apply_suppression(rev, &rev_q));
                                merged_q.push(rev_q);
                            }
                        }
                        FailMode::PolyN => {
                            // combine reads with N-linker
                            fwd.extend((0..LINKER_LEN).map(|_| b'N'));
                            fwd.extend(rev.iter());

                            merged_seqs.push(fwd);
                            merged_q.push(match (fwd_q, rev_q) {
                                (None, None) => None,
                                (Some(fq), Some(rq)) => {
                                    let mut qual = fq.to_vec();
                                    qual.extend((0..LINKER_LEN).map(|_| 1));
                                    qual.extend(rq.iter());

                                    Some(qual)
                                }
                                _ => unreachable!(),
                            });
                        }
                        FailMode::Both => {
                            merged_seqs.push(fwd);
                            merged_seqs.push(rev);
                            merged_q.push(fwd_q);
                            merged_q.push(rev_q);
                        }
                    }
                }
                _ => return Err(From::from("problem aligning")),
            },
        };
    }

    let out_fh: Box<Write> = if out_fname == String::from("-") {
        Box::new(stdout())
    } else {
        Box::new(File::create(out_fname)?)
    };
    let out_id = "foo";

    if has_qual {
        let mut writer = FqWriter::new(out_fh);
        for (merged, qual) in merged_seqs.iter().zip(merged_q.iter()) {
            let q_conv: Vec<u8> = qual
                .as_ref()
                .unwrap()
                .iter()
                .map(|q| *q + PHRED_OFFSET)
                .collect();
            writer.write(out_id, None, merged.as_ref(), q_conv.as_ref())?;
        }
    } else {
        let mut writer = FaWriter::new(out_fh);

        for merged in merged_seqs.iter() {
            writer.write(out_id, None, merged.as_ref())?;
        }
    };

    Ok(())
}
