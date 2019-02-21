use bio::io::fasta::{Reader as FaReader, Record as FaRecord, Records as FaRecords};
use bio::io::fastq::{Reader as FqReader, Record as FqRecord, Records as FqRecords};
use bio::utils::TextSlice;
use std::io;
use std::io::Read;
use *;

#[derive(Debug, PartialEq)]
pub enum PosAnnot {
    PhredQual(Vec<u8>),
}

#[derive(Debug, PartialEq)]
pub enum DataFmt {
    FASTA,
    FASTQ,
}

pub trait SeqRecord {
    fn id(&self) -> &str;
    fn seq(&self) -> TextSlice;
    fn pos_annots(&self) -> Vec<PosAnnot>;
}

/// NOTE: to test: this should be more efficient if we switch to (1) allocating a single boxed record as
/// part of the struct, and then (2) using reader.read to overwrite the contents each time
struct FaWrapper<R: Read>(FaRecords<R>);
impl<R: Read> Iterator for FaWrapper<R> {
    type Item = Result<Box<dyn SeqRecord>, io::Error>;
    fn next(&mut self) -> Option<Result<Box<dyn SeqRecord>, io::Error>> {
        self.0.next().and_then(|rec| {
            Some(match rec {
                Ok(x) => {
                    let y: Box<dyn SeqRecord> = Box::new(x);
                    Ok(y)
                }
                Err(e) => Err(e),
            })
        })
    }
}

struct FqWrapper<R: Read>(FqRecords<R>);
impl<R: Read> Iterator for FqWrapper<R> {
    type Item = Result<Box<dyn SeqRecord>, io::Error>;
    fn next(&mut self) -> Option<Result<Box<dyn SeqRecord>, io::Error>> {
        self.0.next().and_then(|rec| {
            Some(match rec {
                Ok(x) => {
                    let y: Box<dyn SeqRecord> = Box::new(x);
                    Ok(y)
                }
                Err(e) => Err(e),
            })
        })
    }
}

pub fn open<R>(
    input: R,
    fmt: DataFmt,
) -> Box<Iterator<Item = Result<Box<dyn SeqRecord>, io::Error>>>
where
    R: Read + 'static,
{
    match fmt {
        DataFmt::FASTA => Box::new(FaWrapper(FaReader::new(input).records())),
        DataFmt::FASTQ => Box::new(FqWrapper(FqReader::new(input).records())),
    }
}

impl SeqRecord for FaRecord {
    fn id(&self) -> &str {
        self.id()
    }
    fn seq(&self) -> TextSlice {
        self.seq()
    }
    fn pos_annots(&self) -> Vec<PosAnnot> {
        vec![]
    }
}

impl SeqRecord for FqRecord {
    fn id(&self) -> &str {
        self.id()
    }
    fn seq(&self) -> TextSlice {
        self.seq()
    }
    fn pos_annots(&self) -> Vec<PosAnnot> {
        let q = self.qual().iter().map(|q| q - PHRED_OFFSET).collect();
        vec![PosAnnot::PhredQual(q)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FASTQ_FILE: &'static [u8] = b"@id desc
ACCGTAGGCTGA
+
IIIIIIJJJJJJ
";
    const FASTA_FILE: &'static [u8] = b">id desc
ACCGTAGGCTGA
CCGTAGGCTGAA
CGTAGGCTGAAA
GTAGGCTGAAAA
CCCC
>id2
ATTGTTGTTTTA
ATTGTTGTTTTA
ATTGTTGTTTTA
GGGG
";

    #[test]
    fn test_seqio_fastq() {
        let reader = FqReader::new(FASTQ_FILE);
        for rec in reader.records() {
            let r = rec.unwrap();
            assert_eq!(
                r.pos_annots()[0],
                PosAnnot::PhredQual(r.qual().iter().map(|q| q - PHRED_OFFSET).collect())
            );
        }
    }
    #[test]
    fn test_seqio_fasta() {
        let reader = FaReader::new(FASTA_FILE);
        for rec in reader.records() {
            let r = rec.unwrap();
            assert_eq!(r.pos_annots().len(), 0);
        }
    }

    #[test]
    fn test_seqio_both() {
        let mut v: Vec<Box<SeqRecord>> = vec![];

        let reader = FaReader::new(FASTA_FILE);
        for rec in reader.records() {
            v.push(Box::new(rec.unwrap()));
        }
        let reader = FqReader::new(FASTQ_FILE);
        for rec in reader.records() {
            v.push(Box::new(rec.unwrap()));
        }
    }

    #[test]
    fn test_seqio_open() {
        let refr: &[u8] = b"ACCGTAGGCTGACCGTAGGCTGAACGTAGGCTGAAAGTAGGCTGAAAACCCC";
        let rec = open(FASTA_FILE, DataFmt::FASTA).next().unwrap().unwrap();
        assert_eq!(rec.seq(), refr);
        assert_eq!(rec.id(), String::from("id"));
    }
}
