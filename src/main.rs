use bio::alignment::pairwise::*;

fn main(){

    if std::env::args().len() != 3 {
        eprintln!("Usage: {} <fasta file> <columns>", std::env::args().next().unwrap());
        std::process::exit(1);
    }

    // Read file from argv
    let infile = std::env::args().nth(1).unwrap();
    let n_columns = std::env::args().nth(2).unwrap().parse::<usize>().unwrap();
    let reader = jseqio::reader::DynamicFastXReader::from_file(&infile).unwrap();
    let db = reader.into_db().unwrap();

    let x = db.get(0).unwrap().seq;
    let y = db.get(1).unwrap().seq;

    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    // gap open score: -5, gap extension score: -1
    let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
    let alignment = aligner.global(x, y);
    eprintln!("{}", alignment.pretty(x, y, n_columns));
}

