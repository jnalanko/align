use bio::alignment::pairwise::*;
use clap::{Command, Arg, ArgAction};

fn main(){

    let cli = Command::new("align")
        .arg_required_else_help(true)
        .arg(Arg::new("fasta")
            .help("Input fasta file")
            .required(true)
            .index(1))
        .arg(Arg::new("columns")
            .help("Number of columns to print")
            .required(true)
            .index(2))
        .arg(Arg::new("global")
            .short('g')
            .long("global")
            .help("Use global alignment (default is local)")
            .action(ArgAction::SetTrue));

    let matches = cli.get_matches();
    let infile = matches.get_one::<String>("fasta").unwrap();
    let n_columns = matches.get_one::<String>("columns").unwrap().parse::<usize>().unwrap();
    let global = matches.get_flag("global");

    let reader = jseqio::reader::DynamicFastXReader::from_file(&infile).unwrap();
    let db = reader.into_db().unwrap();

    let x = db.get(0).unwrap().seq;
    let y = db.get(1).unwrap().seq;

    let x_rc_vec = x.iter().rev().map(|c| bio::alphabets::dna::complement(*c)).collect::<Vec<_>>();
    let x_rc = x_rc_vec.as_slice();

    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    // gap open score: -5, gap extension score: -1
    let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);

    let alignment = match global {
        true => aligner.global(x, y),
        false => aligner.local(x, y),
    };
    let alignment_rc = match global {
        true => aligner.global(x_rc, y),
        false => aligner.local(x_rc, y),
    };
    
    if alignment_rc.score > alignment.score {
        eprintln!("Aligned reverse complement to forward");
        eprintln!("{}", alignment_rc.pretty(x, y, n_columns));
    }
    else {
        eprintln!("Aligned forward to forward");
        eprintln!("{}", alignment.pretty(x, y, n_columns));
    }
}

