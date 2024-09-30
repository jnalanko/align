use bio::alignment::pairwise::*;
use clap::{Command, Arg, ArgAction};

fn main(){

    let cli = Command::new("align")
        .arg_required_else_help(true)
        .arg(Arg::new("fasta")
            .help("Input fasta file containing the two sequences to align.")
            .required(true))
        .arg(Arg::new("columns")
            .long("columns")
            .short('c')
            .help("Number of columns to print.")
            .default_value("80")
        )
        .arg(Arg::new("length-only")
            .help("Print only the alignment length in each query. Only makes sense with local alignment.")
            .long("length-only")
            .conflicts_with("global")
            .conflicts_with("columns")
            .action(ArgAction::SetTrue))
        .arg(Arg::new("global")
            .short('g')
            .long("global")
            .help("Use global alignment (default is local).")
            .action(ArgAction::SetTrue));

    let matches = cli.get_matches();
    let infile = matches.get_one::<String>("fasta").unwrap();
    let n_columns = matches.get_one::<String>("columns").unwrap().parse::<usize>().unwrap();
    let global = matches.get_flag("global");
    let length_only = matches.get_flag("length-only");

    let reader = jseqio::reader::DynamicFastXReader::from_file(&infile).unwrap();
    let db = reader.into_db().unwrap();

    assert!(db.sequence_count() % 2 == 0);

    for pair_idx in 0..(db.sequence_count() / 2) {
        let x = db.get(2*pair_idx).unwrap().seq;
        let y = db.get(2*pair_idx + 1).unwrap().seq;

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

        if length_only {
            let better_len = std::cmp::max(alignment.x_aln_len(), alignment_rc.x_aln_len());
            println!("{}", better_len);
        }
        else {
            #[allow(clippy::collapsible_else_if)] // Clearer this way
            if alignment_rc.score > alignment.score {
                println!("Aligned reverse complement to forward");
                println!("{}", alignment_rc.pretty(x, y, n_columns));
            }
            else {
                println!("Aligned forward to forward");
                println!("{}", alignment.pretty(x, y, n_columns));
            }
        }
    }

}

