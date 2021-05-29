"""Console script for pytxi."""
import argparse
import sys

import pytxi

def main():
    """Console script for pytxi."""
    parser = argparse.ArgumentParser(
        description="Transcript-level to gene-level quantification"
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"pytxi: v{pytxi.__version__}",
    )
    parser.add_argument(
        "input", help="input file(s) in salmon or kallisto format", nargs="+"

    )
    parser.add_argument(
        "outdir", help="name of output directory"
    )
    args = parser.parse_args()


    txi = pytxi.TxImport()
    txi.import_files(args.input)
    txi.abundance.to_csv(f"{outdir}/abundance.tsv", sep="\t")
    txi.counts.to_csv(f"{outdir}/counts.tsv", sep="\t")
    txi.length.to_csv(f"{outdir}/length.tsv", sep="\t")

if __name__ == "__main__":
    sys.exit(main())
