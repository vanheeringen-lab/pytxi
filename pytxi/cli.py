"""Console script for pytxi."""
import argparse
import sys
import os

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
    parser.add_argument("outdir", help="name of output directory")
    parser.add_argument(
        "-s",
        "--species",
        help="species tax_id or genomepy genome name (default is human, 9606)",
        default=None,
    )
    parser.add_argument(
        "--tx2gene",
        help="tx2gene file (default is lookup from mygene.info)",
        default=None,
    )
    if len(sys.argv) == 1:
        return parser.print_help()
    args = parser.parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    txi = pytxi.TxImport()
    txi.import_files(args.input, tx2gene=args.tx2gene, species=args.species)
    txi.abundance.to_csv(f"{outdir}/abundance.tsv", sep="\t")
    txi.counts.to_csv(f"{outdir}/counts.tsv", sep="\t")
    txi.length.to_csv(f"{outdir}/length.tsv", sep="\t")


if __name__ == "__main__":
    sys.exit(main())
