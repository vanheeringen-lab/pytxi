import sys

import pandas as pd
import numpy as np
from loguru import logger

from genomepy.annotation import query_mygene
from genomepy import Genome

__version__ = "0.1.1"


FMTS = {
    "kallisto": (
        ["length", "eff_length", "est_counts", "tpm"],
        ["tpm", "est_counts", "eff_length"],
    ),
    "salmon": (
        ["Length", "EffectiveLength", "TPM", "NumReads"],
        ["TPM", "NumReads", "EffectiveLength"],
    ),
}


KALLISTO_COLS = ["length", "eff_length", "est_counts", "tpm"]
SALMON_COLS = ["Length", "EffectiveLength", "TPM", "NumReads"]


class TxImport:
    tx2gene = None
    abundance = None
    length = None
    counts = None

    def __init__(self):
        pass

    @staticmethod
    def _parse_species(species):
        if species is None:
            logger.info("Using default tax_id 9606 (human)")
            return 9606

        try:
            # tax_id
            int(species)
            logger.info(f"Using tax_id {species}")
            return int(species)
        except (ValueError, TypeError):
            pass

        try:
            # genomepy.Genome object
            tax_id = species.tax_id
            logger.info(f"Using tax_id {tax_id} from genome {species.name}")
            return tax_id
        except AttributeError:
            pass

        try:
            # genomepy genome name/fasta
            g = Genome(species)
            logger.info(f"Using tax_id {g.tax_id} from genome {species.name}")
            return g.tax_id
        except FileNotFoundError:
            logger.error(
                f"Provided species is not a tax_id and I cannot find a genome with the name {species}"
            )
            logger.error("Don't know what to do now :(")
            sys.exit()

    def set_tx2gene(self, tx2gene=None, transcript_ids=None, species=None):
        if tx2gene:
            logger.info("Using provided tx2gene file")
            result = pd.read_csv(tx2gene, index_col=0)
            result.columns = ["symbol"]
        else:
            logger.info("Mapping transcripts to genes using mygene.info")
            tax_id = self._parse_species(species)
            if not isinstance(transcript_ids, pd.Series):
                transcript_ids = pd.Series(transcript_ids)

            transcripts = transcript_ids.str.replace(
                r"\.[\d_]+$", "", regex=True
            )
            result = query_mygene(transcripts, tax_id, "symbol")
            result = result[["symbol"]]
        self.tx2gene = result

    def import_files(
        self, fnames, sample_names=None, tx2gene=None, species=None
    ):
        """Convert transcriptome-level quantification to gene_level quantification.

        Parameters
        ----------
        fnames : list
            List of file names.
        sample_names : list, optional
            Use these sample names. If not specified, the name of the directories
            is used as sample name.
        tx2gene : str, optional
            Filename of transcript to gene mapping. Should contain the transcript_id
            in the first column and the gene_id in the second column. If not specified,
            genomepy is used to query mygene to automatically determine the mapping.
        species : str, optional
            Species to use for mygene query. Human is set as default. Can be taxonomy
            id (int or string), a genomepy genome name, or a Genome instance.
        """
        if not sample_names:
            sample_names = [fname.split("/")[-2] for fname in fnames]

        dfs = [pd.read_table(fname, index_col=0) for fname in fnames]

        for filetype, (cols, use_cols) in FMTS.items():
            if dfs[0].shape[1] == len(cols) and list(dfs[0].columns) == cols:
                logger.info(f"Detected {filetype} files")
                tpm_col, counts_col, length_col = use_cols
                break
        else:
            logger.error("Unknown filetype")
            logger.error(dfs[0].columns)
            sys.exit()

        tpm = pd.concat([df[tpm_col] for df in dfs], axis=1)
        tpm.columns = sample_names

        self.set_tx2gene(tx2gene, tpm.index, species)

        counts = pd.concat([df[counts_col] for df in dfs], axis=1)
        counts.columns = sample_names

        length = pd.concat([df[length_col] for df in dfs], axis=1)
        length.columns = sample_names

        default_len = length.join(self.tx2gene).groupby("symbol").mean().mean(1)

        abundance = tpm.join(self.tx2gene).groupby("symbol").sum()
        counts = counts.join(self.tx2gene).groupby("symbol").sum()

        length = length.mul(tpm)
        length[length == 0] = np.nan
        length = (
            length.join(self.tx2gene).groupby("symbol").sum().div(abundance)
        )

        f = length.mean(1).isna()
        for col in length.columns:
            length.loc[f, col] = default_len.loc[f]

        mean_length = length.mean(1)
        for col in length.columns:
            length.loc[length[col].isna(), col] = mean_length[
                length[col].isna()
            ]

        self.abundance = abundance
        self.length = length
        self.counts = counts

    def mean(self, fname=None, sample_names=None):
        if not sample_names:
            sample_names = self.counts.columns

        df = self.abundance[sample_names].mean(1).to_frame("tpm")
        df["length"] = self.length[sample_names].mean(1)
        df["counts"] = self.counts[sample_names].mean(1)

        if not fname:
            return df

        df.to_csv(fname, sep="\t", float_format="%0.4f")
