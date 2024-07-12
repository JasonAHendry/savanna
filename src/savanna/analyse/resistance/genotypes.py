from typing import Dict

import pandas as pd


class SampleGenotypes:
    """
    Encapsulate genotype calls for a sample in an easy to
    query format

    genotypes = {
        "crt": {"K76T": 2, ...}
        "dhfr": {"N51I": 2, "S108N": 0, ...}
        ...
    }

    TODO:
    - Could also be useful for creating summaries of the genotype information

    """

    UNPHASED_GT_TO_INT = {"0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2, "./.": -1}

    def __init__(self, name: str, genotypes: Dict[str, Dict[str, int]]) -> None:
        """Initialise field, sanity check, compute some summaries"""
        # Sanity check
        if not isinstance(name, str):
            raise ValueError("name must be a string.")
        if not isinstance(genotypes, Dict):
            raise ValueError("genotypes must be a dictioary.")

        # Main attributes
        self.name = name
        self.genotypes = genotypes

        # summaries
        self.genes = list(genotypes.keys())

    def __str__(self):
        return f"SampleGenotypes(sample_name={self.name}, genes={self.genes})"

    def __repr__(self):
        return self.__str__()

    @classmethod
    def from_sample_dataframe(cls, sample_name: str, sample_df: pd.DataFrame):
        """
        Instantiate from a dataframe
        """

        genotypes = {}
        for gene, gene_df in sample_df.groupby("gene"):
            genotypes[gene] = {
                aa: cls.UNPHASED_GT_TO_INT[gt]
                for aa, gt in zip(gene_df["aa_change"], gene_df["gt"])
            }

        return cls(sample_name, genotypes)
