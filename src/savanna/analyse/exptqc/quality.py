import warnings
import pandas as pd
import numpy as np


class ExperimentQualityControl:
    """
    Perform quality control evaluation of an experiment,
    using amplicon coverage data as input

    TODO:
     - Improve docstrings
     - Add some writing methods

    """

    params = {
        "min_cov": 100,
        "min_frac_passing": 0.8,
        "max_per_contamination": 5,
        "min_per_expt_passing": 50,
    }

    def __init__(self, bedcov_df: pd.DataFrame, metadata: pd.DataFrame) -> None:

        # Store amplicon coverage data
        self.bedcov_df = bedcov_df.query("barcode != 'unclassified'")

        # Processes metadata
        self.metadata = metadata
        self._identify_negative_controls()
        self._identify_positive_controls()
        self.metadata.insert(
            4, "is_sample", self.metadata["is_negative"] & self.metadata["is_positive"]
        )
        self.n_barcode = self.metadata.shape[0]
        self.n_samples = self.n_barcode - self.n_negative - self.n_positive

        # Merge
        n = self.bedcov_df.shape[0]
        self.merged_df = pd.merge(
            left=self.bedcov_df, right=self.metadata, on="barcode"
        )
        assert n == self.merged_df.shape[0], "Lost samples during merging."

        # Compute quality
        self.sample_df = self._create_sample_summary()
        self.amplicon_df = self._create_amplicon_summary()
        self.expt_dict = self._create_experiment_summary()

    def _add_indicator_column(self, column_name: str, indicators: str):
        if column_name in self.metadata.columns:
            return

        values = [
            any(indicator in sample_id for indicator in indicators)
            for sample_id in self.metadata["sample_id"]
        ]

        self.metadata.insert(3, column_name, values)

    def _identify_negative_controls(self):
        """
        Identify negative controls within an experiment

        """
        self._add_indicator_column(
            column_name="is_negative", indicators=["NTC", "Water"]
        )
        self.n_negative = self.metadata["is_negative"].sum()

    def _identify_positive_controls(self):
        """
        Identify positive controls within an experiment

        TODO: Note that this is very P.f. specific...

        """
        self._add_indicator_column(
            column_name="is_positive", indicators=["3d7", "3D7", "Dd2", "HB3", "IPC"]
        )
        self.n_positive = self.metadata["is_positive"].sum()

    def _create_sample_summary(self):
        """
        Create a per-sample quality control table

        """

        sample_df = (
            self.merged_df.groupby("barcode")
            .agg(
                sample_id=pd.NamedAgg("sample_id", np.unique),
                is_positive=pd.NamedAgg("is_positive", all),
                is_negative=pd.NamedAgg("is_negative", all),
                n_amplicons=pd.NamedAgg("name", len),
                n_amplicons_pass_cov=pd.NamedAgg(
                    "mean_cov", lambda x: sum(x >= self.params["min_cov"])
                ),
                amplicon_mean_cov=pd.NamedAgg("mean_cov", np.mean),
                amplicon_med_cov=pd.NamedAgg("mean_cov", np.median),
            )
            .reset_index()
        )

        # Per sample
        if self.n_negative == 0:
            self.negative_max_cov = 0
            warnings.warn("No negative controls found. Assuming no contamination.")
        else:
            self.negative_max_cov = sample_df.query("is_negative")[
                "amplicon_mean_cov"
            ].max()

        # Compute an estimate of percentage contamination
        sample_df.insert(
            sample_df.shape[1],
            "amplicon_per_contamination",
            100 * self.negative_max_cov / sample_df["amplicon_mean_cov"],
        )

        # Make final assessment of pass / fail
        sample_df.insert(
            sample_df.shape[1],
            "sample_pass_cov",
            sample_df["n_amplicons_pass_cov"] / sample_df["n_amplicons"]
            >= self.params["min_frac_passing"],
        )
        sample_df.insert(
            sample_df.shape[1],
            "sample_pass_contamination",
            sample_df["amplicon_per_contamination"]
            <= self.params["max_per_contamination"],
        )
        sample_df.insert(
            sample_df.shape[1],
            "sample_pass",
            sample_df["sample_pass_contamination"] & sample_df["sample_pass_cov"],
        )

        return sample_df

    def _create_amplicon_summary(self):
        """"""
        amplicon_df = (
            self.merged_df.groupby("name")
            .agg(
                n_samples=pd.NamedAgg("is_negative", lambda x: len(x) - sum(x)),
                n_samples_pass_cov=pd.NamedAgg(
                    "mean_cov", lambda x: sum(x >= self.params["min_cov"])
                ),
                sample_mean_cov=pd.NamedAgg("mean_cov", np.mean),
                sample_med_cov=pd.NamedAgg("mean_cov", np.median),
            )
            .reset_index()
            .rename({"name": "amplicon"}, axis=1)
        )
        amplicon_df.insert(
            amplicon_df.shape[1],
            "amplicon_cov_pass",
            amplicon_df["n_samples_pass_cov"] / amplicon_df["n_samples"]
            >= self.params["min_frac_passing"],
        )

        return amplicon_df

    def _create_experiment_summary(self):
        """
        Create a final dictionary summarising experiment results

        """

        # TODO: For now, limiting to field samples and positive controls
        # - In future may want to split this
        ndf = self.sample_df.query("not is_negative")
        N = ndf.shape[0]

        per_samples_pass = 100 * ndf["sample_pass"].sum() / N

        return {
            "n_barcodes": int(self.metadata.shape[0]),
            "n_negative_cntrls": int(self.n_negative),
            "n_positive_cntrls": int(self.n_positive),
            "n_samples": int(self.n_samples),
            "n_samples_pass_cov": int(ndf["sample_pass_cov"].sum()),
            "n_samples_pass_contamination": int(ndf["sample_pass_contamination"].sum()),
            "n_samples_pass": int(ndf["sample_pass"].sum()),
            "per_samples_pass": float(per_samples_pass),
            "per_contamination_mean": float(ndf["amplicon_per_contamination"].mean()),
            "expt_pass": bool(per_samples_pass >= self.params["min_per_expt_passing"]),
        }
