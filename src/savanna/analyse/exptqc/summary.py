from typing import List

import json
import pandas as pd

from savanna.analyse._interfaces import SummaryAnalysis
from savanna.analyse.exptqc.quality import ExperimentQualityControl


class SampleQualityControl(SummaryAnalysis):

    name = "exptqc"

    def _define_inputs(self) -> List[str]:
        """
        Define input files for sample quality control analysis
        """

        self.bedcov_csv = (
            f"{self.expt_dirs.approach_dir}/summary/bedcov/summary.bedcov.csv"
        )

        return [self.bedcov_csv]

    def _define_outputs(self) -> List[str]:

        self.sample_csv = f"{self.output_dir}/summary.sample_qc.csv"
        self.amplicon_csv = f"{self.output_dir}/summary.amplicon_qc.csv"
        self.expt_qc_json = f"{self.output_dir}/summary.experiment_qc.json"

        return [self.sample_csv, self.amplicon_csv, self.expt_qc_json]

    def _run(self):
        """
        Run the quality control analysis

        """

        bedcov_df = pd.read_csv(self.bedcov_csv)

        qc = ExperimentQualityControl(bedcov_df, self.metadata.df)
        qc.sample_df.to_csv(self.sample_csv, index=False)
        qc.amplicon_df.to_csv(self.amplicon_csv, index=False)

        with open(self.expt_qc_json, "w") as file:
            json.dump(qc.expt_dict, file)

        with open(
            f"{self.output_dir}/{'PASS.log' if qc.expt_dict['expt_pass'] else 'FAIL.log'}", "w"
        ) as file:
            for key, value in qc.expt_dict.items():
                file.write(f"{key}: {round(value, 2)}\n")

    def _plot(self):
        pass
