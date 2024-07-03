import json
import pandas as pd

import matplotlib.pyplot as plt

from .interface import ExperimentAnalysis
from savanna.run.barcode.interface import BarcodeAnalysis
from savanna.run.barcode.map import (
    BarcodeMapToPfalciparum,
    BarcodeMapUnmappedToHSapiens,
)


class ExperimentMapToPfalciparum(ExperimentAnalysis):
    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return BarcodeMapToPfalciparum(barcode_name, self.expt_dirs)

    def _summarise(self):
        pass

    def _plot(self):
        pass


class ExperimentMapUnmappedToHSapiens(ExperimentAnalysis):
    """
    In this case, we really are just wrapping the mapping steps
    We don't produce any summary files or plots of the BAM files themselves

    """

    def _get_barcode_analysis(self, barcode_name: str) -> BarcodeAnalysis:
        return BarcodeMapUnmappedToHSapiens(barcode_name, self.expt_dirs)

    def _summarise(self):
        pass

    def _plot(self):
        pass


# class ExperimentMapToPfalciparum(ExperimentAnalysis):
#     """
#     In this case, we really are just wrapping the mapping steps
#     We don't produce any summary files or plots of the BAM files themselves

#     I believe _run() can be metaprogrammed
#     _run(BarcodeAnalysis)

#     Then _summarise() should only operate for thoose where
#     run() has returned successfully.

#     """

#     def _run(self):
#         self.outputs = []
#         for b in self.metadata.barcodes:
#             analysis = BarcodeMapToPfalciparum(b, self.expt_dirs)
#             success = analysis.run()
#             if not success:
#                 print(f"MAPPING FAILED FOR {b}!")

#     def _summarise(self):
#         pass

#     def _plot(self):
#         pass


# class ExperimentMapUnmappedToHSapiens(ExperimentAnalysis):
#     """
#     In this case, we really are just wrapping the mapping steps
#     We don't produce any summary files or plots of the BAM files themselves

#     """
#     def _run(self):
#         self.outputs = []
#         for b in self.metadata.barcodes:
#             analysis = BarcodeMapUnmappedToHSapiens(b, self.expt_dirs)
#             success = analysis.run()
#             if not success:
#                 print(f"MAPPING FAILED FOR {b}!")

#     def _summarise(self):
#         pass

#     def _plot(self):
#         pass
