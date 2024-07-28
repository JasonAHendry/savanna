import re
import warnings
import pandas as pd
from typing import List
from .exceptions import MetadataFormatError


def get_csv_delimiter(csv_path: str, delimiters: List[str] = [",", ";"]):
    """
    Determine which delimiter is being used in a CSV

    I couldn't find a solution just using `pd.read_csv`; using
    a regex sep=[,;] does not work.

    """

    with open(csv_path, "r") as csv:
        for header in csv:
            break

        used = [delimiter for delimiter in delimiters if header.find(delimiter) != -1]

        if not used:
            return None

        if not len(used) == 1:
            raise MetadataFormatError(
                f"Found multiple delimiters ({' '.join(used)}) in header: {header.strip()}."
            )

    return used[0]


def check_barcode_format(barcode: str, try_to_fix: bool = True) -> str:
    """
    Check that the format of a barcode is as expected, and optionally
    try and fix if it is not

    """

    if barcode == "unclassified":
        return barcode

    # Some hard-coded settings
    MAX_BARCODE = 96
    EXPECTED = "barcode[0-9]{2}$"
    EXAMPLE = "barcode01"

    if not isinstance(barcode, str):
        barcode = str(barcode)

    if not re.match(EXPECTED, barcode):

        if not try_to_fix:
            raise MetadataFormatError(
                f"Barcode '{barcode}' has bad format: must conform to '{EXAMPLE}'."
            )

        # Raise a warning
        warnings.warn(
            f"Barcode '{barcode}' has bad format: must conform to '{EXAMPLE}'. Trying to fix..."
        )

        nums = re.findall("[0-9]+", barcode)

        if not nums:
            raise MetadataFormatError(
                f"Barcode '{barcode}' has bad format: must conform to '{EXAMPLE}'."
            )

        if len(nums) > 1:
            raise MetadataFormatError(f"Multiple numbers found in barcode: {barcode}.")

        barcode_int = int(nums[0])
        if barcode_int > 96:
            raise MetadataFormatError(
                f"Barcode '{barcode}' exceeds maximum of {MAX_BARCODE}."
            )

        barcode = f"barcode{barcode_int:02d}"

    return barcode


class MetadataTableParser:
    """
    Parse the `metadata_csv` table, and make sure that it is formatted
    correctly

    """

    REQUIRED_COLUMNS = ["barcode", "sample_id"]
    UNIQ_COLUMNS = ["barcode"] # allowing duplicated 'sample_id' for now.

    def __init__(self, metadata_csv: str, include_unclassified: bool = False):
        """
        Load and sanity check the metadata table

        """

        self.csv = metadata_csv
        self.df = pd.read_csv(self.csv, delimiter=get_csv_delimiter(self.csv))

        self._check_for_columns()
        self._check_entries_unique()
        self._check_all_barcodes()

        self.barcodes = self.df["barcode"].tolist()
        if include_unclassified:
            self.barcodes.append("unclassified")

    @property
    def n_barcodes(self):
        return  len(self.barcodes)

    def _check_for_columns(self):
        """
        Check the correct columns are present

        """

        for c in self.REQUIRED_COLUMNS:
            if not c in self.df.columns:
                raise MetadataFormatError(f"Metadata must contain column called {c}!")

    def _check_entries_unique(self):
        """
        Check entires of the required columns are unique

        """

        for c in self.UNIQ_COLUMNS:
            all_entries = self.df[c].tolist()
            observed_entries = []
            for entry in all_entries:
                if entry in observed_entries:
                    raise MetadataFormatError(
                        f"Column {c} must contain only unique entires, but {entry} is duplicated."
                    )
                observed_entries.append(entry)

    def _check_all_barcodes(self) -> List[str]:
        self.df["barcode"] = [check_barcode_format(b) for b in self.df["barcode"]]
