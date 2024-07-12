from typing import Dict
from enum import Enum
from abc import ABC, abstractmethod


# ================================================================
# Possible drug resistance predictions
#
# ================================================================


class ResistanceStatus(str, Enum):
    R = "Resistant"
    S = "Sensitive"
    U = "Undetermined"
    M = "Missing"


# ================================================================
# Interface for predicting drug resistance
#
# ================================================================


class ResistanceRule(ABC):
    def __init__(self):
        self.name = self._get_name()

    @abstractmethod
    def _get_name(self) -> str:
        pass

    @abstractmethod
    def predict(self, sample_gt: Dict[str, Dict[str, int]]) -> ResistanceStatus:
        pass


# ================================================================
# Concrete drug resistance prediction rules
#
# ================================================================

# ----------------------------------------------------------------
# Antimalarials for Plasmodium falciparum
# - Using MalariaGEN's Pf7 classification rules
# ----------------------------------------------------------------


class Chloroquine(ResistanceRule):
    """
    Predict resistance to Chloroquine
    """

    classifier = {
        0: ResistanceStatus.S,
        1: ResistanceStatus.U,
        2: ResistanceStatus.R,
        -1: ResistanceStatus.M,
    }

    def _get_name(self):
        return "CQ-R"

    def predict(self, sample_gt: Dict):
        gt = sample_gt.get("crt").get("K76T")

        return self.classifier[gt]


class Pyrimethamine(ResistanceRule):
    """
    Predict resistance to Pyrimethamine monotherapy
    """

    classifier = {
        0: ResistanceStatus.S,
        1: ResistanceStatus.U,
        2: ResistanceStatus.R,
        -1: ResistanceStatus.M,
    }

    def _get_name(self):
        return "P-R"

    def predict(self, sample_gt: Dict):

        gt = sample_gt.get("dhfr").get("S108N")

        return self.classifier[gt]


class Sulfadoxine(ResistanceRule):
    """
    Predict resistance to Sulfadoxine monotherapy
    """

    classifier = {
        0: ResistanceStatus.S,
        1: ResistanceStatus.U,
        2: ResistanceStatus.R,
        -1: ResistanceStatus.M,
    }

    def _get_name(self):
        return "S-R"

    def predict(self, sample_gt: Dict):
        gt = sample_gt.get("dhps").get("A437G")

        return self.classifier[gt]


class SulfadoxinePyrimethamine(ResistanceRule):
    """
    Predict resistance to SP treatment

    TODO:
    - Why, according to MalariaGEN, is only DHFR required?
    - Confusing because can be sensitive to S but resistant to combination
    """

    RESISTANCE_MUTS = ["N51I", "C59R", "S108N"]

    def _get_name(self):
        return "SP-R"

    def predict(self, sample_gt: Dict):
        # Collect all genotypes
        gts = [sample_gt.get("dhfr").get(m) for m in self.RESISTANCE_MUTS]

        # Prioritise reporting missing
        if (-1 in gts) or (None in gts):
            return ResistanceStatus.M
        if 0 in gts:
            return ResistanceStatus.S
        if 1 in gts:
            return ResistanceStatus.U
        if sum(gts) == 2 * len(gts):
            return ResistanceStatus.R

        raise RuntimeError(f"Cannot determine resistance status for {self._get_name()}")


class Artemisinin(ResistanceRule):
    """
    Define calling mechanism for resistance to artemisinin monotherapy using WHO mutations

    """

    VALIDATED = [
        "F446I",
        "N458Y",
        "M476I",
        "Y493H",
        "R539T",
        "I543T",
        "P553L",
        "R561H",
        "P574L",
        "C580Y",
    ]
    CANDIDATE = [
        "P441L",
        "G449A",
        "C469F",
        "C469Y",
        "A481V",
        "R515K",
        "P527H",
        "N537I",
        "N537D",
        "G538V",
        "V568G",
        "R622I",
        "A675V",
    ]
    POTENTIAL = ["K479I", "G533A", "R575K", "M579I", "D584V", "P667T", "F673I", "H719N"]

    RESISTANCE_MUTS = VALIDATED + CANDIDATE

    def _get_name(self):
        return "ART-R"

    def predict(self, sample_gt: Dict):
        # Collect all genotypes
        gts = [sample_gt.get("kelch13").get(m) for m in self.RESISTANCE_MUTS]
        print(gts)

        # Prioritise reporting missing
        if 2 in gts:
            return ResistanceStatus.R
        if sum(gts) == 0:
            return ResistanceStatus.S
        if 1 in gts:
            return ResistanceStatus.U
        if all(gt == -1 for gt in gts):
            return ResistanceStatus.M

        raise RuntimeError(f"Cannot determine resistance status for {self._get_name()}")


# ================================================================
# Drug collections
#
# ================================================================


ANTIMALARIALS = [
    Chloroquine(),
    Pyrimethamine(),
    Sulfadoxine(),
    SulfadoxinePyrimethamine(),
]  # , Artemisinin()]
