from .ento import EntoPipeline
from .plasmo import PlasmoPipeline
from .plasmohuman import PlasmoHumanPipeline

PIPELINE_COLLECTION = {
    "plasmo": PlasmoPipeline,
    "ento": EntoPipeline,
    "plasmohuman": PlasmoHumanPipeline
}

__all__ = [PIPELINE_COLLECTION]
