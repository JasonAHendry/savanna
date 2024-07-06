from .ento import EntoPipeline
from .plasmo import PlasmoPipeline

PIPELINE_COLLECTION = {
    "plasmo": PlasmoPipeline,
    "ento": EntoPipeline,
}

__all__ = [PIPELINE_COLLECTION]
