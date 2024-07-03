class MetadataFormatError(Exception):
    """Error in format / contents of a metadata file"""

    pass


class BEDFormatError(Exception):
    """Error in the format or contents of a BED file"""

    pass


class ReferenceGenomeMissingError(Exception):
    """Reference genome has not been downloaded"""

    pass


class InputDoesNotExit(Exception):
    """Input necessary for a given analysis does not exist"""

    pass
