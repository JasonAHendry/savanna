import os
import re
import logging
import subprocess
import pandas as pd
from dataclasses import dataclass
from savanna.util.dirs import produce_dir
from savanna.download.references import Reference


log = logging.getLogger()


# ================================================================
# Parsing output strings from `bcftools csq`
#
# ================================================================


@dataclass
class Consequence:
    csq: str
    target_gene: str
    transcript: str
    biotype: str

    # Optional, assigned only if within coding sequence
    strand: str = None
    aa_change: str = None
    nt_change: str = None

    # Optional, assigned only if aa_change exists
    concise_aa_change: str = None
    aa_pos: int = -1  # too keep column int in pandas

    def __post_init__(self):
        """
        Clean `aa_change` and `aa_pos`

        """
        if self.aa_change is None:
            return

        AAs = "ARNDCEQGHILKMFPSTWYV"
        stop = "\*"
        match = re.match(
            f"([0-9]+)([{AAs}|{stop}])(?:>[0-9]+([{AAs}|{stop}]))?", self.aa_change
        )

        if match is None:
            log.warn(f"Unable to parse AA change for: {self.aa_change}.")
            return

        aa_pos, from_aa, to_aa = match.groups()
        self.aa_pos = int(aa_pos)

        # Handle synonymous case
        if to_aa is None:
            to_aa = from_aa

        self.concise_aa_change = f"{from_aa}{aa_pos}{to_aa}"

    @classmethod
    def from_string(cls, csq_string: str):
        """
        Parse from the output string

        What to do if "double"?
        Or @
        Or *

        """
        if csq_string == ".":  # intergenic
            return cls(".", ".", ".", ".")

        if csq_string.startswith("@"):  # compound variety, recorded elsewhere
            return cls(".", ".", ".", f"compound{csq_string}")

        consequences = csq_string.split(",")
        if len(consequences) > 1:
            log.warn(
                f"Found multiple consequences of variant: {csq_string}! Keeping only first."
            )

        fields = consequences[0].split("|")
        assert len(fields) >= 4, f"Failed for {csq_string}"
        assert len(fields) <= 7, f"Failed for {csq_string}"

        return cls(*fields)


# ================================================================
# Annotating a VCF with information about amplicons & effects
#
# ================================================================


class VariantAnnotator:
    AMP_HEADER = (
        "##INFO=<ID=AMP_ID,Number=1,Type=String,Description=Amplicon identifier>"
    )

    def __init__(
        self, vcf_path: str, bed_path: str, reference: Reference, output_dir: str
    ) -> None:
        self.vcf_path = vcf_path
        self.bed_path = bed_path
        self.reference = reference  # The GFF path needs to be GFF3 compliant
        self.output_dir = produce_dir(output_dir)

        # Outputs
        self.output_vcf = f"{output_dir}/{os.path.basename(self.vcf_path).replace('.vcf.gz', '.annotated.vcf.gz')}"
        self.output_tsv = self.output_vcf.replace(".vcf.gz", ".tsv")

    def _get_wsaf_command(self, input_vcf: str = "-", output_vcf: str = "") -> str:
        """
        Compute the WSAF for each variant based on allelic depths
        """
        cmd = f"bcftools +fill-tags"
        if output_vcf:
            cmd += f" -Oz -o {output_vcf}"
        cmd += f" {input_vcf}"
        cmd += " -- -t FORMAT/WSAF=1-FORMAT/AD/FORMAT/DP"

        return cmd

    def _get_annotate_command(self, input_vcf: str = "-", output_vcf: str = "") -> str:
        """
        Create a string representing command required to annotate variants with
        their amplicon position
        """
        cmd = "bcftools annotate"
        cmd += f" -a {self.bed_path}"
        cmd += " -c CHROM,FROM,TO,AMP_ID"
        cmd += f" -H '{self.AMP_HEADER}'"
        cmd += " -Oz"
        if output_vcf:
            cmd += f" -o {output_vcf}"
        cmd += f" {input_vcf}"

        return cmd

    def _get_csq_command(self, input_vcf: str = "-", output_vcf: str = "") -> str:
        """
        Create a string representing command required
        to compute variant consequences
        """
        cmd = "bcftools csq"
        cmd += f" -f {self.reference.fasta_path}"
        cmd += f" -g {self.reference.gff_standard_path}"
        cmd += " --phase a"
        cmd += " -Oz"
        if output_vcf:
            cmd += f" -o {output_vcf}"
        cmd += f" {input_vcf}"

        return cmd

    def run(self):
        """
        Annotate variant calls using `bcftools csq`
        """
        cmd_tags = self._get_wsaf_command(
            input_vcf=self.vcf_path,
        )
        cmd_annot = self._get_annotate_command(
            input_vcf="-",
        )
        cmd_csq = self._get_csq_command(input_vcf="-", output_vcf=self.output_vcf)
        cmd = f"{cmd_tags} | {cmd_annot} | {cmd_csq}"
        subprocess.run(cmd, shell=True, check=True)

    def _convert_to_tsv(self):
        """
        Convert the annotated VCF file to a small TSV file
        in preparation for plotting
        """

        # Define fixed and per-sample fields
        # Note that sample name also gets added.
        fixed = {
            "chrom": "CHROM",
            "pos": "POS",
            "ref": "REF",
            "alt": "ALT",
            "qual": "QUAL",
            "consequence": "BCSQ",
            "amplicon": "AMP_ID",
        }
        called = {"gt": "GT", "gq": "GQ", "dp": "DP", "wsaf": "WSAF{0}"}

        # Write header
        sep = "\t"
        cmd_header = f"printf 'sample\t{sep.join(list(fixed) + list(called))}\n' > {self.output_tsv}"
        sep += "%"

        # Iterate and query for each sample
        cmd_query = f"for sample in `bcftools query {self.output_vcf} -l`; do"
        cmd_query += "  bcftools query -s $sample"
        cmd_query += f' -f "$sample\t%{sep.join(fixed.values())}\t[%{sep.join(called.values())}]\n"'
        cmd_query += f" {self.output_vcf} >> {self.output_tsv};"
        cmd_query += " done;"

        cmd = f"{cmd_header} && {cmd_query}"
        subprocess.run(cmd, shell=True, check=True)

    def _parse_consequences(self):
        """
        Parse the consequenc string in the TSV
        """
        df = pd.read_csv(self.output_tsv, sep="\t")
        csqs = [Consequence.from_string(c) for c in df["consequence"]]
        if csqs:
            mut_type, aa_change, aa_pos, strand = zip(
                *[(c.csq, c.concise_aa_change, c.aa_pos, c.strand) for c in csqs]
            )
        else:
            print(f"No mutations passed quality control for {self.barcode_name}.")
            mut_type, aa_change, aa_pos, strand = None, None, None, None

        df.insert(6, "mut_type", mut_type)
        df.insert(7, "aa_change", aa_change)
        df.insert(8, "aa_pos", aa_pos)
        df.insert(9, "strand", strand)
        df.drop("consequence", axis=1, inplace=True)
        df.to_csv(self.output_tsv, sep="\t", index=False)

    def convert_to_tsv(self):
        """
        This is more preparing to merge across barcodes
        """
        self._convert_to_tsv()
        self._parse_consequences()
