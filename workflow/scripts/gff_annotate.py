from typing import Tuple

from pysam import (
    VariantFile,
    TabixFile,
    asGFF3,
    VariantRecord,
)
import re


def parse_attribute(s: str) -> Tuple[str, str]:
    t = s.split("=", maxsplit=1)
    return t[0], t[1]


BND_RE = re.compile(r""".*([]\[])((?P<seqname>.+):(?P<position>[0-9]+))([]\[])?.*""")


def parse_bnd_alt(s: str) -> Tuple[str, int]:
    return BND_RE.search(s)["seqname"], int(BND_RE.search(s)["position"])


gff_in: TabixFile = TabixFile(snakemake.input["annotation"], parser=asGFF3())

with VariantFile(snakemake.input["variants"]) as bcf_in:
    header = bcf_in.header
    info_genes = r"""##INFO=<ID=GENES,Number=.,Type=String,Description="genes overlapping the variant region">"""
    header.add_line(info_genes)
    info_exons = r"""##INFO=<ID=NUM_EXONS,Number=1,Type=Integer,Description="number of exons contained completely in the variant region">"""
    header.add_line(info_exons)

    with VariantFile(snakemake.output["variants"], "w", header=header) as bcf_out:
        for record in bcf_in:
            record: VariantRecord = record
            start, end = record.start, record.stop
            if "BND" in record.info["SVTYPE"]:
                mate_chrom, mate_pos = parse_bnd_alt(record.alts[0])
                if mate_chrom == record.chrom:
                    end = mate_pos
            start, end = min(start, end), max(start, end)
            num_exons = 0
            try:
                genes = set()
                exons = set()
                try:
                    for annot in gff_in.fetch(f"chr{record.chrom}", start, end):
                        attributes = dict(map(parse_attribute, annot.attributes.split(";")))
                        genes.add(attributes["gene_name"])
                        if "exon_number" in attributes and (
                            annot.start >= start and annot.end <= end
                        ):
                            exons.add(attributes["exon_number"])
                except ValueError:
                    for annot in gff_in.fetch(f"{record.chrom}", start, end):
                        attributes = dict(map(parse_attribute, annot.attributes.split(";")))
                        genes.add(attributes["gene_name"])
                        if "exon_number" in attributes and (
                            annot.start >= start and annot.end <= end
                        ):
                            exons.add(attributes["exon_number"])
                if len(genes) > 0:
                    record.info["GENES"] = list(genes)
                record.info["NUM_EXONS"] = len(set(exons))
            except ValueError:
                # Usually means: no annotations for a specific region
                pass
            finally:
                bcf_out.write(record)
