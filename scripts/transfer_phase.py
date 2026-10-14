"""Copy phased genotypes from a phased VCF onto the matching records of a gVCF.

Phasing tools (WhatsHap, longphase) write variant-only VCFs, while PhaseSpace
needs the gVCF's reference blocks for evidence. This keeps every gVCF record
and only rewrites GT order, the phased flag and PS where the phased VCF has the
same position, REF and called allele sequences for the sample. Allele indexes
may differ between the files (for example after a trio merge), so alleles are
matched by sequence.

    python scripts/transfer_phase.py SAMPLE calls.g.vcf.gz phased.vcf.gz out.g.vcf.gz
"""
import sys

import pysam


def main(sample, gvcf_path, phased_path, out_path):
    phase = {}
    with pysam.VariantFile(phased_path) as phased:
        for record in phased:
            call = record.samples[sample]
            gt = call.get("GT")
            if not call.phased or gt is None or None in gt or gt[0] == gt[1]:
                continue
            phase[(record.contig, record.pos, record.ref)] = (tuple(record.alleles[i] for i in gt), call.get("PS"))
    applied = 0
    with pysam.VariantFile(gvcf_path) as reader:
        header = reader.header
        if "PS" not in header.formats:
            header.add_meta("FORMAT", items=[("ID", "PS"), ("Number", 1), ("Type", "Integer"), ("Description", "Phase set")])
        with pysam.VariantFile(out_path, "wz", header=header) as writer:
            for record in reader:
                call = record.samples[sample]
                gt = call.get("GT")
                match = phase.get((record.contig, record.pos, record.ref))
                if match and gt is not None and None not in gt:
                    ordered, ps = match
                    called = [record.alleles[i] for i in gt]
                    if sorted(called) == sorted(ordered):
                        call["GT"] = tuple(record.alleles.index(allele) for allele in ordered)
                        call.phased = True
                        if ps is not None:
                            call["PS"] = ps
                        applied += 1
                writer.write(record)
    pysam.tabix_index(out_path, preset="vcf", force=True)
    print(f"Applied phase to {applied} of {len(phase)} phased heterozygous calls.")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
