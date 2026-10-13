"""Generate a tiny, wholly synthetic diploid locus with no personal data."""
from pathlib import Path

import pysam


def write_demo(directory: str | Path) -> dict[str, Path]:
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    paths = {"fasta": directory / "reference.fa", "annotation": directory / "transcripts.gtf",
             "vcf": directory / "sample.g.vcf.gz", "plain_vcf": directory / "sample.vcf"}
    if any(p.exists() for p in paths.values()) or (directory / "sample.g.vcf").exists():
        raise FileExistsError("Demo input directory already contains PhaseSpace inputs; choose an empty directory.")
    codons = ["ATG"] + ["GCT"] * 30 + ["TAA"]
    codons[7] = "CAA"
    codons[20] = "TGG"
    cds = "".join(codons)
    genome = list("N" * 400)
    intervals = [(100, 136), (200, 230), (300, 330)]
    offset = 0
    for start, end in intervals:
        genome[start:end] = cds[offset:offset + end - start]
        offset += end - start
    sequence = "".join(genome)
    paths["fasta"].write_text(">chrDemo\n" + "\n".join(sequence[i:i + 60] for i in range(0, 400, 60)) + "\n")
    pysam.faidx(str(paths["fasta"]))
    paths["annotation"].write_text("".join(
        f'chrDemo\tPhaseSpace\tCDS\t{start + 1}\t{end}\t.\t+\t0\tgene_id "DEMO1"; transcript_id "DEMO1.1"; gene_name "DEMO1";\n'
        for start, end in intervals))
    header = '''##fileformat=VCFv4.3
##source=PhaseSpaceSyntheticDemo
##contig=<ID=chrDemo,length=400>
##ALT=<ID=NON_REF,Description="Any other allele">
##INFO=<ID=END,Number=1,Type=Integer,Description="Inclusive block end">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum block depth">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC
'''
    variant_rows = [
        "chrDemo\t122\tstop_A\tC\tT\t60\tPASS\t.\tGT:GQ:DP\t0/1:60:30\n",
        "chrDemo\t226\tstop_B\tG\tA\t60\tPASS\t.\tGT:GQ:DP\t0/1:60:30\n",
    ]
    paths["plain_vcf"].write_text(header + "".join(variant_rows))
    blocks = [(100, 121), (122, 136), (200, 225), (226, 230), (300, 310), (316, 330)]
    rows = [(int(row.split('\t')[1]), row) for row in variant_rows]
    for start, end in blocks:
        rows.append((start + 1, f"chrDemo\t{start + 1}\t.\t{sequence[start]}\t<NON_REF>\t.\tPASS\tEND={end}\tGT:GQ:DP:MIN_DP\t0/0:60:30:28\n"))
    raw = directory / "sample.g.vcf"
    raw.write_text(header + "".join(row for _, row in sorted(rows)))
    pysam.tabix_compress(str(raw), str(paths["vcf"]), force=False)
    pysam.tabix_index(str(paths["vcf"]), preset="vcf", force=False)
    return paths
