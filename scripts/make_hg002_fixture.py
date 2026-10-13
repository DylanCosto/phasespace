"""Cut small, real HG002 regression fixtures from public chr20 data.

The fixture keeps true GRCh38 chr20 coordinates. The reference is chr20 with
every base outside the selected transcripts' coding spans masked to N, so it
compresses to a few hundred kilobytes. Sources are listed in the generated
README; see docs/validation.md for how they were obtained.

    python scripts/make_hg002_fixture.py \\
        --fasta GRCh38.chr20.fa --gtf gencode.v46.annotation.gtf.gz \\
        --gff3 Homo_sapiens.GRCh38.112.chromosome.20.gff3.gz \\
        --vcf pacbio=pacbio-chr20.HG002.g.vcf.gz --vcf deepvariant=dv.HG002.g.vcf.gz \\
        --vcf giab=giab-v4.2.1.HG002.chr20.vcf.gz --out tests/data/hg002_chr20
"""
import argparse
import gzip
from pathlib import Path
import re

import pysam

# One MANE Select transcript per gene, each chosen for a behaviour seen in real data.
TRANSCRIPTS = {
    "ENST00000371998.8": "NCOA3: heterozygous in-frame deletion, forward strand",
    "ENST00000372726.5": "TP53TG5: homozygous in-frame deletion, reverse strand",
    "ENST00000245960.10": "CDC25B: three unphased synonymous SNVs; CDS changes but protein never does",
    "ENST00000217428.7": "SPINT3: three unphased SNVs including missense, reverse strand",
    "ENST00000382398.4": "DEFB126: homozygous frameshift that reads through the stop codon",
}
PAD = 100


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _attribute(text, key):
    match = re.search(rf'{key}[ =]"?([^";]+)"?', text)
    return match.group(1) if match else None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--gff3", required=True)
    parser.add_argument("--vcf", action="append", required=True, help="label=path")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    gtf_lines, gene_ids, spans = [], set(), {}
    for line in _open(args.gtf):
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        transcript = _attribute(fields[8], "transcript_id")
        if transcript in TRANSCRIPTS:
            gtf_lines.append(line)
            gene_ids.add(_attribute(fields[8], "gene_id"))
            if fields[2] in ("CDS", "stop_codon"):
                left, right = spans.get(transcript, (int(fields[3]) - 1, int(fields[4])))
                spans[transcript] = (min(left, int(fields[3]) - 1), max(right, int(fields[4])))
    gene_lines = [line for line in _open(args.gtf) if not line.startswith("#") and line.split("\t")[2] == "gene" and _attribute(line.split("\t")[8], "gene_id") in gene_ids]
    (args.out / "gencode.v46.gtf").write_text("".join(gene_lines + gtf_lines))

    unversioned = {name.split(".")[0] for name in TRANSCRIPTS}
    genes = {gene.split(".")[0] for gene in gene_ids}
    gff3 = ["##gff-version 3\n"]
    for line in _open(args.gff3):
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        identifier = _attribute(fields[8], "ID") or ""
        parent = _attribute(fields[8], "Parent") or ""
        if identifier.split(":")[-1] in genes | unversioned or parent.split(":")[-1] in unversioned:
            gff3.append("chr" + line if not line.startswith("chr") else line)
    (args.out / "ensembl112.gff3").write_text("".join(gff3))

    regions = sorted((left - PAD, right + PAD) for left, right in spans.values())
    # Records may start before a region (gVCF blocks); their REF bases stay unmasked.
    ref_spans = list(regions)
    for item in args.vcf:
        label, path = item.split("=", 1)
        with pysam.VariantFile(path) as reader:
            target = args.out / f"{label}.HG002.vcf.gz"
            with pysam.VariantFile(str(target), "wz", header=reader.header) as writer:
                seen = set()
                for left, right in regions:
                    for record in reader.fetch("chr20", left, right):
                        key = (record.pos, record.ref, record.alts, record.stop)
                        if key not in seen:
                            seen.add(key)
                            writer.write(record)
                            ref_spans.append((record.start, record.start + len(record.ref)))
        pysam.tabix_index(str(target), preset="vcf", force=True)

    with pysam.FastaFile(args.fasta) as fasta:
        sequence = list("N" * fasta.get_reference_length("chr20"))
        for left, right in ref_spans:
            sequence[left:right] = fasta.fetch("chr20", left, right).upper()
    plain = args.out / "GRCh38.chr20.masked.fa"
    sequence = "".join(sequence)
    plain.write_text(">chr20\n" + "\n".join(sequence[i:i + 60] for i in range(0, len(sequence), 60)) + "\n")
    pysam.tabix_compress(str(plain), str(plain) + ".gz", force=True)
    plain.unlink()
    pysam.faidx(str(plain) + ".gz")

    (args.out / "README.md").write_text(
        "# HG002 chr20 regression fixture\n\n"
        "Real GIAB HG002 calls for five GRCh38 chr20 transcripts, cut by `scripts/make_hg002_fixture.py`.\n"
        "Coordinates are unchanged; reference bases outside the coding spans (plus 100 bp) are masked to N.\n\n"
        "| Transcript | Why it is here |\n|---|---|\n"
        + "".join(f"| {name} | {why} |\n" for name, why in TRANSCRIPTS.items())
        + "\nSources (all public):\n\n"
        "- `pacbio.HG002.vcf.gz`: DeepVariant 1.10 DeepTrio PacBio chr20 gVCF, "
        "`gs://deepvariant/case-study-outputs/1.10.0/deeptrio/pacbio-chr20/HG002.child.g.vcf.gz`.\n"
        "- `deepvariant.HG002.vcf.gz`: DeepVariant 1.10 DeepTrio gVCF published under "
        "`case-study-outputs/1.10.0/deeptrio/ont-r104-chr20/HG002.child.g.vcf.gz`. Despite the folder name, the file "
        "also holds exome-scale calls on other chromosomes, so its sequencing platform is not asserted here.\n"
        "- `giab.HG002.vcf.gz`: GIAB HG002 NISTv4.2.1 GRCh38 benchmark VCF (unphased, variant-only).\n"
        "- `GRCh38.chr20.masked.fa.gz`: chr20 of the GRCh38 no-alt analysis set used by DeepVariant.\n"
        "- `gencode.v46.gtf`: GENCODE v46 comprehensive annotation. `ensembl112.gff3`: Ensembl 112 chr20 GFF3, "
        "contig renamed from `20` to `chr20`.\n"
    )


if __name__ == "__main__":
    main()
