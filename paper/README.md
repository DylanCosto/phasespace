# JOSS submission checklist

Draft paper for the [Journal of Open Source Software](https://joss.theoj.org). Before submitting:

- [ ] Replace the ORCID and affiliation placeholders in `paper.md`, and write the acknowledgements.
- [ ] Check every reference in `paper.bib` (author lists, volumes, DOIs).
- [ ] Decide whether to cite GEMPROT (PubMed 30508040) and ProHap Explorer as related work; their references are not yet in `paper.bib`.
- [ ] Release a version (tag `v0.2.0`), archive it on Zenodo, and add the DOI to `CITATION.cff`.
- [ ] Preview the PDF: the [JOSS GitHub Action](https://github.com/marketplace/actions/open-journals-pdf-generator) or `docker run --rm -v "$PWD:/data" openjournals/inara -o pdf paper/paper.md`.
- [ ] Submit at https://joss.theoj.org/papers/new.
