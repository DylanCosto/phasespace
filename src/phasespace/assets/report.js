/* All metadata and sequence content enters the document as text, never markup. */
(function () {
  "use strict";

  // NCBI translation table 1, matching analysis.py; codons with N translate to X.
  const BASES = "TCAG";
  const AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  const CODONS = {};
  for (let i = 0; i < 64; i += 1) CODONS[BASES[i >> 4] + BASES[(i >> 2) & 3] + BASES[i & 3]] = AMINO_ACIDS[i];

  function filterScenarios(scenarios, constraints) {
    return scenarios.filter(function (scenario) {
      return constraints.every(function (constraint) {
        const a = scenario.orientations[constraint.a];
        const b = scenario.orientations[constraint.b];
        return (a === 0 || a === 1) && (b === 0 || b === 1) && (a ^ b) === constraint.parity;
      });
    });
  }

  function groupScenarios(scenarios) {
    const groups = new Map();
    scenarios.forEach(function (scenario) {
      if (!groups.has(scenario.outcome_key)) groups.set(scenario.outcome_key, []);
      groups.get(scenario.outcome_key).push(scenario);
    });
    return Array.from(groups.values());
  }

  /** Add a relationship; one already present with the same direction gains the
   *  new source, and one with the opposite direction is replaced. Returns a new list. */
  function mergeConstraint(constraints, constraint) {
    const same = function (c) { return (c.a === constraint.a && c.b === constraint.b) || (c.a === constraint.b && c.b === constraint.a); };
    const existing = constraints.find(same);
    if (!existing) return constraints.concat([{ a: constraint.a, b: constraint.b, parity: constraint.parity, sources: constraint.sources.slice() }]);
    return constraints.map(function (c) {
      if (!same(c)) return c;
      if (c.parity !== constraint.parity) return { a: constraint.a, b: constraint.b, parity: constraint.parity, sources: constraint.sources.slice() };
      return { a: c.a, b: c.b, parity: c.parity, sources: c.sources.concat(constraint.sources.filter(function (s) { return c.sources.indexOf(s) < 0; })) };
    });
  }

  /** Withdraw one source; relationships nothing else supports are dropped. Returns a new list. */
  function withdrawSource(constraints, source) {
    return constraints.map(function (c) { return { a: c.a, b: c.b, parity: c.parity, sources: c.sources.filter(function (s) { return s !== source; }) }; })
      .filter(function (c) { return c.sources.length; });
  }

  /** Arrangements grouped by their pair of predicted proteins, in first-seen order. */
  function groupByProtein(scenarios) {
    const groups = new Map();
    scenarios.forEach(function (scenario) {
      const key = scenario.protein_outcome_key || scenario.outcome_key;
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key).push(scenario);
    });
    return Array.from(groups.values());
  }

  function sharedProteinScale(scenarios, referenceLength) {
    return scenarios.reduce(function (longest, scenario) {
      return scenario.haplotypes.reduce(function (length, haplotype) {
        return Math.max(length, haplotype.protein_length);
      }, longest);
    }, Math.max(1, referenceLength || 0));
  }

  /** Translate complete codons up to, not including, the first stop. */
  function translate(cds) {
    let protein = "";
    for (let i = 0; i + 3 <= cds.length; i += 3) {
      const amino = CODONS[cds.slice(i, i + 3)] || "X";
      if (amino === "*") break;
      protein += amino;
    }
    return protein;
  }

  /** Rebuild one copy's CDS from the reference and its called edits. */
  function rebuildCds(reference, haplotype, variantLookup) {
    if (typeof haplotype.cds === "string") return haplotype.cds;
    const edits = (haplotype.variant_ids || []).map(function (id, index) {
      return { variant: variantLookup.get(id), allele: haplotype.alleles[index] };
    });
    edits.sort(function (a, b) { return b.variant.cds_start - a.variant.cds_start; });
    return edits.reduce(function (cds, edit) {
      return cds.slice(0, edit.variant.cds_start) + edit.allele.toUpperCase() + cds.slice(edit.variant.cds_end);
    }, reference);
  }

  /** Restore the complete analysis schema from a report's compact encoding. */
  function expandAnalysis(data) {
    const reference = data.transcript.reference_cds;
    const lookup = new Map((data.variants || []).map(function (variant) { return [variant.id, variant]; }));
    const expanded = Object.assign({}, data, {
      scenarios: (data.scenarios || []).map(function (scenario) {
        return Object.assign({}, scenario, {
          haplotypes: scenario.haplotypes.map(function (haplotype) {
            if (typeof haplotype.cds === "string") return haplotype;
            const cds = rebuildCds(reference, haplotype, lookup);
            const copy = {};
            Object.keys(haplotype).forEach(function (key) {
              copy[key] = haplotype[key];
              if (key === "alleles") { copy.cds = cds; copy.protein = translate(cds); }
            });
            return copy;
          })
        });
      })
    });
    delete expanded.report_encoding;
    return expanded;
  }

  function isReferenceCopy(haplotype, reference) {
    return (haplotype.consequences || []).indexOf("reference") >= 0 || (typeof haplotype.cds === "string" && haplotype.cds === reference);
  }

  const UNCHANGED_PROTEIN = ["reference", "synonymous", "protein_sequence_unchanged"];

  /** A copy whose predicted protein matches the reference (its CDS may differ synonymously). */
  function keepsReferenceProtein(haplotype, reference) {
    if (isReferenceCopy(haplotype, reference)) return true;
    const consequences = haplotype.consequences || [];
    return consequences.length > 0 && consequences.every(function (key) { return UNCHANGED_PROTEIN.indexOf(key) >= 0; });
  }

  /** How many compatible arrangements keep at least one copy with the reference protein. */
  function phaseVerdict(scenarios, reference) {
    const counts = scenarios.map(function (scenario) {
      return scenario.haplotypes.filter(function (haplotype) { return keepsReferenceProtein(haplotype, reference); }).length;
    });
    const keep = counts.filter(function (count) { return count > 0; }).length;
    const tone = !counts.length ? "empty" : keep === counts.length ? "good" : keep === 0 ? "none" : "mixed";
    const cdsChanged = scenarios.some(function (scenario) {
      return scenario.haplotypes.some(function (haplotype) { return !isReferenceCopy(haplotype, reference); });
    });
    return { tone: tone, total: counts.length, keep: keep, min: Math.min.apply(null, counts), max: Math.max.apply(null, counts), cdsChanged: cdsChanged };
  }

  /** The arrangement reached by moving one block's alleles to the other copy. */
  function flipTarget(scenarios, blockIds, current, blockId) {
    if (blockIds.length < 2) return null;
    const target = Object.assign({}, current.orientations);
    // Block 1 is fixed to remove copy-label symmetry; flipping it equals flipping the rest.
    blockIds.forEach(function (id) {
      if (blockId === blockIds[0] ? id !== blockId : id === blockId) target[id] ^= 1;
    });
    return scenarios.find(function (scenario) {
      return blockIds.every(function (id) { return scenario.orientations[id] === target[id]; });
    }) || null;
  }

  /** A short coding-change label. Substitutions are exact; indels are not 3'-normalized. */
  function changeLabel(variant, allele) {
    const ref = variant.cds_ref.toUpperCase();
    const alt = allele.toUpperCase();
    const start = variant.cds_start;
    const span = function (a, b) { return a === b ? String(a) : a + "_" + b; };
    if (ref.length === alt.length) return ref.length === 1 ? "c." + (start + 1) + ref + ">" + alt : "c." + span(start + 1, start + ref.length) + "delins" + alt;
    if (alt.length < ref.length && ref.startsWith(alt)) return "c." + span(start + alt.length + 1, start + ref.length) + "del";
    if (alt.length < ref.length && ref.endsWith(alt)) return "c." + span(start + 1, start + ref.length - alt.length) + "del";
    if (alt.length > ref.length && alt.startsWith(ref)) return "c." + (start + ref.length) + "_" + (start + ref.length + 1) + "ins" + alt.slice(ref.length);
    if (alt.length > ref.length && alt.endsWith(ref)) return "c." + start + "_" + (start + 1) + "ins" + alt.slice(0, alt.length - ref.length);
    return "c." + span(start + 1, start + ref.length) + "delins" + alt;
  }

  const MAX_LABEL_ROWS = 3;

  /**
   * Assign each label (sorted by center) a row, or -1 when it cannot be placed.
   * Labels in a row keep a 6 px gap; a stem never crosses a label in a lower row,
   * and a label never sits across a stem rising from a lower pin to a higher row.
   */
  function placeLabels(items, width, maxRows) {
    const placed = [];
    return items.map(function (item) {
      const left = Math.max(0, Math.min(width - item.width, item.center - item.width / 2));
      const right = left + item.width;
      for (let row = 0; row < maxRows; row += 1) {
        const clash = placed.some(function (spot) {
          if (spot.row === row) return left < spot.right + 6 && spot.left < right + 6;
          if (spot.row < row) return item.center > spot.left - 3 && item.center < spot.right + 3;
          return spot.center > left - 3 && spot.center < right + 3;
        });
        if (!clash) {
          placed.push({ left: left, right: right, row: row, center: item.center });
          return { row: row, left: left };
        }
      }
      return { row: -1, left: left };
    });
  }

  // Export the exact report logic for independent tests.
  if (typeof module !== "undefined" && module.exports) {
    module.exports = {
      filterScenarios: filterScenarios, groupScenarios: groupScenarios, groupByProtein: groupByProtein, sharedProteinScale: sharedProteinScale,
      mergeConstraint: mergeConstraint, withdrawSource: withdrawSource,
      translate: translate, rebuildCds: rebuildCds, expandAnalysis: expandAnalysis, phaseVerdict: phaseVerdict,
      flipTarget: flipTarget, changeLabel: changeLabel, placeLabels: placeLabels
    };
  }
  if (typeof document === "undefined") return;

  const data = JSON.parse(document.getElementById("phasespace-data").textContent);
  const transcript = data.transcript;
  const reference = transcript.reference_cds;
  const blocks = data.blocks || [];
  const blockIds = blocks.map(function (block) { return block.id; });
  const scenarios = data.scenarios || [];
  const variants = data.variants || [];
  const spliceRegion = data.splice_region || [];
  const readEvidence = data.read_evidence || null;
  const parentEvidence = data.inheritance_evidence || null;
  const infoText = function (variant) { return Object.keys(variant.info || {}).map(function (key) { return key + "=" + variant.info[key]; }).join("  ·  "); };
  const originText = function (variant) {
    const origin = parentEvidence && parentEvidence.variants[variant.id];
    if (!origin) return "";
    if (origin.status === "informative") return "\nFrom father: " + variant.alleles[origin.paternal_allele] + " · from mother: " + variant.alleles[1 - origin.paternal_allele];
    return "\nParent of origin: " + ({ uninformative: "unclear (both parents heterozygous)", conflict: "Mendelian conflict", unknown: "parents not genotyped" }[origin.status] || origin.status);
  };
  const readText = function (variant) {
    const counts = readEvidence && readEvidence.variants[variant.id];
    if (!counts) return "";
    return "\nReads: " + variant.alleles.map(function (allele, index) { return (allele.length > 8 ? allele.slice(0, 7) + "…" : allele) + " ×" + counts.allele_reads[index]; }).join(" · ");
  };
  const hasInfo = variants.some(function (variant) { return Object.keys(variant.info || {}).length > 0; });
  const spliceIds = new Set(spliceRegion.map(function (hit) { return hit.id; }));
  const proteinScale = sharedProteinScale(scenarios, data.reference_protein_length);
  const variantLookup = new Map(variants.map(function (variant) { return [variant.id, variant]; }));
  const blockLookup = new Map(blocks.map(function (block) { return [block.id, block]; }));
  const scenarioNumber = new Map(scenarios.map(function (scenario, index) { return [scenario.id, index + 1]; }));
  const proteinKey = function (scenario) { return scenario.protein_outcome_key || scenario.outcome_key; };
  const outcomeNumber = new Map(groupByProtein(scenarios).map(function (group, index) { return [proteinKey(group[0]), index + 1]; }));
  const sequenceCache = new Map();
  const state = { constraints: [], selectedId: scenarios.length ? scenarios[0].id : null, compatible: [], showAllOutcomes: false };
  const OUTCOME_PREVIEW = 24;
  const byId = function (id) { return document.getElementById(id); };
  const number = function (value) { return Number(value).toLocaleString("en-US"); };
  const plural = function (count, word, many) { return number(count) + " " + (count === 1 ? word : (many || word + "s")); };
  const text = function (id, value) { byId(id).textContent = String(value); };
  const blockName = function (id) { const index = blockIds.indexOf(id); return index >= 0 ? "Block " + (index + 1) : String(id); };
  const SVG_NS = "http://www.w3.org/2000/svg";
  const CONSEQUENCES = {
    reference: ["Reference CDS", "good"], synonymous: ["Synonymous", "good"], protein_sequence_unchanged: ["Protein unchanged", "good"],
    stop_gained: ["Premature stop", "critical"], frameshift: ["Frameshift", "critical"], start_lost: ["Start lost", "critical"], stop_lost: ["Stop lost", "critical"],
    missense: ["Missense", "warning"], inframe_insertion: ["In-frame insertion", "warning"], inframe_deletion: ["In-frame deletion", "warning"],
    frame_restored: ["Frame restored", "warning"], uncertain_codon: ["Uncertain codon", "warning"]
  };
  const DISRUPTIVE = ["frameshift", "stop_gained", "start_lost", "stop_lost"];

  function element(tag, className, content) {
    const node = document.createElement(tag);
    if (className) node.className = className;
    if (content !== undefined && content !== null) node.textContent = String(content);
    return node;
  }

  function icon(paths, viewBox) {
    const svg = document.createElementNS(SVG_NS, "svg");
    svg.setAttribute("viewBox", viewBox || "0 0 24 24");
    svg.setAttribute("aria-hidden", "true");
    paths.forEach(function (d) { const path = document.createElementNS(SVG_NS, "path"); path.setAttribute("d", d); svg.append(path); });
    return svg;
  }
  const ICONS = {
    check: ["M5 12.5l4.5 4.5L19 7.5"],
    alert: ["M12 8v5", "M12 16.5v.5", "M10.3 3.9 2.6 17.5A2 2 0 0 0 4.3 20.5h15.4a2 2 0 0 0 1.7-3L13.7 3.9a2 2 0 0 0-3.4 0z"],
    split: ["M6 4v6a4 4 0 0 0 4 4h8", "M6 20v-6", "M15 11l3 3-3 3"],
    dot: ["M12 12h.01"],
    flip: ["M4 5h7l-2-2", "M12 11H5l2 2"]
  };

  function option(value, label) {
    const node = element("option", "", label);
    node.value = String(value);
    return node;
  }

  function tip(node, value) { node.dataset.tip = value; return node; }

  function consequenceLabel(key) { return (CONSEQUENCES[key] || [String(key).replace(/_/g, " ")])[0]; }
  function consequenceTone(key) { return (CONSEQUENCES[key] || [null, "neutral"])[1]; }

  /** Allele carried by one copy, from the analysis or (for older results) its orientation. */
  function alleleOnCopy(scenario, copyIndex, variant) {
    const haplotype = scenario.haplotypes[copyIndex];
    const index = (haplotype.variant_ids || []).indexOf(variant.id);
    if (index >= 0 && haplotype.alleles) return haplotype.alleles[index];
    const bit = variant.block ? scenario.orientations[variant.block] || 0 : 0;
    return variant.cds_alleles[copyIndex ^ bit];
  }

  function sequencesFor(haplotype) {
    if (!sequenceCache.has(haplotype)) {
      const cds = rebuildCds(reference, haplotype, variantLookup);
      sequenceCache.set(haplotype, { cds: cds, protein: typeof haplotype.protein === "string" ? haplotype.protein : translate(cds) });
    }
    return sequenceCache.get(haplotype);
  }

  /** Coding intervals in transcript (5'→3') order with their CDS offsets. */
  const exons = (function () {
    const ordered = (transcript.cds_intervals || []).slice();
    if (transcript.strand === "-") ordered.reverse();
    let offset = 0;
    return ordered.map(function (interval, index) {
      const exon = { index: index, left: interval[0], right: interval[1], start: offset, end: offset + interval[1] - interval[0] };
      offset = exon.end;
      return exon;
    });
  })();
  const cdsLength = Math.max(1, reference.length);
  const pct = function (position) { return Math.max(0, Math.min(100, 100 * position / cdsLength)); };

  function toCds(exon, a, b) {
    return transcript.strand === "-" ? [exon.start + exon.right - b, exon.start + exon.right - a] : [exon.start + a - exon.left, exon.start + b - exon.left];
  }

  /** Base counts for each evidence status, in CDS coordinates. */
  const evidence = (function () {
    const pieces = [];
    const counts = { reference: 0, variant: 0, unknown: 0 };
    (data.evidence || []).forEach(function (item) {
      const status = item.status === "reference" || item.status === "variant" ? item.status : "unknown";
      exons.forEach(function (exon) {
        const a = Math.max(exon.left, item.start);
        const b = Math.min(exon.right, item.end);
        if (b <= a) return;
        const range = toCds(exon, a, b);
        counts[status] += b - a;
        pieces.push({ start: range[0], end: range[1], status: status, detail: item.detail, genomic: [a, b] });
      });
    });
    const total = exons.reduce(function (sum, exon) { return sum + exon.end - exon.start; }, 0);
    // Missing evidence still counts as unknown, never as reference.
    counts.unknown = Math.max(counts.unknown, total - counts.reference - counts.variant);
    return { pieces: pieces, counts: counts, total: total };
  })();

  /* ---------- Header, verdict and stats ---------- */

  function renderIdentity() {
    text("gene-title", transcript.gene || transcript.id);
    text("transcript-id", transcript.id);
    text("sample-label", "Sample " + data.sample);
    const meta = byId("locus-meta");
    const start = Math.min.apply(null, exons.map(function (exon) { return exon.left; }));
    const end = Math.max.apply(null, exons.map(function (exon) { return exon.right; }));
    [transcript.contig + ":" + number(start + 1) + "–" + number(end), (transcript.strand === "-" ? "Reverse" : "Forward") + " strand",
      number(reference.length) + " coding bases", plural(exons.length, "coding exon")].forEach(function (value) {
      meta.append(element("span", "", value));
    });
    byId("synthetic-label").hidden = data.synthetic !== true;
    text("total-arrangements", number(scenarios.length));
    text("arrangements-note", blocks.length > 1 ? "ways to place " + number(blocks.length) + " independent blocks on two copies" : blocks.length ? "one block, so its phase is fixed" : "no heterozygous calls to place");
    text("total-protein-outcomes", number(outcomeNumber.size));
    const codingOutcomes = new Set(scenarios.map(function (scenario) { return scenario.outcome_key; })).size;
    text("protein-note", codingOutcomes > outcomeNumber.size ? "protein pairs; " + number(codingOutcomes) + " if synonymous differences count" : "different protein pairs across arrangements");
    const heterozygous = variants.filter(function (variant) { return variant.gt[0] !== variant.gt[1]; }).length;
    text("variant-stat", number(variants.length));
    text("variant-note", number(heterozygous) + " heterozygous · " + number(variants.length - heterozygous) + " on both copies");
    const called = evidence.counts.reference + evidence.counts.variant;
    const share = evidence.total ? Math.round(100 * called / evidence.total) : 0;
    text("coverage-stat", share + "%");
    byId("coverage-meter").style.width = share + "%";
    tip(byId("coverage-meter").parentElement.parentElement, number(called) + " of " + number(evidence.total) + " coding bases have a passing reference or variant call.");
    text("source-meta", (data.source ? "Source: " + data.source : "Source metadata unavailable") + (data.version ? "  ·  PhaseSpace " + data.version : ""));
    if (data.version) text("version-label", data.version);
    document.title = (transcript.gene || transcript.id) + " · " + data.sample + " · PhaseSpace";
    if ((data.warnings || []).length) {
      byId("warning-panel").hidden = false;
      data.warnings.forEach(function (warning) { byId("warnings").append(element("li", "", warning)); });
    }
  }

  function renderVerdict() {
    const verdict = phaseVerdict(state.compatible, reference);
    const box = byId("verdict");
    box.dataset.tone = verdict.tone;
    box.querySelector(".verdict-icon").replaceChildren(icon(verdict.tone === "good" ? ICONS.check : verdict.tone === "mixed" ? ICONS.split : ICONS.alert));
    const sources = Array.from(new Set([].concat.apply([], state.constraints.map(function (c) { return c.sources; }))));
    const linkedBy = Array.from(new Set([].concat.apply([], (data.linked_blocks || []).map(function (group) { return group.sources; }))));
    const base = linkedBy.length ? "With blocks linked by " + linkedBy.map(function (source) { return source === "parents" ? "inheritance" : source; }).join(" and ") : "From the calls alone";
    text("verdict-kicker", !sources.length ? base : "With " + sources.map(function (source) {
      return { reads: "read evidence", parents: "inheritance", hypothesis: "your hypotheses" }[source];
    }).join(" and ") + " applied");
    let title, detail;
    if (verdict.tone === "empty") {
      title = "No arrangement fits the applied relationships";
      detail = "They contradict each other. Remove one or reset to see supported arrangements again.";
    } else if (verdict.tone === "good" && verdict.min === 2) {
      title = "Both copies keep the reference protein sequence";
      detail = verdict.cdsChanged ? "The included calls change the coding sequence only synonymously." : "No included call changes the CDS of either copy.";
    } else if (verdict.tone === "good") {
      title = "Every arrangement keeps a copy with the reference protein";
      detail = verdict.total === 1 ? "One copy's predicted protein matches the reference." : "Whatever the phase, at least one copy's predicted protein matches the reference.";
    } else if (verdict.tone === "none") {
      title = "No copy keeps the reference protein";
      detail = verdict.total === 1 ? "Both copies carry protein-changing edits." : "Both copies carry protein-changing edits in all " + number(verdict.total) + " arrangements.";
    } else {
      title = "Phase decides whether a copy keeps the reference protein";
      detail = number(verdict.keep) + " of " + number(verdict.total) + " arrangements " + (verdict.keep === 1 ? "keeps" : "keep") + " one; the rest change the protein on both copies.";
    }
    text("verdict-title", title);
    text("verdict-detail", detail);
    const caveats = byId("verdict-caveats");
    caveats.replaceChildren();
    if (verdict.tone === "empty") return;
    const add = function (label, hint) { caveats.append(tip(element("li", "", label), hint)); };
    if (evidence.counts.unknown > 0) add(plural(evidence.counts.unknown, "coding base") + " without a call", "Results cover the included calls only; variants could hide where there is no reference or variant call.");
    if (spliceRegion.length) add(plural(spliceRegion.length, "splice-region variant") + " not modelled", "Within 3 exonic or 8 intronic bases of a coding-exon junction. Splicing effects are outside the model.");
    if ((data.linked_blocks || []).length) add("Blocks linked by evidence before analysis", "--link-evidence merged blocks using clear read or inheritance evidence; these arrangements assume it.");
  }

  /* ---------- Phase evidence (reads, parents) ---------- */

  const EVIDENCE_SOURCES = [
    { key: "reads", name: "Reads", data: readEvidence, verb: "read-supported",
      none: "No read spans two blocks",
      meta: function (e) { return (e.supported.length ? plural(e.supported.length, "clear relationship") + " · " : "") + plural(e.fragments, "informative fragment") + " · " + e.files.join(", "); },
      within: function (agree, total) { return "Within blocks, reads agree with the input phase in " + number(agree) + " of " + number(total) + " observations."; },
      rule: function (e) { return "Clear means at least " + e.thresholds.min_count + " fragments with " + Math.round(100 * e.thresholds.min_fraction) + "% agreeing."; } },
    { key: "parents", name: "Inheritance", data: parentEvidence, verb: "inherited", none: "No two blocks have calls of clear parental origin",
      meta: function (e) {
        const informative = Object.keys(e.variants).filter(function (id) { return e.variants[id].status === "informative"; }).length;
        return (e.supported.length ? plural(e.supported.length, "clear relationship") + " · " : "") + "parent of origin clear for " + plural(informative, "call") + (e.conflicts ? " · " + plural(e.conflicts, "Mendelian conflict") : "") + " · father " + e.father + ", mother " + e.mother;
      },
      within: function (agree, total) { return "Within blocks, inheritance agrees with the input phase for " + number(agree) + " of " + number(total) + " call pairs."; },
      rule: function () { return "Clear means every informative pair of calls agrees."; } }
  ];
  const VERDICT_PHRASES = { good: "A copy keeps the reference protein", none: "No copy keeps the reference protein",
    mixed: "Phase still decides", conflicting: "Its relationships contradict each other", empty: "No arrangement fits" };
  const VERDICT_TONES = { good: "good", none: "none", mixed: "mixed", conflicting: "none", empty: "none" };

  function pairRow(pair) {
    const total = pair.same + pair.opposite;
    const direction = pair.same >= pair.opposite ? "Same" : "Opposite";
    const majority = Math.max(pair.same, pair.opposite);
    const row = element("div", "pair-row");
    const top = element("div", "pair-top");
    top.append(element("span", "", blockName(pair.block_a) + " ↔ " + blockName(pair.block_b)), element("strong", "", direction + " · " + number(majority) + " of " + number(total)));
    const bar = element("div", "agreement-bar");
    const strong = element("span", "agree"); strong.style.flex = String(majority || 0.0001);
    const weak = element("span", "dissent"); weak.style.flex = String(total - majority || 0.0001);
    bar.append(strong, weak);
    tip(row, direction + " orientation: " + number(pair.same) + " same, " + number(pair.opposite) + " opposite.\n" + pair.variant_pairs.map(function (vp) { return vp.a + " ↔ " + vp.b + ": " + vp.same + " same, " + vp.opposite + " opposite"; }).join("\n"));
    row.append(top, bar);
    return row;
  }

  function renderEvidenceRows() {
    const present = EVIDENCE_SOURCES.filter(function (source) { return source.data; });
    const container = byId("evidence-rows");
    container.replaceChildren();
    if (!present.length) {
      const hint = element("p", "evidence-hint");
      const code = function (value) { return element("code", "", value); };
      if (blocks.length > 1) {
        hint.append(document.createTextNode("Phase is unresolved here. Rerun with "), code("--reads sample.bam"), document.createTextNode(" (reads spanning the calls) or "),
          code("--parents FATHER,MOTHER"), document.createTextNode(" (which parent passed on each allele) to see which arrangement the data support."));
      } else {
        hint.append(document.createTextNode("Phase does not change this result. "), code("--reads"), document.createTextNode(" and "), code("--parents"), document.createTextNode(" add read and inheritance evidence when it does."));
      }
      container.append(hint);
      return;
    }
    present.forEach(function (source) {
      const e = source.data;
      const row = element("article", "evidence-row");
      row.dataset.source = source.key;
      const head = element("div", "evidence-head");
      head.append(element("span", "evidence-name", source.name));
      const linked = (data.linked_blocks || []).some(function (group) { return group.sources.indexOf(source.key) >= 0; });
      let outcome, tone;
      if (e.supported.length) { outcome = VERDICT_PHRASES[e.verdict_with]; tone = VERDICT_TONES[e.verdict_with]; }
      else if (linked) { outcome = "Used to link blocks before analysis"; tone = "good"; }
      else if (!e.pairs.length) { outcome = blocks.length > 1 ? source.none : "Nothing left to resolve"; tone = "muted"; }
      else { outcome = "Evidence too weak or mixed"; tone = "muted"; }
      head.append(element("span", "evidence-result tone-" + tone, outcome));
      const withinTotal = e.within_blocks.agree + e.within_blocks.disagree;
      const meta = element("p", "evidence-meta", source.meta(e) + (withinTotal ? " · agrees with the phase shown in " + number(e.within_blocks.agree) + " of " + number(withinTotal) : ""));
      const actions = element("div", "evidence-actions");
      const apply = element("button", "btn btn-small evidence-apply");
      apply.type = "button";
      apply.disabled = !e.supported.length;
      apply.addEventListener("click", function () {
        if (sourceApplied(source.key)) removeSource(source.key);
        else e.supported.forEach(function (c) { applyConstraint({ a: c.a, b: c.b, parity: c.parity, sources: [source.key] }); });
        update();
      });
      const details = element("details", "evidence-details");
      details.append(element("summary", "", e.pairs.length ? "Show " + plural(e.pairs.length, "block pair") : "Details"));
      const body = element("div", "evidence-body");
      if (!e.pairs.length) body.append(element("p", "help", "No two blocks are both covered by clear evidence."));
      e.pairs.forEach(function (pair) { body.append(pairRow(pair)); });
      const within = e.within_blocks.agree + e.within_blocks.disagree;
      if (within) body.append(element("p", "help", source.within(e.within_blocks.agree, within)));
      body.append(element("p", "help", source.rule(e)));
      details.append(body);
      actions.append(apply, details);
      row.append(head, meta, actions);
      container.append(row);
    });
    const agreement = byId("evidence-agreement");
    if (present.length === 2 && present.every(function (source) { return source.data.supported.length; })) {
      const parity = new Map(present[0].data.supported.map(function (c) { return [c.a + "|" + c.b, c.parity]; }));
      const shared = present[1].data.supported.filter(function (c) { return parity.has(c.a + "|" + c.b); });
      const conflicts = shared.filter(function (c) { return parity.get(c.a + "|" + c.b) !== c.parity; });
      agreement.hidden = !shared.length;
      agreement.className = "evidence-agreement" + (conflicts.length ? " disagree" : "");
      agreement.textContent = conflicts.length
        ? "Reads and inheritance disagree on " + plural(conflicts.length, "relationship") + "; check the calls before relying on either."
        : "Reads and inheritance agree on all " + plural(shared.length, "relationship") + " they share.";
    }
  }

  function sourceApplied(key) {
    const source = EVIDENCE_SOURCES.find(function (s) { return s.key === key; });
    return Boolean(source && source.data && source.data.supported.length) && source.data.supported.every(function (c) {
      return state.constraints.some(function (k) { return k.sources.indexOf(key) >= 0 && k.parity === c.parity && ((k.a === c.a && k.b === c.b) || (k.a === c.b && k.b === c.a)); });
    });
  }

  function refreshEvidenceRows() {
    document.querySelectorAll(".evidence-row").forEach(function (row) {
      const applied = sourceApplied(row.dataset.source);
      const button = row.querySelector(".evidence-apply");
      button.textContent = applied ? "Applied · undo" : "Apply";
      button.classList.toggle("is-applied", applied);
      row.classList.toggle("applied", applied);
    });
  }

  /* ---------- Phase blocks ---------- */

  function renderBlocks(selected) {
    const strip = byId("block-strip");
    strip.replaceChildren();
    blocks.forEach(function (block, index) {
      const card = element("div", "block-card");
      const name = element("div", "block-name");
      name.append(document.createTextNode(blockName(block.id)));
      name.append(element("small", "", block.label.replace(/^Unphased .*/, "unphased").replace(/^Phase set /, "PS ").replace(/^Phased calls without PS.*/, "implicit phase set")
        .replace(/^Linked by (.*) \((\d+) blocks\)$/, "linked · $1")));
      tip(name, block.label + (block.linked_from ? "\nMerged: " + block.linked_from.join(", ") : "") + "\n" + block.variant_ids.join("\n"));
      const flip = element("button", "flip");
      flip.type = "button";
      flip.append(icon(["M3 5h8l-2.5-2.5", "M13 11H5l2.5 2.5"], "0 0 16 16"), document.createTextNode("Flip"));
      const target = flipTarget(scenarios, blockIds, selected, block.id);
      const allowed = target && state.compatible.indexOf(target) >= 0;
      flip.disabled = !allowed;
      flip.setAttribute("aria-label", "Flip " + blockName(block.id) + " to the other copy");
      if (!target) tip(flip, "With one block, flipping only swaps the copy labels.");
      else if (!allowed) tip(flip, "That arrangement is excluded by your hypotheses.");
      else tip(flip, "Move this block's alleles to the other copy" + (index === 0 ? " (equivalent to flipping every other block)" : ""));
      flip.addEventListener("click", function () { if (target) { state.selectedId = target.id; renderSelected(); } });
      card.append(name, flip);
      const sides = element("div", "block-sides");
      [0, 1].forEach(function (copyIndex) {
        const side = element("div", "block-side");
        const dot = element("span", "side-dot", copyIndex ? "B" : "A");
        dot.style.setProperty("--side-color", copyIndex ? "var(--copy-b)" : "var(--copy-a)");
        const labels = block.variant_ids.map(function (id) { return variantLookup.get(id); }).filter(Boolean).map(function (variant) {
          const allele = alleleOnCopy(selected, copyIndex, variant);
          return allele.toUpperCase() === variant.cds_ref.toUpperCase() ? null : changeLabel(variant, allele);
        }).filter(Boolean);
        const list = element("span", "side-variants" + (labels.length ? "" : " empty"), labels.length ? labels.join(", ") : "reference alleles");
        if (labels.length) tip(list, labels.join("\n"));
        side.append(dot, list);
        sides.append(side);
      });
      card.append(sides);
      strip.append(card);
    });
    const fixed = variants.filter(function (variant) { return !variant.block; });
    if (fixed.length) {
      const card = element("div", "block-card");
      const name = element("div", "block-name");
      name.append(document.createTextNode("Fixed"), element("small", "", "on both copies"));
      card.append(name, element("span", "block-fixed", plural(fixed.length, "variant")));
      const list = element("div", "block-sides");
      const labels = fixed.map(function (variant) { return changeLabel(variant, variant.cds_alleles[0]); });
      const row = element("div", "block-side");
      const variantsText = element("span", "side-variants", labels.join(", "));
      tip(variantsText, "Homozygous calls sit on both copies in every arrangement.\n" + labels.join("\n"));
      row.append(variantsText);
      list.append(row);
      card.append(list);
      strip.append(card);
    }
    if (!strip.children.length) strip.append(element("p", "no-edits", "No heterozygous phase blocks. Both copies are determined by the included calls."));
  }

  /* ---------- Copy tracks ---------- */

  function copyVerdict(haplotype) {
    const consequences = haplotype.consequences || [];
    if (isReferenceCopy(haplotype, reference)) return { tone: "is-reference", label: "Reference-identical CDS", icon: ICONS.check };
    const disruptive = DISRUPTIVE.filter(function (key) { return consequences.indexOf(key) >= 0; });
    if (disruptive.length) return { tone: "is-disrupted", label: consequenceLabel(disruptive[0]), icon: ICONS.alert };
    if (consequences.indexOf("synonymous") >= 0 || consequences.indexOf("protein_sequence_unchanged") >= 0) return { tone: "", label: "Protein unchanged", icon: ICONS.dot };
    return { tone: "is-altered", label: "Altered protein", icon: ICONS.dot };
  }

  function renderCopy(scenario, copyIndex) {
    const haplotype = scenario.haplotypes[copyIndex];
    const other = scenario.haplotypes[1 - copyIndex];
    const copy = element("article", "copy " + (copyIndex ? "copy-b" : "copy-a"));
    const head = element("div", "copy-head");
    const title = element("h3", "", "Copy " + (copyIndex ? "B" : "A"));
    const verdict = copyVerdict(haplotype);
    const pill = element("span", "copy-verdict " + verdict.tone);
    pill.append(icon(verdict.icon), document.createTextNode(verdict.label));
    const headline = element("span", "copy-headline");
    headline.append(element("strong", "", number(haplotype.protein_length) + " aa"), document.createTextNode(" of " + number(data.reference_protein_length) + " in the reference"));
    head.append(element("span", "copy-badge", copyIndex ? "B" : "A"), title, pill, headline);
    copy.append(head);
    const changes = (haplotype.protein_changes || []).filter(function (change) { return change !== "p.="; });
    if (changes.length || (haplotype.protein_changes || []).length) {
      const list = element("div", "protein-changes");
      list.append(element("span", "chips-caption", "Protein"));
      (changes.length ? changes : ["p.= (unchanged)"]).slice(0, 8).forEach(function (change) { list.append(element("span", "chip", change)); });
      if (changes.length > 8) list.append(tip(element("span", "chips-caption", "+" + (changes.length - 8) + " more"), changes.slice(8).join("\n")));
      copy.append(list);
    }

    const label = element("div", "track-label");
    label.append(element("span", "", "Coding sequence"), element("span", "", number(reference.length) + " nt reference · " + plural(exons.length, "exon")));
    copy.append(label);
    const track = element("div", "cds-track");
    const pins = element("div", "pin-layer");
    const bar = element("div", "exon-bar");
    bar.setAttribute("role", "img");
    exons.slice(1).forEach(function (exon) {
      const gap = element("span", "exon-gap");
      gap.style.left = pct(exon.start) + "%";
      bar.append(gap);
    });
    const ids = haplotype.variant_ids || [];
    const described = [];
    ids.forEach(function (id) {
      const variant = variantLookup.get(id);
      if (!variant) return;
      const allele = alleleOnCopy(scenario, copyIndex, variant);
      const name = changeLabel(variant, allele);
      described.push(name);
      const center = pct((variant.cds_start + Math.max(variant.cds_end, variant.cds_start + 1)) / 2);
      const pin = element("span", "pin");
      pin.style.left = center + "%";
      const labelNode = tip(element("span", "pin-label" + (spliceIds.has(variant.id) ? " splice" : ""), name), name + "\n" + variant.id + "\nGT " + variant.gt.join(variant.phased ? "|" : "/") + (variant.block ? " · " + blockName(variant.block) : " · on both copies") + (spliceIds.has(variant.id) ? "\nSplice region: effect on splicing not modelled" : "") + (infoText(variant) ? "\n" + infoText(variant) : "") + readText(variant) + originText(variant));
      pin.append(labelNode, element("span", "pin-stem"));
      pins.append(pin);
      const head = element("span", "pin-head");
      head.style.left = center + "%";
      bar.append(head);
    });
    (other.variant_ids || []).filter(function (id) { return ids.indexOf(id) < 0; }).forEach(function (id) {
      const variant = variantLookup.get(id);
      if (!variant) return;
      const ghost = tip(element("span", "ghost"), changeLabel(variant, alleleOnCopy(scenario, 1 - copyIndex, variant)) + " is on the other copy; this copy has the reference allele.");
      ghost.style.left = pct((variant.cds_start + Math.max(variant.cds_end, variant.cds_start + 1)) / 2) + "%";
      bar.append(ghost);
    });
    bar.setAttribute("aria-label", "Copy " + (copyIndex ? "B" : "A") + " coding sequence" + (described.length ? " with " + described.join(", ") : " without alternate alleles"));
    track.append(pins, bar);
    const axis = element("div", "exon-axis");
    axis.append(element("span", "start", "c.1"), element("span", "end", "c." + number(reference.length)));
    track.append(axis);
    if (described.length) {
      // Shown by layoutPins only when some labels do not fit above the track.
      const chips = element("div", "change-chips");
      chips.hidden = true;
      chips.append(element("span", "chips-caption", plural(described.length, "change") + " on this copy:"));
      described.forEach(function (name) { chips.append(element("span", "chip", name)); });
      track.append(chips);
    }
    copy.append(track);
    if (!ids.length) copy.insertBefore(element("p", "no-edits", "No alternate coding alleles on this copy."), track);

    const proteinLabel = element("div", "track-label");
    proteinLabel.append(element("span", "", "Predicted protein"), element("span", "", haplotype.first_stop == null ? "no stop codon within the CDS" : "first stop at codon " + number(haplotype.first_stop)));
    const protein = element("div", "protein-track");
    protein.setAttribute("role", "img");
    protein.setAttribute("aria-label", "Predicted translation of " + haplotype.protein_length + " amino acids, reference " + data.reference_protein_length);
    const fillWidth = 100 * haplotype.protein_length / proteinScale;
    const rail = element("span", "protein-rail");
    const fill = tip(element("span", "protein-fill"), number(haplotype.protein_length) + " amino acids before the first stop");
    fill.style.width = Math.max(0, Math.min(100, fillWidth)) + "%";
    const refTick = tip(element("span", "ref-tick"), "Reference protein: " + number(data.reference_protein_length) + " aa");
    refTick.style.left = Math.min(100, 100 * data.reference_protein_length / proteinScale) + "%";
    protein.append(rail, fill, refTick);
    if (haplotype.first_stop != null) {
      const stop = tip(element("span", "stop-mark"), "First stop at codon " + number(haplotype.first_stop));
      stop.style.left = Math.max(1.2, Math.min(98.8, fillWidth)) + "%";
      protein.append(stop);
    }
    const proteinAxis = element("div", "protein-axis");
    proteinAxis.append(element("span", "", "0"), element("span", "", number(proteinScale) + " aa"));
    copy.append(proteinLabel, protein, proteinAxis);

    const pills = element("ul", "pills");
    (haplotype.consequences || []).forEach(function (key) { pills.append(element("li", "pill sev-" + consequenceTone(key), consequenceLabel(key))); });
    if (pills.children.length) copy.append(pills);

    const sequences = element("details", "sequence-details");
    sequences.append(element("summary", "", "Sequences"));
    sequences.addEventListener("toggle", function () {
      if (!sequences.open || sequences.querySelector(".sequence-content")) return;
      const sequence = sequencesFor(haplotype);
      const content = element("div", "sequence-content");
      content.append(element("h4", "", "Coding sequence · transcript orientation · " + number(sequence.cds.length) + " nt"), element("pre", "", sequence.cds),
        element("h4", "", "Predicted translation · " + number(sequence.protein.length) + " aa"), element("pre", "", sequence.protein || "No amino acids before the first stop"));
      sequences.append(content);
    });
    copy.append(sequences);
    return copy;
  }

  /** Lay out pin labels: no label overlaps another, and no stem crosses a label. */
  function layoutPins() {
    document.querySelectorAll(".pin-layer").forEach(function (layer) {
      const width = layer.clientWidth;
      const pins = Array.from(layer.children).map(function (pin) {
        pin.classList.remove("collapsed");
        return { pin: pin, label: pin.firstChild, center: parseFloat(pin.style.left) / 100 * width, width: pin.firstChild.offsetWidth };
      }).sort(function (a, b) { return a.center - b.center; });
      const spots = placeLabels(pins, width, MAX_LABEL_ROWS);
      let rows = 0;
      pins.forEach(function (item, index) {
        const spot = spots[index];
        if (spot.row < 0) { item.pin.classList.add("collapsed"); item.pin.style.height = "0px"; return; }
        rows = Math.max(rows, spot.row + 1);
        item.label.style.left = (spot.left - item.center) + "px";
        item.pin.style.height = (30 + spot.row * 26) + "px";
      });
      layer.style.height = rows ? (30 + (rows - 1) * 26) + "px" : "0px";
      const chips = layer.parentElement.querySelector(".change-chips");
      if (chips) chips.hidden = !spots.some(function (spot) { return spot.row < 0; });
    });
  }

  /* ---------- Outcomes ---------- */

  function renderOutcomes(groups, selected) {
    const grid = byId("outcome-options");
    grid.replaceChildren();
    const visible = state.showAllOutcomes ? groups : groups.slice(0, OUTCOME_PREVIEW);
    visible.forEach(function (group) {
      const first = group[0];
      const card = element("button", "outcome-card");
      card.type = "button";
      card.dataset.outcomeKey = proteinKey(first);
      card.setAttribute("aria-pressed", String(proteinKey(first) === proteinKey(selected)));
      const top = element("div", "outcome-top");
      const coding = new Set(group.map(function (scenario) { return scenario.outcome_key; })).size;
      top.append(element("strong", "", "Outcome " + outcomeNumber.get(proteinKey(first))), element("span", "", plural(group.length, "arrangement")));
      if (coding > 1) tip(top, plural(coding, "coding-sequence variant") + " of this protein pair (synonymous differences).");
      const glyph = element("div", "glyph");
      first.haplotypes.forEach(function (haplotype, copyIndex) {
        const row = element("div", "glyph-row");
        row.style.setProperty("--g-color", copyIndex ? "var(--copy-b)" : "var(--copy-a)");
        const rail = element("span", "glyph-rail");
        const fill = element("span");
        fill.style.width = Math.min(100, 100 * haplotype.protein_length / proteinScale) + "%";
        const tick = element("i");
        tick.style.left = Math.min(100, 100 * data.reference_protein_length / proteinScale) + "%";
        rail.append(fill, tick);
        row.append(element("b", "", copyIndex ? "B" : "A"), rail, element("span", "", number(haplotype.protein_length) + " aa"));
        glyph.append(row);
      });
      const flags = element("div", "outcome-flags");
      const kept = first.haplotypes.filter(function (haplotype) { return keepsReferenceProtein(haplotype, reference); }).length;
      flags.append(element("span", kept ? "sev-good" : "sev-critical", kept ? plural(kept, "reference protein") : "No reference protein"));
      const seen = new Set();
      first.haplotypes.forEach(function (haplotype) {
        (haplotype.consequences || []).forEach(function (key) {
          if (DISRUPTIVE.indexOf(key) >= 0 && !seen.has(key)) { seen.add(key); flags.append(element("span", "sev-critical", consequenceLabel(key))); }
        });
      });
      const described = first.haplotypes.map(function (haplotype, copyIndex) {
        const changes = (haplotype.protein_changes || []).filter(function (change) { return change !== "p.="; });
        return (copyIndex ? "B: " : "A: ") + (changes.length ? changes.slice(0, 3).join(", ") + (changes.length > 3 ? " …" : "") : "reference protein");
      });
      card.append(top, glyph, flags, element("p", "outcome-changes", described.join("  ·  ")));
      card.setAttribute("aria-label", "Outcome " + outcomeNumber.get(proteinKey(first)) + ", " + plural(group.length, "arrangement") + ", copy proteins " + first.haplotypes.map(function (haplotype) { return haplotype.protein_length + " amino acids"; }).join(" and "));
      card.addEventListener("click", function () { state.selectedId = first.id; renderSelected(); });
      grid.append(card);
    });
    const more = byId("show-all-outcomes");
    more.hidden = groups.length <= OUTCOME_PREVIEW;
    more.textContent = state.showAllOutcomes ? "Show fewer outcomes" : "Show all " + number(groups.length) + " outcomes";
  }

  /* ---------- Evidence ---------- */

  function renderEvidence() {
    const summary = byId("evidence-summary");
    const stack = element("div", "stack");
    stack.setAttribute("role", "img");
    const legend = element("div", "stack-legend");
    const names = { reference: "Reference call", variant: "Variant call", unknown: "Unknown" };
    stack.setAttribute("aria-label", ["reference", "variant", "unknown"].map(function (status) { return number(evidence.counts[status]) + " " + names[status].toLowerCase(); }).join(", "));
    ["reference", "variant", "unknown"].forEach(function (status) {
      const count = evidence.counts[status];
      const share = evidence.total ? 100 * count / evidence.total : 0;
      if (count > 0) {
        const part = tip(element("span", "ev-" + status), names[status] + ": " + number(count) + " bases (" + share.toFixed(1) + "%)");
        part.style.flex = String(count);
        stack.append(part);
      }
      const item = element("div");
      item.append(element("i", "ev-" + status), element("strong", "", number(count)), element("span", "", names[status].toLowerCase() + " bases"), element("small", "", share < 1 && count ? "<1%" : Math.round(share) + "%"));
      legend.append(item);
    });
    summary.append(stack, legend);

    const map = byId("evidence-map");
    const head = element("div", "map-head");
    head.append(element("span", "", "By coding position"), element("span", "", "c.1 → c." + number(reference.length) + " · " + (transcript.strand === "-" ? "reverse strand, shown 5′→3′" : "forward strand")));
    const bar = element("div", "map-bar");
    bar.setAttribute("role", "img");
    bar.setAttribute("aria-label", "Evidence along the coding sequence; see the counts above.");
    evidence.pieces.forEach(function (piece) {
      const seg = element("span", "seg ev-" + piece.status);
      seg.style.left = pct(piece.start) + "%";
      seg.style.width = Math.max(0.25, pct(piece.end) - pct(piece.start)) + "%";
      tip(seg, "c." + number(piece.start + 1) + "–" + number(piece.end) + "  ·  " + transcript.contig + ":" + number(piece.genomic[0] + 1) + "–" + number(piece.genomic[1]) + "\n" + names[piece.status] + (piece.detail ? "\n" + piece.detail : ""));
      bar.append(seg);
    });
    exons.slice(1).forEach(function (exon) {
      const gap = element("span", "exon-gap");
      gap.style.left = pct(exon.start) + "%";
      bar.append(gap);
    });
    const labels = element("div", "exon-labels");
    exons.forEach(function (exon) {
      const width = pct(exon.end) - pct(exon.start);
      if (exons.length > 1 && width < 3) return;
      const label = tip(element("span", "", exons.length > 12 || width < 8 ? String(exon.index + 1) : "Exon " + (exon.index + 1)), "Coding exon " + (exon.index + 1) + "\nc." + number(exon.start + 1) + "–" + number(exon.end) + "\n" + transcript.contig + ":" + number(exon.left + 1) + "–" + number(exon.right));
      label.style.left = pct((exon.start + exon.end) / 2) + "%";
      labels.append(label);
    });
    map.append(head, bar, labels);
  }

  function renderSplice() {
    if (!spliceRegion.length) return;
    byId("splice-panel").hidden = false;
    text("splice-count", plural(spliceRegion.length, "splice-region variant"));
    const tbody = byId("splice-body");
    spliceRegion.forEach(function (hit) {
      const row = element("tr");
      const side = hit.where === "intron" ? "intronic" : hit.where === "exon" ? "exonic" : "spans junction";
      [[hit.id, "mono"], [hit.gt.join(hit.phased ? "|" : "/"), "mono"], [side, ""], [hit.distance + " bp from " + transcript.contig + ":" + number(hit.junction), ""],
        [spliceIds.has(hit.id) && variants.some(function (v) { return v.id === hit.id; }) ? "Yes, coding change modelled" : "No, intronic", ""]].forEach(function (cell) {
        row.append(element("td", cell[1], cell[0]));
      });
      tbody.append(row);
    });
  }

  function renderVariants() {
    text("variant-count", plural(variants.length, "called coding variant"));
    byId("info-head").hidden = !hasInfo;
    const tbody = byId("variants-body");
    if (!variants.length) {
      const row = element("tr");
      const cell = element("td", "", "No called coding variants in the supported input evidence."); cell.colSpan = 10; row.append(cell); tbody.append(row);
    }
    variants.forEach(function (variant) {
      const row = element("tr");
      const changes = Array.from(new Set(variant.cds_alleles.filter(function (allele) { return allele.toUpperCase() !== variant.cds_ref.toUpperCase(); }).map(function (allele) { return changeLabel(variant, allele); })));
      const block = variant.block || variant.block_id;
      const alone = (variant.effects || []).map(function (effect) { return effect.protein_change; }).join(", ");
      [[variant.id + (spliceIds.has(variant.id) ? " · splice region" : ""), "mono"], [changes.join(", "), "mono"], [alone || "—", "mono"], [transcript.contig + ":" + number(variant.pos + 1), "mono"], [variant.ref + " → " + Array.from(new Set(variant.alleles)).join(" / "), "mono"],
        [variant.gt.join(variant.phased ? "|" : "/"), "mono"], [variant.phase_set == null ? "—" : variant.phase_set, "mono"], [block ? blockName(block) : "Both copies", ""],
        [variant.gq == null ? "—" : number(variant.gq), "num"], [variant.dp == null ? "—" : number(variant.dp), "num"]].forEach(function (cell) {
        row.append(element("td", cell[1], cell[0]));
      });
      if (hasInfo) row.append(element("td", "info-cell", infoText(variant) || "—"));
      tbody.append(row);
    });
  }

  /* ---------- Hypotheses ---------- */

  function renderConstraints() {
    const list = byId("active-constraints");
    list.replaceChildren();
    text("constraint-count", state.constraints.length);
    byId("reset-constraints").disabled = !state.constraints.length;
    byId("hypothesis-badge").hidden = !state.constraints.length;
    if (!state.constraints.length) list.append(element("li", "no-constraints", "None yet. Every supported arrangement is shown."));
    state.constraints.forEach(function (constraint, index) {
      const row = element("li");
      const label = element("span", "", blockName(constraint.a) + " ↔ " + blockName(constraint.b));
      const evidenceSources = constraint.sources.filter(function (source) { return source !== "hypothesis"; });
      label.append(element("span", "relation" + (evidenceSources.length ? " from-reads" : ""), (constraint.parity === 0 ? "same" : "opposite") + (evidenceSources.length ? " · " + evidenceSources.join(", ") : "")));
      const remove = element("button", "remove-constraint", "×");
      remove.type = "button";
      remove.setAttribute("aria-label", "Remove hypothesis between " + blockName(constraint.a) + " and " + blockName(constraint.b));
      remove.addEventListener("click", function () { state.constraints.splice(index, 1); update(); });
      row.append(label, remove);
      list.append(row);
    });
  }

  function applyConstraint(constraint) { state.constraints = mergeConstraint(state.constraints, constraint); }

  function removeSource(source) { state.constraints = withdrawSource(state.constraints, source); }



  function renderInformative() {
    const pairs = (data.informative_pairs || []).filter(function (pair) { return pair.score > 0 && blockLookup.has(pair.block_a) && blockLookup.has(pair.block_b); });
    if (!pairs.length) return;
    byId("informative-section").hidden = false;
    const best = pairs[0].score;
    pairs.slice(0, 5).forEach(function (pair) {
      const button = element("button", "informative-button");
      button.type = "button";
      tip(button, pair.description || "Select these blocks for a hypothetical orientation relationship");
      const bar = element("span", "bar");
      const fill = element("span");
      fill.style.width = (100 * pair.score / best) + "%";
      bar.append(fill);
      button.append(element("span", "", blockName(pair.block_a) + " ↔ " + blockName(pair.block_b)), element("span", "score", plural(pair.score, "outcome pair")), bar);
      button.addEventListener("click", function () {
        byId("block-a").value = String(pair.block_a); byId("block-b").value = String(pair.block_b); validatePair(); byId("relationship-same").focus();
      });
      byId("informative-pairs").append(button);
    });
  }

  /* ---------- Selection ---------- */

  function renderSelected() {
    const selected = state.compatible.find(function (scenario) { return scenario.id === state.selectedId; });
    if (!selected) return;
    byId("scenario-select").value = String(selected.id);
    const index = state.compatible.indexOf(selected);
    byId("previous-scenario").disabled = index === 0;
    byId("next-scenario").disabled = index === state.compatible.length - 1;
    renderBlocks(selected);
    const tracks = byId("copy-tracks");
    tracks.replaceChildren(renderCopy(selected, 0), renderCopy(selected, 1));
    layoutPins();
    byId("outcome-options").querySelectorAll(".outcome-card").forEach(function (card) { card.setAttribute("aria-pressed", String(card.dataset.outcomeKey === proteinKey(selected))); });
  }

  function update() {
    renderConstraints();
    state.compatible = filterScenarios(scenarios, state.constraints);
    const groups = groupByProtein(state.compatible);
    text("compatible-count", number(state.compatible.length));
    text("outcome-count", number(groups.length));
    renderVerdict();
    refreshEvidenceRows();
    text("result-explanation", state.compatible.length
      ? "Arrangements that give the same pair of proteins are grouped; synonymous differences are noted on hover. Select one to view it above."
      : "The applied relationships exclude every supported arrangement.");
    byId("empty-state").hidden = state.compatible.length > 0;
    byId("results-content").hidden = !state.compatible.length;
    byId("outcome-options").hidden = !state.compatible.length;
    const select = byId("scenario-select");
    select.replaceChildren();
    state.compatible.forEach(function (scenario) {
      select.append(option(scenario.id, "Arrangement " + number(scenarioNumber.get(scenario.id)) + " · outcome " + outcomeNumber.get(proteinKey(scenario))));
    });
    select.disabled = state.compatible.length < 2;
    if (!state.compatible.some(function (scenario) { return scenario.id === state.selectedId; })) state.selectedId = state.compatible.length ? state.compatible[0].id : null;
    const selected = state.compatible.find(function (scenario) { return scenario.id === state.selectedId; });
    if (selected) renderOutcomes(groups, selected);
    else byId("show-all-outcomes").hidden = true;
    renderSelected();
  }

  function reset() { state.constraints = []; update(); }

  function validatePair() { byId("add-constraint").disabled = byId("block-a").value === byId("block-b").value; }

  /* ---------- Tooltip ---------- */

  function setupTooltip() {
    const tooltip = byId("tooltip");
    let active = null;
    function show(target, x, y) {
      active = target;
      tooltip.textContent = target.dataset.tip;
      tooltip.hidden = false;
      const box = tooltip.getBoundingClientRect();
      const left = Math.max(8, Math.min(window.innerWidth - box.width - 8, x - box.width / 2));
      const top = y - box.height - 12 < 8 ? y + 18 : y - box.height - 12;
      tooltip.style.left = left + "px";
      tooltip.style.top = top + "px";
    }
    function hide() { active = null; tooltip.hidden = true; }
    document.addEventListener("pointermove", function (event) {
      const target = event.target.closest ? event.target.closest("[data-tip]") : null;
      if (!target) { if (active) hide(); return; }
      show(target, event.clientX, event.clientY);
    });
    document.addEventListener("focusin", function (event) {
      const target = event.target.closest ? event.target.closest("[data-tip]") : null;
      if (!target) return;
      const box = target.getBoundingClientRect();
      show(target, box.left + box.width / 2, box.top);
    });
    document.addEventListener("focusout", hide);
    document.addEventListener("scroll", hide, { passive: true });
    document.addEventListener("keydown", function (event) { if (event.key === "Escape") hide(); });
  }

  /* ---------- Wiring ---------- */

  renderIdentity();
  renderEvidence();
  renderVariants();
  renderSplice();
  blocks.forEach(function (block) {
    ["block-a", "block-b"].forEach(function (id) { const item = option(block.id, blockName(block.id)); item.title = block.label; byId(id).append(item); });
  });
  if (blocks.length > 1) byId("block-b").value = String(blocks[1].id);
  else {
    byId("constraint-form").hidden = true;
    text("constraint-help", blocks.length
      ? "There is only one phase block, and it is already internally oriented. Swapping the copy labels does not create a different coding outcome."
      : "No heterozygous phase blocks are available to relate. Both copies follow directly from the included calls.");
  }
  byId("block-a").addEventListener("change", validatePair);
  byId("block-b").addEventListener("change", validatePair);
  byId("constraint-form").addEventListener("submit", function (event) {
    event.preventDefault();
    const a = byId("block-a").value;
    const b = byId("block-b").value;
    if (a === b || !blockLookup.has(a) || !blockLookup.has(b)) return;
    const parity = Number(document.querySelector("input[name='relationship']:checked").value);
    applyConstraint({ a: a, b: b, parity: parity, sources: ["hypothesis"] });
    update();
  });
  byId("reset-constraints").addEventListener("click", reset);
  byId("empty-reset").addEventListener("click", reset);
  byId("show-all-outcomes").addEventListener("click", function () { state.showAllOutcomes = !state.showAllOutcomes; update(); });
  byId("scenario-select").addEventListener("change", function () {
    const selected = state.compatible.find(function (scenario) { return String(scenario.id) === byId("scenario-select").value; });
    if (selected) state.selectedId = selected.id;
    renderSelected();
  });
  function step(delta) {
    const index = state.compatible.findIndex(function (scenario) { return scenario.id === state.selectedId; });
    const next = state.compatible[index + delta];
    if (next) { state.selectedId = next.id; renderSelected(); }
  }
  byId("previous-scenario").addEventListener("click", function () { step(-1); });
  byId("next-scenario").addEventListener("click", function () { step(1); });
  byId("download-json").addEventListener("click", function () {
    const exported = { analysis: expandAnalysis(data), exploration: { hypothetical: true, constraints: state.constraints, selected_scenario_id: state.selectedId } };
    const blob = new Blob([JSON.stringify(exported, null, 2) + "\n"], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const link = element("a");
    link.href = url;
    link.download = "phasespace-" + String(data.sample).replace(/[^a-zA-Z0-9_-]+/g, "_").slice(0, 80) + ".json";
    document.body.append(link); link.click(); link.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
  });
  let frame = 0;
  window.addEventListener("resize", function () { cancelAnimationFrame(frame); frame = requestAnimationFrame(layoutPins); });
  setupTooltip();
  renderInformative();
  renderEvidenceRows();
  validatePair();
  update();
})();
