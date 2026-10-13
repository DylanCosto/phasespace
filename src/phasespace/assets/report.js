/* All metadata and sequence content enters the document as text, never markup. */
(function () {
  "use strict";

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

  // Export the exact filtering implementation for an independent logic test.
  if (typeof module !== "undefined" && module.exports) {
    module.exports = { filterScenarios: filterScenarios, groupScenarios: groupScenarios };
  }
  if (typeof document === "undefined") return;

  const data = JSON.parse(document.getElementById("phasespace-data").textContent);
  const transcript = data.transcript;
  const blocks = data.blocks || [];
  const scenarios = data.scenarios || [];
  const variants = data.variants || [];
  const variantLookup = new Map(variants.map(function (variant) { return [variant.id, variant]; }));
  const blockLookup = new Map(blocks.map(function (block) { return [block.id, block]; }));
  const state = { constraints: [], selectedId: scenarios.length ? scenarios[0].id : null, compatible: [] };
  const byId = function (id) { return document.getElementById(id); };
  const number = function (value) { return Number(value).toLocaleString("en-US"); };
  const text = function (id, value) { byId(id).textContent = String(value); };
  const blockName = function (id) { const index = blocks.findIndex(function (block) { return block.id === id; }); return index >= 0 ? "Block " + (index + 1) : String(id); };

  function element(tag, className, content) {
    const node = document.createElement(tag);
    if (className) node.className = className;
    if (content !== undefined && content !== null) node.textContent = String(content);
    return node;
  }

  function option(value, label) {
    const node = element("option", "", label);
    node.value = String(value);
    return node;
  }

  function describeVariant(variant) {
    const alternate = (variant.alleles || []).filter(function (allele) { return allele !== variant.ref; });
    return "c." + number(variant.cds_start + 1) + " · " + variant.ref + " → " + Array.from(new Set(alternate)).join(" / ");
  }

  function renderEvidence() {
    const intervals = transcript.cds_intervals || [];
    const container = byId("evidence-map");
    if (!intervals.length) {
      container.append(element("p", "small-note", "No coding intervals available."));
      return;
    }
    const start = Math.min.apply(null, intervals.map(function (interval) { return interval[0]; }));
    const end = Math.max.apply(null, intervals.map(function (interval) { return interval[1]; }));
    const width = Math.max(1, end - start);
    const mapX = function (value) { return ((value - start) / width) * 1000; };
    const label = element("div", "map-label");
    label.append(element("span", "", "Coding intervals · genomic coordinates"), element("span", "", (transcript.strand === "-" ? "← Reverse" : "Forward →") + " strand"));
    container.append(label);
    const ns = "http://www.w3.org/2000/svg";
    const svg = document.createElementNS(ns, "svg");
    svg.setAttribute("viewBox", "0 0 1000 52");
    svg.setAttribute("class", "evidence-svg");
    svg.setAttribute("role", "img");
    svg.setAttribute("aria-label", "Coding evidence across " + transcript.contig + ":" + (start + 1) + "–" + end + ". Unknown intervals are hatched, reference calls teal, and variant calls ochre.");
    const defs = document.createElementNS(ns, "defs");
    const pattern = document.createElementNS(ns, "pattern");
    pattern.setAttribute("id", "unknown-hatch");
    pattern.setAttribute("width", "7"); pattern.setAttribute("height", "7"); pattern.setAttribute("patternUnits", "userSpaceOnUse");
    const patternBackground = document.createElementNS(ns, "rect");
    patternBackground.setAttribute("width", "7"); patternBackground.setAttribute("height", "7"); patternBackground.setAttribute("fill", "#e4e9e3");
    const stripe = document.createElementNS(ns, "path");
    stripe.setAttribute("d", "M-1,1 L1,-1 M0,7 L7,0 M6,8 L8,6"); stripe.setAttribute("stroke", "#c7d2ca"); stripe.setAttribute("stroke-width", "2");
    pattern.append(patternBackground, stripe); defs.append(pattern); svg.append(defs);
    const line = document.createElementNS(ns, "line");
    line.setAttribute("x1", "0"); line.setAttribute("x2", "1000"); line.setAttribute("y1", "25"); line.setAttribute("y2", "25"); line.setAttribute("stroke", "#cdd8d0"); line.setAttribute("stroke-width", "2"); svg.append(line);
    function rectangle(a, b, status, detail) {
      const rect = document.createElementNS(ns, "rect");
      const x = Math.max(0, Math.min(1000, mapX(a)));
      rect.setAttribute("x", String(x));
      rect.setAttribute("y", "10");
      rect.setAttribute("width", String(Math.min(1000 - x, Math.max(status === "variant" ? 3 : 1, mapX(b) - x))));
      rect.setAttribute("height", "30");
      rect.setAttribute("rx", "2");
      rect.setAttribute("fill", status === "reference" ? "#13786c" : status === "variant" ? "#b06b26" : "url(#unknown-hatch)");
      const title = document.createElementNS(ns, "title");
      title.textContent = transcript.contig + ":" + number(a + 1) + "–" + number(b) + " · " + status + (detail ? " · " + detail : "");
      rect.append(title); svg.append(rect);
    }
    intervals.forEach(function (interval) { rectangle(interval[0], interval[1], "unknown", "No confident reference or variant call provided"); });
    const counts = { reference: 0, variant: 0, unknown: 0 };
    (data.evidence || []).forEach(function (evidence) {
      const status = evidence.status === "reference" || evidence.status === "variant" ? evidence.status : "unknown";
      intervals.forEach(function (interval) {
        const a = Math.max(interval[0], evidence.start);
        const b = Math.min(interval[1], evidence.end);
        if (b > a) {
          counts[status] += b - a;
          rectangle(a, b, status, evidence.detail);
        }
      });
    });
    container.append(svg);
    const axis = element("div", "axis-labels");
    axis.append(element("span", "", transcript.contig + ":" + number(start + 1)), element("span", "", number(end)));
    container.append(axis);
    const total = intervals.reduce(function (sum, interval) { return sum + interval[1] - interval[0]; }, 0);
    // Missing evidence still counts as unknown, never as reference.
    counts.unknown = Math.max(counts.unknown, total - counts.reference - counts.variant);
    const key = element("div", "evidence-key");
    key.append(element("span", "", number(total) + " coding bases"), element("span", "", number(counts.reference) + " reference-called"), element("span", "", number(counts.variant) + " variant-called"), element("span", "", number(counts.unknown) + " unknown"));
    container.append(key);
  }

  function renderVariants() {
    text("variant-count", number(variants.length) + " called coding variant" + (variants.length === 1 ? "" : "s"));
    const tbody = byId("variants-body");
    if (!variants.length) {
      const row = element("tr");
      const cell = element("td", "", "No called coding variants in the supported input evidence."); cell.colSpan = 8; row.append(cell); tbody.append(row);
    }
    variants.forEach(function (variant) {
      const row = element("tr");
      const block = variant.block || variant.block_id;
      [variant.id, transcript.contig + ":" + number(variant.pos + 1), variant.ref + " → " + Array.from(new Set(variant.alleles)).join(" / "), variant.gt.join(variant.phased ? "|" : "/"), variant.phase_set == null ? "—" : variant.phase_set, block ? blockName(block) : "Fixed", variant.gq == null ? "—" : variant.gq, variant.dp == null ? "—" : variant.dp].forEach(function (value, index) {
        row.append(element("td", index < 6 ? "mono" : "", value));
      });
      tbody.append(row);
    });
  }

  function renderConstraints() {
    const list = byId("active-constraints");
    list.replaceChildren();
    text("constraint-count", state.constraints.length);
    byId("reset-constraints").disabled = !state.constraints.length;
    byId("hypothesis-badge").hidden = !state.constraints.length;
    if (!state.constraints.length) list.append(element("li", "no-constraints", "All supported arrangements are available."));
    state.constraints.forEach(function (constraint, index) {
      const row = element("li");
      const relation = constraint.parity === 0 ? "same orientation" : "opposite orientation";
      row.append(element("span", "", blockName(constraint.a) + " ↔ " + blockName(constraint.b) + " · " + relation));
      const remove = element("button", "remove-constraint", "×");
      remove.type = "button";
      remove.setAttribute("aria-label", "Remove hypothesis between " + blockName(constraint.a) + " and " + blockName(constraint.b));
      remove.addEventListener("click", function () { state.constraints.splice(index, 1); update(); });
      row.append(remove); list.append(row);
    });
  }

  function renderCopy(haplotype, index, maxProteinLength) {
    const copy = element("article", "copy-track");
    const heading = element("div", "copy-heading");
    const label = element("h3", "copy-label");
    label.append(element("span", "copy-dot"), document.createTextNode("Copy " + (index + 1)));
    const isReference = haplotype.cds === transcript.reference_cds;
    heading.append(label, element("span", "copy-status" + (isReference ? " reference-identical" : ""), isReference ? "Reference-identical CDS" : "Altered CDS"));
    copy.append(heading);
    const description = element("div", "track-description");
    description.append(element("span", "", "Coding variant map · reference coordinates"), element("span", "", number(haplotype.cds.length) + " nt"));
    copy.append(description);
    const ribbon = element("div", "ribbon cds-ribbon");
    ribbon.setAttribute("role", "img");
    ribbon.setAttribute("aria-label", "Coding variant positions on copy " + (index + 1));
    const variantLabels = element("div", "variant-labels");
    (haplotype.variant_ids || []).forEach(function (id) {
      const variant = variantLookup.get(id);
      if (!variant) return;
      const marker = element("span", "variant-mark");
      const start = Math.max(0, Math.min(99.7, 100 * variant.cds_start / Math.max(1, transcript.reference_cds.length)));
      marker.style.left = start + "%";
      marker.style.width = Math.min(100 - start, Math.max(.3, 100 * (variant.cds_end - variant.cds_start) / Math.max(1, transcript.reference_cds.length))) + "%";
      marker.title = variant.id + " · " + describeVariant(variant);
      ribbon.append(marker);
      const chip = element("span", "variant-label", id);
      chip.title = describeVariant(variant);
      variantLabels.append(chip);
    });
    if (!(haplotype.variant_ids || []).length) variantLabels.append(element("span", "no-edits", "No listed alternate coding alleles on this copy"));
    copy.append(ribbon, variantLabels);
    const protein = element("div", "protein-track");
    const proteinLabel = element("div", "track-description");
    proteinLabel.append(element("span", "", "Predicted translation"), element("span", "", number(haplotype.protein_length) + " aa" + (haplotype.first_stop == null ? " · no stop in supplied CDS" : " · first stop at " + number(haplotype.first_stop))));
    const proteinRibbon = element("div", "ribbon");
    proteinRibbon.setAttribute("role", "img");
    proteinRibbon.setAttribute("aria-label", "Predicted translation of " + haplotype.protein_length + " amino acids" + (haplotype.first_stop == null ? " without a stop within supplied CDS" : "; first stop at position " + haplotype.first_stop));
    const fill = element("span", "protein-fill");
    const lengthPercent = Math.max(0, Math.min(100, 100 * haplotype.protein_length / maxProteinLength));
    fill.style.width = lengthPercent + "%";
    proteinRibbon.append(fill);
    if (haplotype.first_stop != null) {
      const stop = element("span", "stop-mark");
      stop.style.left = Math.min(99.8, lengthPercent) + "%";
      stop.title = "First stop at amino-acid position " + haplotype.first_stop;
      proteinRibbon.append(stop);
    }
    const axis = element("div", "track-axis");
    axis.append(element("span", "", "0"), element("span", "", number(maxProteinLength) + " aa · shared scale"));
    protein.append(proteinLabel, proteinRibbon, axis); copy.append(protein);
    const consequences = element("ul", "consequences");
    (haplotype.consequences || []).forEach(function (consequence) { consequences.append(element("li", "", String(consequence).replace(/_/g, " "))); });
    copy.append(consequences);
    const sequences = element("details", "sequence-details");
    sequences.append(element("summary", "", "Inspect coding & protein sequences"));
    const sequenceContent = element("div", "sequence-content");
    sequenceContent.append(element("h4", "", "CODING SEQUENCE · TRANSCRIPT ORIENTATION"), element("pre", "", haplotype.cds), element("h4", "", "PREDICTED TRANSLATION"), element("pre", "", haplotype.protein || "No translated amino acids before the first stop"));
    sequences.append(sequenceContent); copy.append(sequences);
    return copy;
  }

  function renderSelected() {
    const selected = state.compatible.find(function (scenario) { return scenario.id === state.selectedId; });
    if (!selected) return;
    byId("scenario-select").value = String(selected.id);
    const index = state.compatible.indexOf(selected);
    byId("previous-scenario").disabled = index === 0;
    byId("next-scenario").disabled = index === state.compatible.length - 1;
    const orientations = byId("orientation-row");
    orientations.replaceChildren();
    blocks.forEach(function (block) { const pill = element("span", "orientation-pill", blockName(block.id) + " · orientation " + selected.orientations[block.id]); pill.title = block.label; orientations.append(pill); });
    if (!blocks.length) orientations.append(element("span", "orientation-pill", "No heterozygous phase blocks"));
    const maxLength = Math.max.apply(null, [1, data.reference_protein_length || 0].concat(selected.haplotypes.map(function (haplotype) { return haplotype.protein_length; })));
    const tracks = byId("copy-tracks");
    tracks.replaceChildren();
    selected.haplotypes.forEach(function (haplotype, copyIndex) { tracks.append(renderCopy(haplotype, copyIndex, maxLength)); });
    byId("outcome-options").querySelectorAll("button").forEach(function (button) { button.setAttribute("aria-pressed", String(button.dataset.outcomeKey === selected.outcome_key)); });
  }

  function update() {
    renderConstraints();
    state.compatible = filterScenarios(scenarios, state.constraints);
    const groups = groupScenarios(state.compatible);
    text("compatible-count", number(state.compatible.length));
    text("outcome-count", number(groups.length));
    text("outcome-summary-count", number(groups.length) + " outcome" + (groups.length === 1 ? "" : "s"));
    const referenceCounts = state.compatible.map(function (scenario) { return scenario.haplotypes.filter(function (haplotype) { return haplotype.cds === transcript.reference_cds; }).length; });
    const min = Math.min.apply(null, referenceCounts);
    const max = Math.max.apply(null, referenceCounts);
    const proteinCount = new Set(state.compatible.map(function (scenario) { return scenario.protein_outcome_key; })).size;
    const proteinExplanation = state.compatible.length && state.compatible.every(function (scenario) { return Boolean(scenario.protein_outcome_key); }) ? proteinCount + " predicted protein outcome" + (proteinCount === 1 ? ". " : "s. ") : "";
    text("result-explanation", state.compatible.length ? proteinExplanation + (min === max ? String(min) : min + "–" + max) + " reference-identical coding " + (min === max && min === 1 ? "copy" : "copies") + " across the compatible arrangements. Unknown evidence still limits interpretation." : "Conflicting orientation hypotheses exclude every supported arrangement.");
    byId("empty-state").hidden = state.compatible.length > 0;
    byId("results-content").hidden = !state.compatible.length;
    const select = byId("scenario-select");
    select.replaceChildren();
    state.compatible.forEach(function (scenario, index) {
      select.append(option(scenario.id, "Arrangement " + (index + 1) + " of " + state.compatible.length + " · " + scenario.id));
    });
    if (!state.compatible.some(function (scenario) { return scenario.id === state.selectedId; })) state.selectedId = state.compatible.length ? state.compatible[0].id : null;
    const options = byId("outcome-options");
    options.replaceChildren();
    groups.forEach(function (group, index) {
      const first = group[0];
      const button = element("button", "outcome-button");
      button.type = "button";
      button.dataset.outcomeKey = first.outcome_key;
      button.setAttribute("aria-pressed", "false");
      button.append(element("strong", "", "Outcome " + (index + 1)), element("span", "", first.haplotypes.map(function (haplotype) { return number(haplotype.protein_length) + " aa"; }).join(" + ")), element("span", "", group.length + " compatible arrangement" + (group.length === 1 ? "" : "s")));
      button.addEventListener("click", function () { state.selectedId = first.id; renderSelected(); });
      options.append(button);
    });
    renderSelected();
  }

  function reset() { state.constraints = []; update(); }

  text("gene-title", transcript.gene || transcript.id);
  text("locus-meta", "Sample " + data.sample + "  /  " + transcript.id + "  /  " + transcript.contig + "  /  " + (transcript.strand === "-" ? "reverse" : "forward") + " strand");
  text("source-meta", data.source ? "Source: " + data.source : "Source metadata unavailable");
  text("total-arrangements", number(scenarios.length));
  text("total-outcomes", number(groupScenarios(scenarios).length));
  byId("synthetic-label").hidden = data.synthetic !== true;
  if (Number.isInteger(data.protein_outcome_count)) {
    byId("protein-summary").hidden = false;
    text("total-protein-outcomes", number(data.protein_outcome_count));
  }
  document.title = "PhaseSpace · " + (transcript.gene || transcript.id) + " · " + data.sample;
  if ((data.warnings || []).length) {
    byId("warning-panel").hidden = false;
    data.warnings.forEach(function (warning) { byId("warnings").append(element("li", "", warning)); });
  }
  renderEvidence();
  renderVariants();
  blocks.forEach(function (block) {
    const tag = element("span", "block-tag");
    tag.append(element("strong", "", blockName(block.id)), document.createTextNode(" · " + block.variant_ids.length + " variant" + (block.variant_ids.length === 1 ? "" : "s")));
    tag.title = block.label + " · " + block.variant_ids.join(", ");
    byId("phase-blocks").append(tag);
    ["block-a", "block-b"].forEach(function (id) { const item = option(block.id, blockName(block.id)); item.title = block.label; byId(id).append(item); });
  });
  if (blocks.length > 1) byId("block-b").value = String(blocks[1].id);
  else {
    byId("constraint-form").hidden = true;
    text("constraint-help", blocks.length ? "One phase block is already internally oriented. Swapping both copy labels does not create a different coding outcome." : "No heterozygous phase blocks are available to constrain. Both copies still depend on the supplied evidence.");
  }
  function validatePair() { byId("add-constraint").disabled = byId("block-a").value === byId("block-b").value; }
  byId("block-a").addEventListener("change", validatePair);
  byId("block-b").addEventListener("change", validatePair);
  byId("constraint-form").addEventListener("submit", function (event) {
    event.preventDefault();
    const a = byId("block-a").value;
    const b = byId("block-b").value;
    if (a === b || !blockLookup.has(a) || !blockLookup.has(b)) return;
    const parity = Number(byId("relationship").value);
    const existing = state.constraints.findIndex(function (constraint) { return (constraint.a === a && constraint.b === b) || (constraint.a === b && constraint.b === a); });
    const constraint = { a: a, b: b, parity: parity };
    if (existing < 0) state.constraints.push(constraint); else state.constraints[existing] = constraint;
    update();
  });
  byId("reset-constraints").addEventListener("click", reset);
  byId("empty-reset").addEventListener("click", reset);
  byId("scenario-select").addEventListener("change", function () {
    const selected = state.compatible.find(function (scenario) { return String(scenario.id) === byId("scenario-select").value; });
    if (selected) state.selectedId = selected.id;
    renderSelected();
  });
  byId("previous-scenario").addEventListener("click", function () {
    const index = state.compatible.findIndex(function (scenario) { return scenario.id === state.selectedId; });
    if (index > 0) state.selectedId = state.compatible[index - 1].id;
    renderSelected();
  });
  byId("next-scenario").addEventListener("click", function () {
    const index = state.compatible.findIndex(function (scenario) { return scenario.id === state.selectedId; });
    if (index < state.compatible.length - 1) state.selectedId = state.compatible[index + 1].id;
    renderSelected();
  });
  const pairs = (data.informative_pairs || []).filter(function (pair) { return pair.score > 0 && blockLookup.has(pair.block_a) && blockLookup.has(pair.block_b); });
  if (pairs.length) {
    byId("informative-section").hidden = false;
    pairs.slice(0, 5).forEach(function (pair) {
      const button = element("button", "informative-button");
      button.type = "button";
      button.title = pair.description || "Select these blocks for a hypothetical orientation relationship";
      button.append(element("span", "", blockName(pair.block_a) + " ↔ " + blockName(pair.block_b)), element("span", "", "Explore →"));
      button.addEventListener("click", function () {
        byId("block-a").value = String(pair.block_a); byId("block-b").value = String(pair.block_b); validatePair(); byId("relationship").focus();
      });
      byId("informative-pairs").append(button);
    });
  }
  byId("download-json").addEventListener("click", function () {
    const exported = { analysis: data, exploration: { hypothetical: true, constraints: state.constraints, selected_scenario_id: state.selectedId } };
    const blob = new Blob([JSON.stringify(exported, null, 2) + "\n"], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const link = element("a");
    link.href = url;
    link.download = "phasespace-" + String(data.sample).replace(/[^a-zA-Z0-9_-]+/g, "_").slice(0, 80) + ".json";
    document.body.append(link); link.click(); link.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
  });
  validatePair();
  update();
})();
