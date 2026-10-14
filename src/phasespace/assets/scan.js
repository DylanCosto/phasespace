/* Scan index: every data-derived string enters the page as text, never markup. */
(function () {
  "use strict";
  const data = JSON.parse(document.getElementById("phasespace-data").textContent);
  const byId = function (id) { return document.getElementById(id); };
  const number = function (value) { return Number(value).toLocaleString("en-US"); };
  const LABELS = { stop_gained: "Premature stop", frameshift: "Frameshift", start_lost: "Start lost", stop_lost: "Stop lost", missense: "Missense",
    inframe_insertion: "In-frame insertion", inframe_deletion: "In-frame deletion", frame_restored: "Frame restored", uncertain_codon: "Uncertain codon" };
  const TIER_NAMES = Object.assign({}, data.tiers, { 5: "Not analysed" });
  const COLUMNS = [
    { key: "tier", label: "Tier" }, { key: "gene", label: "Gene" }, { key: "protein_changes", label: "Protein changes" },
    { key: "heterozygous", label: "Calls", numeric: true }, { key: "arrangements", label: "Arrangements", numeric: true },
    { key: "unknown_cds_bases", label: "No call", numeric: true }
  ];
  const state = { tier: "all", query: "", sort: { key: "tier", dir: 1 }, limit: 300 };
  const rows = data.rows || [];
  const VERDICTS = { good: "Reference protein kept", none: "No reference protein", mixed: "Phase still decides", conflicting: "Reads conflict" };
  if (rows.some(function (row) { return row.verdict_with_reads; })) COLUMNS.splice(3, 0, { key: "verdict_with_reads", label: "With reads" });
  if (rows.some(function (row) { return row.verdict_with_parents; })) COLUMNS.splice(3, 0, { key: "verdict_with_parents", label: "With parents" });

  function element(tag, className, content) {
    const node = document.createElement(tag);
    if (className) node.className = className;
    if (content !== undefined && content !== null) node.textContent = String(content);
    return node;
  }

  function renderHeader() {
    byId("sample-label").textContent = "Sample " + data.sample;
    const scope = number(data.counts.transcripts_scanned) + " transcripts" + (data.counts.genes_on_absent_contigs ? " (" + number(data.counts.genes_on_absent_contigs) + " more on contigs absent from the inputs)" : "");
    [data.source, "Annotation " + data.annotation, scope, data.seconds + " s"].forEach(function (value) {
      byId("scan-meta").append(element("span", "", value));
    });
    byId("source-meta").textContent = "Source: " + data.source + "  ·  PhaseSpace " + data.version;
    byId("version-label").textContent = data.version;
    document.title = "Scan · " + data.sample + " · PhaseSpace";
    const stats = [["With coding variants", data.counts.with_coding_variants, ""], ["Phase decides", data.counts.phase_decides, "tier-1"],
      ["No reference protein, disruptive", data.counts.no_reference_protein_disruptive, "tier-2"], ["No reference protein", data.counts.no_reference_protein, "tier-3"],
      ["Reports written", data.counts.reports, ""]];
    stats.forEach(function (stat) {
      const box = element("div", "stat");
      const label = element("dt", stat[2] ? "stat-tier " + stat[2] : "", stat[0]);
      box.append(label, element("dd", "", number(stat[1])));
      byId("scan-stats").append(box);
    });
    (data.skipped || []).slice(0, 50).forEach(function (item) { byId("skipped").append(element("li", "", item.gene + ": " + item.reason)); });
    byId("skipped-panel").hidden = !(data.skipped || []).length;
  }

  function renderFilters() {
    const counts = { all: rows.length };
    rows.forEach(function (row) { counts[row.tier] = (counts[row.tier] || 0) + 1; });
    const container = byId("tier-filters");
    container.replaceChildren();
    ["all", 1, 2, 3, 4, 5].forEach(function (tier) {
      if (tier !== "all" && !counts[tier]) return;
      const chip = element("button", "tier-chip" + (tier === "all" ? "" : " tier-" + tier));
      chip.type = "button";
      chip.dataset.tier = String(tier);
      chip.setAttribute("aria-pressed", String(String(state.tier) === String(tier)));
      chip.append(document.createTextNode(tier === "all" ? "All" : TIER_NAMES[tier]), element("b", "", number(counts[tier])));
      chip.addEventListener("click", function () { state.tier = String(tier); state.limit = 300; renderFilters(); renderRows(); });
      container.append(chip);
    });
  }

  function renderHead() {
    const head = byId("scan-head");
    head.replaceChildren();
    COLUMNS.forEach(function (column) {
      const th = element("th", column.numeric ? "num" : "");
      th.setAttribute("aria-sort", state.sort.key === column.key ? (state.sort.dir > 0 ? "ascending" : "descending") : "none");
      const button = element("button", "", column.label);
      button.type = "button";
      button.addEventListener("click", function () {
        state.sort = { key: column.key, dir: state.sort.key === column.key ? -state.sort.dir : 1 };
        renderHead(); renderRows();
      });
      th.append(button);
      head.append(th);
    });
  }

  function compare(a, b) {
    const key = state.sort.key;
    const x = a[key], y = b[key];
    let result;
    if (typeof x === "number" || typeof y === "number") result = (x === undefined || x === "" ? Infinity : x) - (y === undefined || y === "" ? Infinity : y);
    else result = String(x || "").localeCompare(String(y || ""));
    if (result === 0 && key === "tier") result = (b.disruptive || "").length - (a.disruptive || "").length || String(a.gene).localeCompare(String(b.gene));
    return result * state.sort.dir;
  }

  function renderRows() {
    const query = state.query.trim().toLowerCase();
    const visible = rows.filter(function (row) {
      return (state.tier === "all" || String(row.tier) === state.tier) && (!query || (row.gene + " " + row.transcript).toLowerCase().indexOf(query) >= 0);
    }).sort(compare);
    const body = byId("scan-body");
    body.replaceChildren();
    visible.slice(0, state.limit).forEach(function (row) {
      const tr = element("tr");
      const cells = {};
      cells.tier = element("td");
      cells.tier.append(element("span", "tier-badge tier-" + row.tier, TIER_NAMES[row.tier]));
      cells.gene = element("td", "gene");
      const name = row.report ? element("a", "", row.gene) : element("span", "gene-name", row.gene);
      if (row.report) name.href = row.report;
      cells.gene.append(name, element("span", "transcript", row.transcript));
      cells.protein_changes = element("td");
      if (row.tier === 5) cells.protein_changes.append(element("span", "status-note", row.status));
      else {
        const pills = element("div", "change-pills");
        const changes = (row.protein_changes || "").split("; ").filter(Boolean);
        changes.slice(0, 4).forEach(function (change) { pills.append(element("span", "pchange", change)); });
        if (changes.length > 4) { const more = element("span", "more", "+" + (changes.length - 4)); more.title = changes.slice(4).join("\n"); pills.append(more); }
        const disruptive = (row.disruptive || "").split(",").filter(Boolean);
        disruptive.forEach(function (key) { pills.append(element("span", "disruptive", LABELS[key] || key)); });
        if (Number(row.splice_region)) pills.append(element("span", "splice", "Splice region ×" + row.splice_region));
        if (!changes.length && !disruptive.length) pills.append(element("span", "", row.changes ? "Synonymous only" : "No protein change"));
        cells.protein_changes.append(pills);
      }
      [["verdict_with_parents", function () { return row.parent_supported_relationships + " inherited relationship(s); " + row.mendelian_conflicts + " Mendelian conflict(s)"; }],
       ["verdict_with_reads", function () { return row.read_supported_relationships + " read-supported relationship(s) from " + row.read_fragments + " informative fragments"; }]].forEach(function (spec) {
        const cell = element("td");
        const value = row[spec[0]];
        if (value) {
          cell.append(element("span", "tier-badge tier-" + ({ good: 4, none: 3, mixed: 1, conflicting: 5 }[value]), VERDICTS[value]));
          cell.title = spec[1]();
        } else cell.textContent = "—";
        cells[spec[0]] = cell;
      });
      cells.heterozygous = element("td", "num", row.variants === undefined || row.variants === "" ? "—" : row.variants + (row.heterozygous ? " · " + row.heterozygous + " het" : ""));
      cells.arrangements = element("td", "num", row.arrangements === undefined || row.arrangements === "" ? "—" : number(row.arrangements));
      cells.unknown_cds_bases = element("td", "num", row.unknown_cds_bases === undefined || row.unknown_cds_bases === "" ? "—" : number(row.unknown_cds_bases));
      COLUMNS.forEach(function (column) { tr.append(cells[column.key]); });
      body.append(tr);
    });
    byId("scan-empty").hidden = visible.length > 0;
    const more = byId("scan-more");
    more.hidden = visible.length <= state.limit;
    more.textContent = "Show " + number(Math.min(300, visible.length - state.limit)) + " more of " + number(visible.length - state.limit);
  }

  renderHeader();
  renderFilters();
  renderHead();
  renderRows();
  byId("search").addEventListener("input", function (event) { state.query = event.target.value; state.limit = 300; renderRows(); });
  byId("scan-more").addEventListener("click", function () { state.limit += 300; renderRows(); });
})();
