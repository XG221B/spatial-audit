"use strict";
const $ = (id) => document.getElementById(id);
const sources = [
  globalThis.SPATIAL_PUBLIC || {},
  globalThis.SPATIAL_SYNTHETIC || {},
];
const reports = sources.flatMap((s) => s.reports || []);
const localReports = new Map();
const catalog = sources[0].catalog || [];
const dims = ["power", "design", "pseudorep", "composition", "nugget"];
const dimNames = [
  "Exact permutation floor",
  "Acquisition design",
  "Pseudoreplication",
  "Multiscale composition",
  "HF-nugget sensitivity",
];
const esc = (x) =>
  String(x ?? "").replace(
    /[&<>"']/g,
    (c) =>
      ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[
        c
      ],
  );
const pct = (x) =>
  Number.isFinite(x) ? (100 * x).toFixed(1) + "%" : "Unassessed";
const num = (x) =>
  Number.isFinite(x)
    ? x.toLocaleString("en-US", { maximumFractionDigits: 3 })
    : "—";
const safeURL = (u) => (/^https?:\/\//i.test(u || "") ? esc(u) : "");
const status = (x) => (["green", "yellow", "red"].includes(x) ? x : "gray");
const dot = (x) =>
  '<span class="dot ' +
  status(x) +
  '" role="img" aria-label="' +
  esc(x || "Unassessed") +
  '"></span>';
const svg = (body, label) =>
  '<svg viewBox="0 0 410 116" role="img" aria-label="' +
  esc(label) +
  '">' +
  body +
  "</svg>";
const line = (points, color) =>
  '<polyline fill="none" stroke="' +
  color +
  '" stroke-width="2.4" points="' +
  points.map((p) => p.join(",")).join(" ") +
  '"/>';
const txt = (x, y, t) =>
  '<text x="' + x + '" y="' + y + '">' + esc(t) + "</text>";
let selected = reports.find((r) => r.kind === "public")?.id || catalog[0]?.id;

function getRows() {
  const kind = $("kind").value,
    query = $("search").value.trim().toLowerCase();
  if (kind === "local") {
    return Array.from(localReports, ([id, r]) => ({ id, platform: r.platform, r }))
      .filter(x => (x.r.id + " " + x.platform).toLowerCase().includes(query));
  }
  const archiveRows = catalog.map(c => ({
    c, id: c.id, platform: c.platform,
    r: reports.find(r => r.kind === "public" && r.id === c.id) ||
       reports.find(r => r.kind === "public" && c.audit_ids?.includes(r.id)),
  }));
  let rows;
  if (kind === "public" || kind === "supplementary")
    rows = archiveRows.filter(x => (x.c.scope === "human") === (kind === "public"));
  else if (kind === "variants")
    rows = reports.filter(r => r.kind === "public" && !archiveRows.some(x => x.r === r))
      .map(r => ({ r, id: r.id, platform: r.platform, c: catalog.find(c => c.audit_ids?.includes(r.id)) }));
  else rows = reports.filter(r => r.kind === "synthetic")
      .map(r => ({ id: r.id, platform: r.platform, r }));
  return rows.filter((x) =>
    (x.id + " " + x.platform + " " + (x.c?.contrast || ""))
      .toLowerCase()
      .includes(query),
  );
}

function render() {
  const rows = getRows();
  const current =
    rows.find((x) => x.id === selected || x.r?.id === selected) || rows[0];
  selected = current?.id;
  $("metrics").innerHTML = [
    [catalog.length, "Indexed archives"],
    [new Set(catalog.map((c) => c.platform)).size, "Platform families"],
    [reports.filter((r) => r.kind === "public").length, "Public audits"],
    [
      reports.filter((r) => r.kind === "synthetic").length,
      "Synthetic controls",
    ],
  ]
    .map(
      (x) =>
        '<div class="metric"><strong>' +
        x[0] +
        "</strong><span>" +
        x[1] +
        "</span></div>",
    )
    .join("");
  $("total").textContent = rows.length + " shown";
  $("rows").innerHTML = rows
    .map(
      (x, i) =>
        '<tr class="' +
        (x === current ? "active" : "") +
        '"><td><button type="button" data-i="' +
        i +
        '">' +
        esc(x.r && $("kind").value === "local" ? x.r.id : x.id) +
        "</button><small>" +
        esc(x.platform) +
        (x.r?.id.endsWith("-subset") ? " · subset" : "") +
        "</small></td>" +
        dims
          .map((d) => "<td>" + dot(x.r?.dimensions[d]?.status) + "</td>")
          .join("") +
        "</tr>",
    )
    .join("");
  $("rows")
    .querySelectorAll("tr")
    .forEach(
      (tr, i) =>
        (tr.onclick = () => {
          selected = rows[i].id;
          render();
        }),
    );
  const out = $("detail");
  if (!current) {
    out.textContent = "No cohorts match the current filter.";
    return;
  }
  const { r, c } = current;
  let html =
    '<div class="tag">' +
    ($("kind").value === "local"
      ? "Local upload · self-declared " + esc(r.kind)
      : r
      ? r.kind === "public"
        ? "Public archive · measured"
        : "Synthetic Benchmark Fixture"
      : "Public archive · catalog only") +
    "</div><h2>" +
    esc(r?.id || current.id) +
    "</h2>";
  if (!r) {
    if (c) html += "<p>" + esc(c.contrast) + "</p>";
    html +=
      "<p>Preparation status · " + esc(c?.preparation_status || c?.unit_note) + "</p>";
  } else {
    const d = r.dimensions;
    const pwr = d.power || {};
    if (Array.isArray(pwr.groups) && pwr.groups.length === 2 && Number.isFinite(pwr.n0) && Number.isFinite(pwr.n1)) {
      html +=
        "<p><b>Audited contrast:</b> " +
        esc(pwr.groups[0]) + " (n=" + num(pwr.n0) + ") vs " +
        esc(pwr.groups[1]) + " (n=" + num(pwr.n1) + ")</p>";
    }
    if (r.provenance?.selection) {
      html += '<p class="muted">' + esc(r.provenance.selection) + "</p>";
    }
    const evidence = (sources[0].evidence || {})[r.id];
    const links = [];
    if ($("kind").value !== "local" && /^[A-Za-z0-9_.-]+$/.test(r.id))
      links.push('<a download href="data/' + r.kind + '/' + esc(r.id) + '.json">Download audit JSON ↓</a>');
    if ($("kind").value !== "local" && r.kind === "public" && /^[A-Za-z0-9_.-]+\.json$/.test(evidence || ""))
      links.push('<a href="data/public/' + esc(evidence) + '">Input receipts, mapping &amp; calibration ↗</a>');
    if (links.length)
      html += '<div class="detail-links">' + links.join("") + "</div>";
    html +=
      '<div class="stat-strip">' +
      '<span class="stat-pill"><b>' + num(r.n_patients) + "</b> biological units</span>" +
      (Number.isFinite(r.n_sections) ? '<span class="stat-pill"><b>' + num(r.n_sections) + "</b> sections</span>" : "") +
      '<span class="stat-pill"><b>' + num(r.n_observations) + "</b> observations</span>" +
      '<span class="stat-pill"><b>' + num(r.n_features) + "</b> features</span>" +
      "</div>" +
      (Array.isArray(d.design.reasons) && d.design.reasons.length ? "<p><b>Design flags:</b> " + esc(d.design.reasons.map(
        (k) => k === "batch" ? (d.design.categorical?.batch?.label || k) : k
      ).join(", ")) + "</p>" : "") +
      '<div class="chart"><h3>BH discovery rate by unit of analysis</h3>' +
      svg(
        [
          ["Pooled spots", d.pseudorep.spot_rate, "#db9674"],
          ["Patients", d.pseudorep.patient_rate, "#34786c"],
        ]
          .map(
            ([name, v, color], i) =>
              txt(2, 29 + 44 * i, name) +
              '<rect x="98" y="' +
              (14 + i * 44) +
              '" width="248" height="18" rx="2" fill="#edf1ef"/><rect x="98" y="' +
              (14 + i * 44) +
              '" width="' +
              248 * v +
              '" height="18" rx="2" fill="' +
              color +
              '"/>' +
              txt(350, 29 + i * 44, pct(v)),
          )
          .join(""),
        "BH call rates",
      ) +
      "</div>";
    const bands = d.composition.bands || [];
    html +=
      '<div class="chart"><h3>Covariate share across physical wavelengths</h3>' +
      svg(
        bands.slice(1).map((b, i) => Number.isFinite(b.share) && Number.isFinite(bands[i].share)
          ? line([[42 + 140 * i, 81 - bands[i].share * 60], [182 + 140 * i, 81 - b.share * 60]], "#387c70")
          : "").join("") +
          bands
            .map(
              (b, i) =>
                txt(22 + 140 * i, 105, num(b.wavelength_um) + " µm") +
                txt(22 + 140 * i, 20, pct(b.share)) +
                (Number.isFinite(b.share)
                  ? '<circle cx="' +
                    (42 + 140 * i) +
                    '" cy="' +
                    (81 - 60 * b.share) +
                    '" r="4.5" fill="#34786c"/>'
                  : ""),
            )
            .join(""),
        "Multiscale share",
      ) +
      (bands.every(b => Number.isInteger(b.interpretable_patients)) ?
        "<small>Patients with interpretable covariate projections by band: " +
        bands.map(b => num(b.interpretable_patients) + "/" + num(r.n_patients)).join(" · ") +
        ". Unassessed bands reflect footprint size, grid resolution, missing covariates, or projection saturation.</small>" : "") +
      "</div>";
    const pts = d.design.points || [],
      groups = [...new Set(pts.map((p) => p.group))];
    const range = (k) => pts.reduce(
      ([lo, hi], p) => [Math.min(lo, p[k]), Math.max(hi, p[k])],
      [Infinity, -Infinity],
    );
    const xx = range("x"),
      yy = range("y"),
      scale = (v, [a, b]) => (b > a ? (v - a) / (b - a) : 0.5);
    const stride = Math.max(1, Math.ceil(pts.length / 2400));
    const plotted = pts.filter((p, i) => i % stride === 0);
    html +=
      '<div class="chart"><h3>Acquisition coordinates</h3><small>' +
      esc(d.design.points_frame) +
      " · µm" + (stride > 1 ? " · showing " + num(plotted.length) + " of " + num(pts.length) + " section centers" : "") + "</small>" +
      svg(
        plotted
          .map(
            (p) =>
              '<circle cx="' +
              (30 + 350 * scale(p.x, xx)) +
              '" cy="' +
              (83 - 65 * scale(p.y, yy)) +
              '" r="5" fill="' +
              (groups.indexOf(p.group) ? "#b16e4b" : "#31786d") +
              '" opacity=".78"><title>' +
              esc([p.patient, p.group, p.x, p.y].join(" · ")) +
              "</title></circle>",
          )
          .join("") +
          txt(15, 110, "x " + num(xx[0]) + " – " + num(xx[1])) +
          txt(212, 110, "y " + num(yy[0]) + " – " + num(yy[1])),
        "Section centers by group",
      ) +
      groups
        .map(
          (g, i) =>
            '<span class="swatch" style="background:' +
            (i ? "#b16e4b" : "#31786d") +
            '"></span>' +
            esc(g),
        )
        .join("") +
      "</div>" +
      "<details><summary>Dimension-by-dimension audit notes</summary><ul>" +
      dims
        .map(
          (k, i) =>
            "<li>" +
            dot(d[k].status) +
            " <b>" +
            dimNames[i] +
            ":</b> " +
            esc(d[k].message) +
            (k === "nugget" && d[k].measurement_scope ? " " + esc(d[k].measurement_scope.interpretation) : "") +
            "</li>",
        )
        .join("") +
      "</ul><p>" +
      esc(r.statement) +
      "</p></details>";
    const power = d.power;
    if (Number.isFinite(power.effective_p_floor) && Number.isFinite(power.p_floor)) {
      html += "<p><b>Permutation resolution:</b> exhaustive " + esc(power.p_floor.toPrecision(3)) +
        " · executed " + esc(power.effective_p_floor.toPrecision(3)) + " (" +
        esc(power.executed_mode) + "). First feasible BH rank under the executed budget: " +
        esc(power.effective_minimum_bh_rank ?? "none") + ".</p>";
    }
  }
  if (c)
    html +=
      "<details><summary>Archive &amp; publication provenance</summary>" +
      (r && c.contrast ? "<p><b>Broader study design:</b> " + esc(c.contrast) + "</p>" : "") +
      "<p>" +
      esc(c.claim) +
      "</p><p>" +
      esc(c.method) +
      '</p><div class="detail-links"><a href="' +
      safeURL(c.source) +
      '" rel="noreferrer">Official archive ↗</a>' +
      (safeURL(c.paper)
        ? '<a href="' +
          safeURL(c.paper) +
          '" rel="noreferrer">Original publication ↗</a>'
        : "") +
      "</div></details>";
  out.innerHTML = html;
}
$("kind").onchange = render;
$("search").oninput = render;
$("upload").onchange = async (e) => {
  try {
    const file = e.target.files[0];
    if (!file) return;
    const r = JSON.parse(await file.text());
    if (
      r.schema_version !== 1 ||
      !["public", "synthetic"].includes(r.kind) ||
      typeof r.id !== "string" ||
      !dims.every((k) =>
        ["green", "red", "yellow"].includes(r.dimensions?.[k]?.status),
      ) ||
      !Array.isArray(r.dimensions.composition.bands) ||
      !Array.isArray(r.dimensions.design.points) ||
      r.dimensions.composition.bands.length !== 3 ||
      !r.dimensions.design.points.every(
        (p) => p && Number.isFinite(p.x) && Number.isFinite(p.y),
      ) ||
      !r.dimensions.composition.bands.every(
        (b) =>
          b &&
          Number.isFinite(b.wavelength_um) &&
          (b.share === null ||
            (Number.isFinite(b.share) && b.share >= 0 && b.share <= 1)),
      ) ||
      !["spot_rate", "patient_rate"].every(
        (k) =>
          Number.isFinite(r.dimensions.pseudorep[k]) &&
          r.dimensions.pseudorep[k] >= 0 &&
          r.dimensions.pseudorep[k] <= 1,
      )
    )
      throw Error("Unsupported audit report");
    const key = r.kind + "/" + r.id;
    localReports.set(key, r);
    if (!$("kind").querySelector('option[value="local"]'))
      $("kind").add(new Option("Local uploads", "local"));
    $("kind").value = "local";
    $("search").value = "";
    selected = key;
    render();
    $("notice").textContent = "Opened locally · " + file.name;
  } catch (err) {
    $("notice").textContent = err.message;
  }
  e.target.value = "";
};
function benchmarks() {
  const arr = sources[0].sota || [];
  const ids = [...new Set(arr.map((b) => b.cohort_id))];
  $("benchselect").innerHTML = ids
    .map((id) => "<option>" + esc(id) + "</option>")
    .join("");
  const draw = () => {
    $("benchplot").innerHTML = arr.length
      ? "<table><thead><tr><th>Method &amp; estimand</th><th>BH calls</th><th>Spot gap</th><th>Wall (s)</th><th>Peak RSS (MiB)</th><th>Status</th></tr></thead><tbody>" +
        arr
          .filter((b) => b.cohort_id === $("benchselect").value)
          .map(
            (b) =>
              "<tr><td><b>" +
              esc(b.method) +
              "</b><small>" +
              esc(b.estimand) +
              "</small></td><td>" +
              pct(b.call_rate) +
              (Number.isFinite(b.call_rate)
                ? '<i class="bar" style="width:' +
                  Math.min(100, b.call_rate * 100) +
                  '%"></i>'
                : "") +
              "</td><td>" +
              (["DESpace2", "nnSVG", "SPARK-X"].includes(b.method) || b.status === "not_applicable" ? "N/A" : pct(b.spot_minus_method)) +
              "</td><td>" +
              num(b.wall_seconds) +
              "</td><td>" +
              num(
                b.peak_rss_bytes == null ? null : b.peak_rss_bytes / 1048576,
              ) +
              '</td><td title="' +
              esc(b.reason) +
              '">' +
              esc(b.status) +
              "</td></tr>",
          )
          .join("") +
        "</tbody></table>"
      : "No measured comparator runs available.";
  };
  $("benchselect").onchange = draw;
  draw();
}
function sandbox() {
  const controls = [
    ["n", "Patients per group", 2, 30, 6, 1],
    ["cells", "Cells per patient", 10, 1000, 150, 10],
    ["share", "Composition share", 0, 1, 0.4, 0.05],
    ["mean", "Baseline mean multiplier", 1, 4, 1.5, 0.1],
    ["icc", "Within-patient correlation (ICC)", 0, 1, 0.15, 0.01],
  ];
  $("sliders").innerHTML = controls
    .map(
      ([id, name, min, max, val, step]) =>
        "<label>" +
        name +
        '<output id="' +
        id +
        'v"></output><input id="' +
        id +
        '" type="range" min="' +
        min +
        '" max="' +
        max +
        '" value="' +
        val +
        '" step="' +
        step +
        '"></label>',
    )
    .join("");
  // Standard normal CDF approximation, A&S 26.2.17.
  const tail = (x) => {
    const t = 1 / (1 + 0.2316419 * x);
    return (
      0.39894228 *
      Math.exp((-x * x) / 2) *
      t *
      (0.31938153 +
        t *
          (-0.356563782 +
            t * (1.781477937 + t * (-1.821255978 + 1.330274429 * t))))
    );
  };
  const draw = () => {
    const [n, c, r, m, icc] = controls.map(([k]) => {
      const v = Number($(k).value);
      $(k + "v").textContent = v;
      return v;
    });
    const floorAt = (k) => {
      let bin = 1;
      for (let i = 1; i <= k; i++) bin *= (k + i) / i;
      return 2 / bin;
    };
    const err = (k) => 2 * tail(1.959963985 / Math.sqrt(1 + (k - 1) * icc));
    const raw = 1 + m,
      corrected = raw - r - m;
    $("readouts").innerHTML = [
      ["Exact permutation floor", floorAt(n).toPrecision(3)],
      ["Naive pooled null rejection", pct(err(c))],
      ["Adjusted / raw energy", corrected.toFixed(2) + " / " + raw.toFixed(2)],
    ]
      .map(
        ([label, v]) => "<span>" + label + "<strong>" + v + "</strong></span>",
      )
      .join("");
    $("curve").innerHTML =
      svg(
        line(
          Array.from({ length: 41 }, (_, i) => [
            25 + i * 9,
            86 - 70 * err(1 + (i * (c - 1)) / 40),
          ]),
          "#b56b47",
        ) +
          line(
            [
              [25, 82.5],
              [385, 82.5],
            ],
            "#377b70",
          ) +
          txt(10, 110, "1 cell") +
          txt(250, 110, c + " cells per patient"),
        "Naive null rejection curve and nominal five-percent reference",
      ) +
      "<small>Orange shows analytic Gaussian null rejection under the selected intrapatient correlation; green marks the nominal 5% level. The illustrative energy model partitions unit biological energy by composition share and adds a Poisson noise term equal to the mean multiplier. The empirical simulation grid below evaluates the full spatial estimator.</small>";
    const floorY = (k) => 16 - 4.3 * Math.log10(floorAt(k));
    const donorX = (k) => 25 + (360 * (k - 2)) / 28;
    $("curve").insertAdjacentHTML(
      "beforeend",
      '<div class="chart"><h3>Permutation resolution by donor count</h3>' +
        svg(
          line(
            Array.from({ length: 29 }, (_, i) => [
              donorX(i + 2),
              floorY(i + 2),
            ]),
            "#31786d",
          ) +
            '<circle cx="' +
            donorX(n) +
            '" cy="' +
            floorY(n) +
            '" r="5" fill="#b56b47"/>' +
            txt(15, 108, "2 per group") +
            txt(215, 108, "30 per group · log10 p scale"),
          "Exact floor by biological sample size",
        ) +
        "</div>" +
        '<div class="chart"><h3>Illustrative fine-scale energy decomposition</h3>' +
        svg(
          [
            ["Raw", raw],
            ["Adjusted", corrected],
          ]
            .map(
              ([label, v], i) =>
                txt(5, 30 + 40 * i, label) +
                '<rect x="92" y="' +
                (14 + i * 40) +
                '" width="' +
                (245 * v) / 5 +
                '" height="18" rx="2" fill="' +
                (i ? "#34786c" : "#bd7957") +
                '"/>' +
                txt(345, 30 + 40 * i, v.toFixed(2)),
            )
            .join(""),
          "Raw and corrected illustrative energy",
        ) +
        "</div>",
    );
  };
  controls.forEach(([k]) => ($(k).oninput = draw));
  draw();
  const phase = sources[1].phase;
  if (phase)
    $("phase").innerHTML =
      '<p><a href="data/synthetic/sample_size_lookup.csv" download>Download full minimum-donor lookup table (CSV) ↓</a></p><img src="data/synthetic/phase_diagram.svg" alt="Synthetic phase diagram for design and confounding"><p>' +
      num(phase.simulated_cohorts) +
      " synthetic cohorts. Raw rejection evaluates total fine-band energy; adjusted rejection evaluates residual energy after composition projection and noise subtraction. Under the null, latent biological amplitude is matched across groups while baseline mean, composition, or batch alignment may vary. " +
      esc(phase.lookup_rule) +
      "</p><table><thead><tr><th>Composition</th><th>Dispersion</th><th>Mean</th><th>Batch alignment</th><th>Minimum N / arm</th></tr></thead><tbody>" +
      phase.minimum_n
        .filter((r) =>
          Object.entries(phase.figure_slice).every(([k, v]) => r[k] === v),
        )
        .map(
          (r) =>
            "<tr>" +
            [
              r.composition_setting,
              r.nb_dispersion,
              r.case_mean_multiplier,
              r.batch_alignment,
              r.minimum_tested_n_per_group,
            ]
              .map((v) => "<td>" + num(v) + "</td>")
              .join("") +
            "</tr>",
        )
        .join("") +
      "</tbody></table>";
}
render();
benchmarks();
sandbox();
