"use strict";
const $ = (id) => document.getElementById(id);
let language = "en";
const zh = {
  "Original-package R runs on the same registered feature panels. Pseudobulk differential expression, domain-interaction models, and within-section spatial-variability tests target distinct null hypotheses. Peak RSS records the sampled process-tree high-water mark.":
    "在相同的预注册特征面板上运行原版 R 包。伪批量差异表达、区域交互模型和切片内空间变异检验各自对应不同的零假设。峰值 RSS 为采样期间进程树的最大驻留内存。",
  "Explore how donor count, within-patient correlation, tissue composition, and baseline abundance shifts affect exact permutation resolution, naive pooled Type I error, and fine-scale energy. Adjust any slider to recalculate.":
    "调整滑块，查看供体数、患者内相关性、组织组成和基线丰度变化对精确置换分辨率、合并分析第一类错误率及细尺度能量的影响。",
  "Spatial Audit · MIT License · Offline-ready static observatory":
    "Spatial Audit · MIT 许可证 · 支持离线浏览",
  "BH call rates": "BH 检出率",
  "Multiscale share": "多尺度协变量占比",
  "Section centers by group": "按组显示切片中心",
  "Naive null rejection curve and nominal five-percent reference":
    "零假设拒绝率曲线与名义 5% 参考线",
  "Exact floor by biological sample size": "生物学样本量与精确置换下限",
  "Raw and corrected illustrative energy": "原始与校正能量示例",
  "Patients with interpretable covariate projections by band: ":
    "各波段可解释协变量投影的患者数：",
  "Unassessed bands reflect footprint size, grid resolution, missing covariates, or projection saturation.":
    "未评估的波段受空间覆盖、网格分辨率、缺失协变量或投影饱和影响。",
  "Preparation status": "准备状态",
  "section centers": "切片中心",
  exhaustive: "穷举",
  executed: "实际执行",
  "First feasible BH rank under the executed budget":
    "当前预算下首个可达的 BH 秩",
  none: "无",
  "Unsupported audit report": "不支持的审计报告",
  "Opened locally": "已在本地打开",
  "No measured comparator runs available.": "暂无实测比较结果。",
  Raw: "原始",
  Adjusted: "校正",
  "1 cell": "1 个细胞",
  " cells per patient": " 个细胞 / 患者",
  "2 per group": "每组 2 人",
  "30 per group · log10 p scale": "每组 30 人 · p 值对数刻度",
  "Download full minimum-donor lookup table (CSV) ↓":
    "下载完整的最小供体数查询表（CSV）↓",
  Composition: "组成",
  Dispersion: "离散度",
  Mean: "均值",
  "Batch alignment": "批次对齐",
  "Minimum N / arm": "每组最小样本量",
  "Orange shows analytic Gaussian null rejection under the selected intrapatient correlation; green marks the nominal 5% level. The illustrative energy model partitions unit biological energy by composition share and adds a Poisson noise term equal to the mean multiplier. The empirical simulation grid below evaluates the full spatial estimator.":
    "橙色表示指定患者内相关性下的高斯零假设解析拒绝率，绿色表示名义 5% 水平。示例模型按组成占比分配单位生物学能量，并加入等于均值倍数的泊松噪声项。下方模拟评估完整的空间估计方法。",
  " synthetic cohorts. Raw rejection evaluates total fine-band energy; adjusted rejection evaluates residual energy after composition projection and noise subtraction. Under the null, latent biological amplitude is matched across groups while baseline mean, composition, or batch alignment may vary. ":
    " 个合成队列。原始拒绝率基于细波段总能量，校正拒绝率基于组成投影与噪声扣除后的残余能量。零假设下两组的潜在生物学振幅相同，基线均值、组成和批次对齐可变化。",
  "Inspect any spatial-audit JSON report locally in your browser":
    "在浏览器中本地查看 spatial-audit JSON 报告",
  Main: "主导航",
  "Audit status legend": "审计状态图例",
  "Dimension 1: Exact permutation floor and feasible BH rank":
    "维度 1：精确置换下限和可达 BH 秩",
  "Dimension 2: Slide, batch, and acquisition-stage collinearity":
    "维度 2：载玻片、批次与采集坐标共线性",
  "Dimension 3: Pooled vs. patient-level BH discovery gap":
    "维度 3：合并分析与患者水平分析的 BH 检出率差值",
  "Dimension 4: Within-band composition or TIC energy share":
    "维度 4：波段内组成或总离子流能量占比",
  "Dimension 5: High-frequency short-lag noise sensitivity":
    "维度 5：高频短程噪声敏感性",
  Cohorts: "队列",
  "R Benchmarks": "R 基准测试",
  "Trap Sandbox": "统计模拟",
  "↑ Open report": "↑ 打开报告",
  "Multi-patient Spatial Omics · Reproducibility & Confounding Audit":
    "多患者空间组学 · 可重复性与混杂审计",
  "Cohort observatory": "队列总览",
  "Cohort register": "队列目录",
  "Search cohorts": "搜索队列",
  "Five-dimensional statistical assessment of biological replication, acquisition collinearity, pseudoreplication, physical-scale covariate sharing, and short-lag noise sensitivity across public archives.":
    "从生物学重复、采集设计共线性、伪重复、物理尺度协变量贡献及短程噪声敏感性五个维度评估公共存档数据。",
  "Low flag": "低风险",
  Review: "需复核",
  "High flag": "高风险",
  Unassessed: "未评估",
  "Indexed archives": "已收录存档",
  "Platform families": "平台类型",
  "Public audits": "公开审计",
  "Synthetic controls": "合成对照",
  "Primary human cohorts": "主要人类队列",
  "Additional public analyses": "其他公开分析",
  "Supplementary archives": "补充存档",
  "Synthetic fixtures": "合成数据",
  "Local uploads": "本地报告",
  "Filter by accession, platform, or contrast…": "按编号、平台或组间对比筛选…",
  "Cohort / platform": "队列 / 平台",
  Power: "功效",
  Design: "设计",
  "Rep.": "伪重复",
  "Comp.": "组成",
  Nugget: "块金效应",
  "Cross-language R benchmarks": "跨语言 R 基准测试",
  "Statistical trap sandbox": "统计陷阱模拟",
  "Analytic illustration": "解析示例",
  "Simulation phase diagram & minimum-donor lookup table":
    "模拟相图与最小供体数量表",
  "Cohort catalog": "队列目录",
  "Input contract": "输入规范",
  "Validation notes": "验证说明",
  "Public archive · measured": "公共存档 · 已测量",
  "Synthetic Benchmark Fixture": "合成基准数据",
  "Public archive · catalog only": "公共存档 · 仅目录",
  "Local upload · self-declared ": "本地报告 · 声明类型 ",
  "No cohorts match the current filter.": "没有符合筛选条件的队列。",
  " shown": " 条记录",
  "Audited contrast:": "审计分组：",
  "Design flags:": "设计风险：",
  "Permutation resolution:": "置换分辨率：",
  "biological units": "生物学单位",
  sections: "切片",
  observations: "观测",
  features: "特征",
  "Download audit JSON ↓": "下载审计 JSON ↓",
  "Input receipts, mapping & calibration ↗": "输入记录、映射与校准 ↗",
  "BH discovery rate by unit of analysis": "不同分析单位的 BH 检出率",
  "Pooled spots": "合并位点",
  Patients: "患者",
  "Covariate share across physical wavelengths": "各物理波长下的协变量能量占比",
  "Acquisition coordinates": "采集坐标",
  "Dimension-by-dimension audit notes": "各维度审计说明",
  "Archive & publication provenance": "存档与文献来源",
  "Broader study design:": "原研究设计：",
  "Official archive ↗": "官方存档 ↗",
  "Original publication ↗": "原始文献 ↗",
  "Exact permutation floor": "精确置换下限",
  "Acquisition design": "采集设计",
  Pseudoreplication: "伪重复",
  "Multiscale composition": "多尺度组成",
  "HF-nugget sensitivity": "高频噪声敏感性",
  "Per-feature Diagnostics": "逐特征诊断",
  "Search features": "搜索特征",
  Feature: "特征",
  "Pooled q": "合并 q 值",
  "Patient q": "患者 q 值",
  "Raw energy Δ": "原始能量差",
  "Debiased energy Δ": "校正能量差",
  "Nugget flag": "噪声敏感",
  Yes: "是",
  No: "否",
  "No matching features": "无匹配特征",
  "Export SVG": "导出 SVG",
  "Patients per group": "每组患者数",
  "Cells per patient": "每位患者的细胞数",
  "Composition share": "组成占比",
  "Baseline mean multiplier": "基线均值倍数",
  "Within-patient correlation (ICC)": "患者内相关系数（ICC）",
  "Naive pooled null rejection": "合并分析的零假设拒绝率",
  "Adjusted / raw energy": "校正 / 原始能量",
  "Permutation resolution by donor count": "供体数量与置换分辨率",
  "Illustrative fine-scale energy decomposition": "细尺度能量分解示例",
  "Method & estimand": "方法与估计目标",
  "BH calls": "BH 检出率",
  "Spot gap": "位点差值",
  "Wall (s)": "用时（秒）",
  "Peak RSS (MiB)": "峰值内存（MiB）",
  Status: "状态",
  Source: "来源",
  "Benchmark cohort": "基准队列",
  Previous: "上一页",
  Next: "下一页",
};
const t = (key) => (language === "zh" ? (zh[key] ?? key) : key);
function applyStaticLanguage() {
  const local = $("kind").querySelector('option[value="local"]');
  if (local) local.textContent = t("Local uploads");
  for (const [node, key] of staticText) node.nodeValue = t(key);
  for (const [element, attr, key] of staticAttrs)
    element.setAttribute(attr, t(key));
  document.documentElement.lang = language === "zh" ? "zh-CN" : "en";
  $("language").setAttribute("aria-pressed", String(language === "zh"));
}
const staticText = [],
  staticAttrs = [];
const walker = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
while (walker.nextNode()) {
  const node = walker.currentNode,
    key = node.nodeValue.trim().replace(/\s+/g, " ");
  if (Object.hasOwn(zh, key)) staticText.push([node, key]);
}
for (const el of document.querySelectorAll(
  "[placeholder], [aria-label], [title]",
))
  for (const attr of ["placeholder", "aria-label", "title"])
    if (Object.hasOwn(zh, el.getAttribute(attr)))
      staticAttrs.push([el, attr, el.getAttribute(attr)]);

function featureTable(root, report) {
  const { pseudorep, nugget } = report.dimensions;
  const noise = new Map(
    (Array.isArray(nugget.features) ? nugget.features : []).map((r) => [
      r.feature,
      r,
    ]),
  );
  const data = (
    Array.isArray(pseudorep.features) ? pseudorep.features : []
  ).map((name, i) => {
    const r = noise.get(name);
    return [
      name,
      pseudorep.spot_q?.[i],
      pseudorep.patient?.q?.[i],
      r?.raw_delta,
      r?.debiased_delta,
      Number.isFinite(r?.debiased_delta) ? r.mean_noise_sensitive : null,
    ];
  });
  const labels = [
    "Feature",
    "Pooled q",
    "Patient q",
    "Raw energy Δ",
    "Debiased energy Δ",
    "Nugget flag",
  ];
  root.innerHTML =
    "<h3>" +
    t("Per-feature Diagnostics") +
    '</h3><input id="feature-search" type="search" aria-label="' +
    t("Search features") +
    '" placeholder="' +
    t("Search features") +
    '"><div class="scroll"><table><thead><tr>' +
    labels
      .map(
        (key) =>
          '<th scope="col"><button type="button">' + t(key) + "</button></th>",
      )
      .join("") +
    '</tr></thead><tbody></tbody></table></div><div class="feature-pages"><button type="button">' +
    t("Previous") +
    '</button><output></output><button type="button">' +
    t("Next") +
    "</button></div>";
  const headings = [...root.querySelectorAll("th")],
    [prev, next] = root.querySelectorAll(".feature-pages button");
  let col = 0,
    direction = 1,
    page = 0;
  function draw() {
    const query = $("feature-search").value.toLowerCase();
    const rows = data
      .filter((row) => String(row[0]).toLowerCase().includes(query))
      .sort((a, b) => {
        const x = a[col],
          y = b[col];
        if (x == null || y == null)
          return Number(x == null) - Number(y == null);
        return (
          direction *
          (typeof x === "string" ? x.localeCompare(y) : Number(x) - Number(y))
        );
      });
    headings.forEach((h, i) =>
      h.setAttribute(
        "aria-sort",
        i === col ? (direction === 1 ? "ascending" : "descending") : "none",
      ),
    );
    const start = page * 50,
      end = Math.min(start + 50, rows.length);
    root.querySelector("tbody").innerHTML = rows
      .slice(start, end)
      .map(
        (row) =>
          "<tr>" +
          row
            .map((v, i) => {
              let display = "—";
              if (i === 0) display = esc(v);
              else if (i === 5 && typeof v === "boolean")
                display = t(v ? "Yes" : "No");
              else if (Number.isFinite(v)) display = v.toPrecision(4);
              return "<td>" + display + "</td>";
            })
            .join("") +
          "</tr>",
      )
      .join("");
    root.querySelector("output").textContent = rows.length
      ? start + 1 + "–" + end + " / " + rows.length
      : t("No matching features");
    prev.disabled = page === 0;
    next.disabled = end === rows.length;
  }
  headings.forEach(
    (h, i) =>
      (h.querySelector("button").onclick = () => {
        direction = i === col ? -direction : 1;
        col = i;
        page = 0;
        draw();
      }),
  );
  $("feature-search").oninput = () => {
    page = 0;
    draw();
  };
  prev.onclick = () => {
    page--;
    draw();
  };
  next.onclick = () => {
    page++;
    draw();
  };
  draw();
}

function chartDownloads(root, id) {
  root.querySelectorAll(".chart svg").forEach((plot, i) => {
    const button = document.createElement("button");
    button.className = "svg-download";
    button.type = "button";
    button.textContent = t("Export SVG");
    button.onclick = () => {
      const copy = plot.cloneNode(true);
      copy.querySelectorAll("text").forEach((node) => {
        node.setAttribute("font-size", "11");
        node.setAttribute("font-family", "sans-serif");
        node.setAttribute("fill", "#3d5254");
      });
      const url = URL.createObjectURL(
        new Blob([new XMLSerializer().serializeToString(copy)], {
          type: "image/svg+xml",
        }),
      );
      const link = document.createElement("a");
      link.href = url;
      link.download = id + "-chart-" + (i + 1) + ".svg";
      link.click();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    };
    plot.after(button);
  });
}

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
  Number.isFinite(x) ? (100 * x).toFixed(1) + "%" : t("Unassessed");
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
  esc(t(label)) +
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
        t(x[1]) +
        "</span></div>",
    )
    .join("");
  $("total").textContent = rows.length + t(" shown");
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
    out.textContent = t("No cohorts match the current filter.");
    return;
  }
  const { r, c } = current;
  let html =
    '<div class="tag">' +
    ($("kind").value === "local"
      ? t("Local upload · self-declared ") + esc(r.kind)
      : r
      ? r.kind === "public"
        ? t("Public archive · measured")
        : t("Synthetic Benchmark Fixture")
      : t("Public archive · catalog only")) +
    "</div><h2>" +
    esc(r?.id || current.id) +
    "</h2>";
  if (!r) {
    if (c) html += "<p>" + esc(c.contrast) + "</p>";
    html +=
      "<p>" + t("Preparation status") + " · " + esc(c?.preparation_status || c?.unit_note) + "</p>";
  } else {
    const d = r.dimensions;
    const pwr = d.power || {};
    if (Array.isArray(pwr.groups) && pwr.groups.length === 2 && Number.isFinite(pwr.n0) && Number.isFinite(pwr.n1)) {
      html +=
        "<p><b>" + t("Audited contrast:") + "</b> " +
        esc(pwr.groups[0]) + " (n=" + num(pwr.n0) + ") vs " +
        esc(pwr.groups[1]) + " (n=" + num(pwr.n1) + ")</p>";
    }
    if (r.provenance?.selection) {
      html += '<p class="muted">' + esc(r.provenance.selection) + "</p>";
    }
    const evidence = (sources[0].evidence || {})[r.id];
    const links = [];
    if ($("kind").value !== "local" && /^[A-Za-z0-9_.-]+$/.test(r.id))
      links.push('<a download href="data/' + r.kind + '/' + esc(r.id) + '.json">' + t("Download audit JSON ↓") + '</a>');
    if ($("kind").value !== "local" && r.kind === "public" && /^[A-Za-z0-9_.-]+\.json$/.test(evidence || ""))
      links.push('<a href="data/public/' + esc(evidence) + '">' + t("Input receipts, mapping & calibration ↗") + '</a>');
    if (links.length)
      html += '<div class="detail-links">' + links.join("") + "</div>";
    html +=
      '<div class="stat-strip">' +
      '<span class="stat-pill"><b>' + num(r.n_patients) + "</b> " + t("biological units") + "</span>" +
      (Number.isFinite(r.n_sections) ? '<span class="stat-pill"><b>' + num(r.n_sections) + "</b> " + t("sections") + "</span>" : "") +
      '<span class="stat-pill"><b>' + num(r.n_observations) + "</b> " + t("observations") + "</span>" +
      '<span class="stat-pill"><b>' + num(r.n_features) + "</b> " + t("features") + "</span>" +
      "</div>" +
      (Array.isArray(d.design.reasons) && d.design.reasons.length ? "<p><b>" + t("Design flags:") + "</b> " + esc(d.design.reasons.map(
        (k) => k === "batch" ? (d.design.categorical?.batch?.label || k) : k
      ).join(", ")) + "</p>" : "") +
      '<div class="chart"><h3>' + t("BH discovery rate by unit of analysis") + '</h3>' +
      svg(
        [
          ["Pooled spots", d.pseudorep.spot_rate, "#db9674"],
          ["Patients", d.pseudorep.patient_rate, "#34786c"],
        ]
          .map(
            ([name, v, color], i) =>
              txt(2, 29 + 44 * i, t(name)) +
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
    const step = 290 / Math.max(1, bands.length - 1);
    const labelEvery = Math.max(1, Math.ceil((bands.length - 1) / 2));
    html +=
      '<div class="chart"><h3>' + t("Covariate share across physical wavelengths") + '</h3>' +
      svg(
        bands.slice(1).map((b, i) => Number.isFinite(b.share) && Number.isFinite(bands[i].share)
          ? line([[42 + step * i, 81 - bands[i].share * 60], [42 + step * (i + 1), 81 - b.share * 60]], "#387c70")
          : "").join("") +
          bands
            .map(
              (b, i) =>
                (i % labelEvery === 0 || i === bands.length - 1 ?
                  txt(22 + step * i, 105, num(b.wavelength_um) + " µm") +
                  txt(22 + step * i, 20, pct(b.share)) : "") +
                (Number.isFinite(b.share)
                  ? '<circle cx="' +
                    (42 + step * i) +
                    '" cy="' +
                    (81 - 60 * b.share) +
                    '" r="4.5" fill="#34786c"><title>' + esc(num(b.wavelength_um) + " µm · " + pct(b.share)) + "</title></circle>"
                  : ""),
            )
            .join(""),
        "Multiscale share",
      ) +
      (bands.every(b => Number.isInteger(b.interpretable_patients)) ?
        "<small>" + t("Patients with interpretable covariate projections by band: ") +
        bands.map(b => num(b.interpretable_patients) + "/" + num(r.n_patients)).join(" · ") +
        ". " + t("Unassessed bands reflect footprint size, grid resolution, missing covariates, or projection saturation.") + "</small>" : "") +
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
      '<div class="chart"><h3>' + t("Acquisition coordinates") + '</h3><small>' +
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
      "<details><summary>" + t("Dimension-by-dimension audit notes") + "</summary><ul>" +
      dims
        .map(
          (k, i) =>
            "<li>" +
            dot(d[k].status) +
            " <b>" +
            t(dimNames[i]) +
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
      html += "<p><b>" + t("Permutation resolution:") + "</b> " + t("exhaustive") + " " + esc(power.p_floor.toPrecision(3)) +
        " · " + t("executed") + " " + esc(power.effective_p_floor.toPrecision(3)) + " (" +
        esc(power.executed_mode) + "). " + t("First feasible BH rank under the executed budget") + ": " +
        esc(power.effective_minimum_bh_rank ?? t("none")) + ".</p>";
    }
  }
  if (c)
    html +=
      "<details><summary>" + t("Archive & publication provenance") + "</summary>" +
      (r && c.contrast ? "<p><b>" + t("Broader study design:") + "</b> " + esc(c.contrast) + "</p>" : "") +
      "<p>" +
      esc(c.claim) +
      "</p><p>" +
      esc(c.method) +
      '</p><div class="detail-links"><a href="' +
      safeURL(c.source) +
      '" rel="noreferrer">' + t("Official archive ↗") + '</a>' +
      (safeURL(c.paper)
        ? '<a href="' +
          safeURL(c.paper) +
          '" rel="noreferrer">' + t("Original publication ↗") + '</a>'
        : "") +
      "</div></details>";
  out.innerHTML = html;
  if (r) {
    const table = document.createElement("section");
    table.className = "feature-diagnostics";
    out.append(table);
    featureTable(table, r);
    chartDownloads(out, r.id);
  }
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
      !["public", "synthetic", "local"].includes(r.kind) ||
      typeof r.id !== "string" ||
      !dims.every((k) =>
        ["green", "red", "yellow"].includes(r.dimensions?.[k]?.status),
      ) ||
      !Array.isArray(r.dimensions.composition.bands) ||
      !Array.isArray(r.dimensions.design.points) ||
      !r.dimensions.composition.bands.length ||
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
      throw Error(t("Unsupported audit report"));
    const key = r.kind + "/" + r.id;
    localReports.set(key, r);
    if (!$("kind").querySelector('option[value="local"]'))
      $("kind").add(new Option(t("Local uploads"), "local"));
    $("kind").value = "local";
    $("search").value = "";
    selected = key;
    render();
    $("notice").textContent = t("Opened locally") + " · " + file.name;
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
      ? "<table><thead><tr><th>" + t("Method & estimand") + "</th><th>" + t("BH calls") + "</th><th>" + t("Spot gap") + "</th><th>" + t("Wall (s)") + "</th><th>" + t("Peak RSS (MiB)") + "</th><th>" + t("Status") + "</th></tr></thead><tbody>" +
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
        t(name) +
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
        ([label, v]) => "<span>" + t(label) + "<strong>" + v + "</strong></span>",
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
          txt(10, 110, t("1 cell")) +
          txt(250, 110, c + t(" cells per patient")),
        "Naive null rejection curve and nominal five-percent reference",
      ) +
      "<small>" + t("Orange shows analytic Gaussian null rejection under the selected intrapatient correlation; green marks the nominal 5% level. The illustrative energy model partitions unit biological energy by composition share and adds a Poisson noise term equal to the mean multiplier. The empirical simulation grid below evaluates the full spatial estimator.") + "</small>";
    const floorY = (k) => 16 - 4.3 * Math.log10(floorAt(k));
    const donorX = (k) => 25 + (360 * (k - 2)) / 28;
    $("curve").insertAdjacentHTML(
      "beforeend",
      '<div class="chart"><h3>' + t("Permutation resolution by donor count") + '</h3>' +
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
            txt(15, 108, t("2 per group")) +
            txt(215, 108, t("30 per group · log10 p scale")),
          "Exact floor by biological sample size",
        ) +
        "</div>" +
        '<div class="chart"><h3>' + t("Illustrative fine-scale energy decomposition") + '</h3>' +
        svg(
          [
            ["Raw", raw],
            ["Adjusted", corrected],
          ]
            .map(
              ([label, v], i) =>
                txt(5, 30 + 40 * i, t(label)) +
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
      '<p><a href="data/synthetic/sample_size_lookup.csv" download>' + t("Download full minimum-donor lookup table (CSV) ↓") + '</a></p><img src="data/synthetic/phase_diagram.svg" alt="Synthetic phase diagram for design and confounding"><p>' +
      num(phase.simulated_cohorts) +
      t(" synthetic cohorts. Raw rejection evaluates total fine-band energy; adjusted rejection evaluates residual energy after composition projection and noise subtraction. Under the null, latent biological amplitude is matched across groups while baseline mean, composition, or batch alignment may vary. ") +
      esc(phase.lookup_rule) +
      "</p><table><thead><tr><th>" + t("Composition") + "</th><th>" + t("Dispersion") + "</th><th>" + t("Mean") + "</th><th>" + t("Batch alignment") + "</th><th>" + t("Minimum N / arm") + "</th></tr></thead><tbody>" +
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
$("language").onclick = () => {
  language = language === "en" ? "zh" : "en";
  const cohort = $("benchselect").value;
  const sliders = [...document.querySelectorAll("#sliders input")].map(el => [el.id, el.value]);
  applyStaticLanguage();
  render(); benchmarks(); sandbox();
  $("benchselect").value = cohort;
  $("benchselect").onchange();
  sliders.forEach(([id, value]) => $(id).value = value);
  if (sliders.length) $(sliders[0][0]).oninput();
};
render();
benchmarks();
sandbox();
