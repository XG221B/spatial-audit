# Journal of Biomedical Informatics (JBI) 投稿操作清单

本清单整理了在 Elsevier Editorial Manager 系统提交 *Journal of Biomedical Informatics*（JBI）的全部步骤与避坑要点，按顺序照做即可（**0 元版面费**）。

---

## 投稿前 2 分钟必做核查
1. 打开科睿唯安官方 **[Master Journal List](https://mjl.clarivate.com/)**，搜索 `Journal of Biomedical Informatics`（ISSN: `1532-0464`），确认状态仍为 **SCIE Active** 且**无 `On Hold` 标记**。
2. 打开 Elsevier JBI 官网，确认仍支持 **Subscription（传统订阅模式，作者 0 元）**。

---

## 需要上传的文件清单（全部已在本地生成完毕）
| Editorial Manager 栏位 | 上传文件 | 说明 |
|---|---|---|
| **Cover Letter** | [COVER_LETTER.md](COVER_LETTER.md)（复制文本或另存 Word/PDF） | 已针对 JBI“通用生物医学信息学方法学”定位定制 |
| **Manuscript** | [paper/main.pdf](../main.pdf)（或连同 `paper/` 下的 `.tex` + `figures/` + `references.bib`） | 已重构为 JBI 完整 IMRaD 长文结构，含结构化摘要与 Statement of Significance |
| **Highlights** | [HIGHLIGHTS.md](HIGHLIGHTS.md)（复制 5 条要点，每条均 ≤ 85 字符） | 严格符合 Elsevier 字数限制 |
| **Statement of Significance** | [STATEMENT_OF_SIGNIFICANCE.md](STATEMENT_OF_SIGNIFICANCE.md) | JBI 特色必填项（3 段式） |
| **Supplementary Material** | [paper/supplement.pdf](../supplement.pdf) | 补充方法、离散扩散推导与完整复现流程 |

---

## 表单填写关键项（防止手滑踩坑）
- **Article Type**：选 `Original Research Article`（或 `Methodological Review` / `Original Research`）。
- **Publishing Option（最关键！）**：务必选择 **Subscription（Traditional publishing model）= 0 元**，**不要勾选 Gold Open Access**。
- **Author & ORCID**：`Chaohao Liu`，ORCID `0009-0003-4341-2494`，Affiliation 填 `Independent Researcher`（若挂学校单位，先确认学院是否要求导师签字）。
- **CRediT Author Statement**：`Chaohao Liu: Conceptualization, Methodology, Software, Validation, Formal analysis, Data Curation, Writing – Original Draft, Writing – Review & Editing, Visualization.`
- **Declaration of Competing Interests**：`The author declares that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.`
- **Declaration of Generative AI and AI-assisted Technologies**：可直接粘贴正文中的 `AI usage disclosure` 声明（Elsevier 要求透明披露，已写好）。
- **Data & Code Availability**：
  - GitHub 仓库：`https://github.com/XG221B/spatial-audit`
  - 在线交互看板：`https://xg221b.github.io/spatial-audit/`
