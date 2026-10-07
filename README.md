# HallO-diGly — diGly limpa Occupancy Shiny App

Modular Shiny application for label-free diGly (ubiquitin remnant) ubiquitinomics
using DIA-NN parquet output and the [`limpa`](https://bioconductor.org/packages/limpa)
framework. It ports the 031–034 `*_limpa_occupancy_*_vs_rest.Rmd` workflows into a
reusable app: any set of cell lines / conditions that can be encoded as a single
grouping factor in the run names, not just the PDAC panel.

Version 2.0 · README revised 2026-10-07

---

## 1. Quick start

```r
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limpa", "limma", "Biostrings", "mixOmics"))

# CRAN
install.packages(c(
  "shiny", "bslib", "shinyFiles", "DT", "arrow",
  "dplyr", "tidyr", "tibble", "stringr", "readr",
  "ggplot2", "ggrepel", "pheatmap", "svglite"
))

shiny::runApp("path/to/digly_shiny")
```

R ≥ 4.3. `mixOmics` is needed for the PCA / PLS-DA QC and export.

**Typical session**

1. **Data** → set paths, token map, select runs, *Load*.
2. **Configuration** → choose target line / mode / thresholds → **click *Apply Configuration*** (required, see §4).
3. **Filtering → Normalization → Occupancy → DE Analysis** → run each tab in order.
4. **Plots / Export / PSP / GSEA .rnk / Proteome integration** as needed.

To skip the pipeline, use **Restore session** with a previously exported RDS.

---

## 2. Module structure

```
app.R                  entry point; wires modules, status indicators, shared reactiveVals
R/
  progress_helpers.R   inline spinners + named progress-bar stepper
  mod_data.R           paths, run scanning, run annotation table, data loading (+ hybrid_path_ui)
  mod_config.R         design, identification filters, thresholds, highlighted genes
  mod_filtering.R      contaminants, GG filters, missingness, site mapping, DPC
  mod_normalization.R  median centring, DPC-Quant, PCA / PLS-DA QC
  mod_occupancy.R      occupancy correction, reproducibility filter, tiering
  mod_de.R             design/contrast, dpcDE ×3 layers, classification, consistency
  mod_plots.R          volcano (make_volcano), heatmap, acquisition QC
  mod_export.R         TSV / RDS / SVG export + browser downloads
  mod_psp.R            PhosphoSitePlus annotation
  mod_rnk.R            GSEA .rnk rank files
  mod_proteome.R       proteome integration (degradative vs non-degradative)
  mod_restore.R        restore a saved analysis RDS
www/logo.png
```

Data flow: `data → cfg → filtering → normalization → occupancy → DE`. Normalization,
occupancy, DE and config results are held in shared `reactiveVal`s so that
**Restore session** can overwrite them and Plots / Export / PSP / .rnk / Proteome
work without re-running.

---

## 3. Tabs

| Tab | What it does |
|-----|--------------|
| **🔄 Restore session** | Upload an RDS written by *Export → Full analysis RDS*. Restores config, quantification, occupancy and DE results (needs `hallo_version` ≥ 2.0 for a full restore; older files restore partially and the log says what is missing). |
| **📂 Data** | Browse or paste paths (diGly parquet, proteome parquet, FASTA, results/figures/limpa directories). Enter the run-name **token map** (`TOKEN=Cell line`, one per line). Scanning starts automatically once both parquet paths exist and again whenever the token map changes (which resets the run tables). Editable run tables (cell line, replicate) with select all / none / unparsed; unparsed rows are highlighted. *Load* reads the selected runs with `readDIANN` (no q-value filtering at read time), pairs diGly and proteome runs by (cell line, replicate), and reads the FASTA. |
| **⚙️ Configuration** | Comparison mode (*one vs rest* or *one vs control*), target / control line (populated from loaded data), precursor q, PTM q, site confidence, min obs / line, drop C-terminal GG, drop ubiquitin genes, FDR and logFC thresholds, diGly and proteome normalisation, `treat()` vs `eBayes()`, genes to highlight. **Settings take effect only after *Apply Configuration*.** |
| **🔬 Filtering** | Contaminants (cRAP prefix match, *Keep GFP* rescue, custom substrings) → precursor q-value → GG (UniMod:121) → PTM q → site confidence → not C-terminal → not ubiquitin gene → proteotypic. Precursor missingness filter (absolute / fraction / none). Proteome: q-value, proteotypic, single-protein groups. GG positions mapped to protein sequences from the FASTA; one row per site per precursor. DPC estimated per layer. |
| **📐 Normalization** | Median centring (DPC re-estimated afterwards), DPC-Quant at site and protein level, summary tables (NPrec, PropObs, mean SE), PCA / PLS-DA scores plot (layer, method, components selectable; NA and zero-variance features removed). |
| **🔗 Occupancy** | Site–protein matching (exact, then leading accession), occupancy = site − protein, per-line reproducibility filter, tier assignment, occupancy distribution plot. |
| **📊 DE Analysis** | Design matrix (`~0 + cell_line`) and contrast; `dpcDE` → `contrasts.fit` → `treat`/`eBayes` at site, occupancy and protein layers; merged results table with classification, direction and consistency check; layer selector, *significant only*, *consistent only*. |
| **🌋 Plots** | Volcano (occupancy / site / protein), heatmap of top significant occupancy sites, depth and acquisition-order QC. Volcano and heatmap are built when *Draw* is clicked (see §6). |
| **💾 Export** | To the Results / Figures / limpa directories: combined results, per-layer full DE, occupancy and site matrices, sample metadata, full analysis RDS; SVG volcano, DPC curves, PCA / PLS-DA, vooma mean-variance trend; optional heatmap PNG. Browser download buttons for the main tables. |
| **📖 PSP Annotation** | Matches sites to the PhosphoSitePlus ubiquitination flat file (licence preamble auto-skipped; **human entries only**) by accession + lysine position. Adds `reported_in_PSP`, `LT_LIT`, `MS_LIT`, `MS_CST`, `SITE_GRP_ID`; lists significant novel sites; exports TSV. |
| **🧬 GSEA .rnk** | Gene-level rank files from the DE results (occupancy / site / protein). Methods: Stouffer-combined t (default), t-statistic, logFC, logFC × −log10 p. Gene symbol or protein group as identifier. Written to `<Results>/GSEA_rnk/` and downloadable. |
| **🧪 Proteome integration** | Gene-level occupancy (Stouffer z of site P-values with the sign of t, BH across genes) combined with a gene-level total-proteome DE table measured **without MG-132** (columns auto-detected). Classes: non-degradative candidate, same direction (strict / relaxed), opposite direction, intermediate, occupancy n.s. Occupancy comes from the session or a site-level DE file; several comparisons can be added and faceted. Optional E3-substrate literature table. Thresholds are live; tables export as CSV, figures as SVG, or all to `<Results>/Proteome_integration/`. The proteome logFC must have the same orientation as the occupancy contrast. |

---

## 4. Input requirements and conventions

**Files**
- DIA-NN parquet reports for the diGly and total-proteome experiments. Needed columns: `Run`, `Precursor.Normalised`, `Protein.Group`, `Genes`, `Modified.Sequence`, `Stripped.Sequence`, `Proteotypic`, `PTM.Q.Value`, `PTM.Site.Confidence`, and the q-value columns (`Q.Value`, `Lib.Q.Value`, `Lib.PG.Q.Value`, `PG.Q.Value`) where available.
- The FASTA used for the search. Every `Protein.Group` accession should exist in it.

**Run names** (the sample sheet is built from them)
- *Cell line*: matched by the token map. **Use upper-case tokens** (`ASPC=AsPC-1`, `CAPAN=Capan-1`, `PANC=PANC03.27`). Tokens are matched case-insensitively in the run name, but the lookup key is upper-cased, so mixed-case tokens such as `Capan=` do not resolve.
- *Replicate*, *well*, *injection order* are parsed with patterns built for names like `..._<rep>_R<X><n>_1_..._<injection>`. For other naming schemes these stay `NA`; set **Replicate** by hand in the table before loading (every selected run needs a unique cell line + replicate, and diGly and proteome runs are paired on that key).
- Exclude QC-failed runs by deselecting them in the run tables. Deselect **after** finishing any cell edits (see §6).

**Design**
- One grouping factor (`cell_line`), ≥ 2 replicates per line recommended. Injection order is stored but cannot be modelled if acquisition was blocked by line.

**Configuration must be applied.** Pipeline tabs read the *applied* configuration. Click **Apply Configuration** after loading data and again after any change. Applying before loading gives an empty target line and empty output paths.

---

## 5. Method summary

**Occupancy** (log2 space, at DPC-Quant level)

```
Occupancy_ij = log2(site_ij) − log2(protein_ij)
SE_occ       = sqrt(SE_site² + SE_prot²)
```

A relative measure of ubiquitination stoichiometry independent of protein abundance;
not absolute stoichiometry. `n.observations` of the difference is the elementwise
minimum of the two layers.

**DE layers**

| Layer | Input | Question |
|-------|-------|----------|
| Site | `q.site$E` | Is the site changing? (confounded by protein level) |
| Occupancy | `q.occ$E` | Is occupancy changing independent of protein? |
| Protein | `q.prot$E` | Is the parent protein changing? |

**Classification** (FDR < threshold at each layer): ubiquitination-specific (occ + site, protein n.s.), both layers change, occupancy only, abundance-driven, site change not occupancy, ns.

**Consistency**: target-minus-reference line means of occupancy must have the same sign against every reference line and a minimum absolute difference ≥ logFC threshold / 2 (trivial in *vs control* mode).

**Two different "minimum observation" filters**
- *Missingness filter* (Filtering): raw precursor intensities, before site mapping.
- *Min obs / line* (Configuration, applied in Occupancy): `n.observations` after DPC-Quant; a site passes if the target **or** any reference line has ≥ N observations.

**`treat()` vs `eBayes()`**: `treat` tests H₀ |logFC| ≤ threshold and is recommended for occupancy; `eBayes` tests logFC ≠ 0.

### Key defaults

| Parameter | Default | Where |
|-----------|---------|-------|
| Precursor q-value | 0.01 | Configuration |
| PTM q-value | 0.05 | Configuration |
| Site confidence | ≥ 0.75 | Configuration |
| Missingness filter | ≥ 1 valid value | Filtering |
| Min obs / line | 2 | Configuration |
| FDR / logFC | 0.05 / 1.0 | Configuration |
| `treat()` | on | Configuration |
| Normalisation | median (both layers) | Configuration |
| Drop C-terminal GG / ubiquitin genes | on / on | Configuration |
| Remove cRAP / Keep GFP | on / off | Filtering |

---

## 6. Known limitations

Found in a code review of this version (not yet fixed unless marked).

**Behaviour that can mislead**
- **Differences from the Rmds.** The Rmds apply q-value cutoffs inside `readDIANN` (per observation, all q-value columns). The app reads with `q.cutoffs = Inf` and filters afterwards on the minimum of the q-value columns stored per precursor (this avoids NA `Modified.Sequence`). Precursor and site counts may therefore differ slightly. Compare one target line against its Rmd before relying on the app for published numbers.
- **Run tables:** editing a cell re-renders the table and re-selects all rows, undoing any deselection (exclusions). Edit first, deselect last.
- **Plots tab:** volcano and heatmap keep the DE results from the moment *Draw* was clicked. Click *Draw* again after re-running DE.
- **Missingness mode "none"** keeps precursors that are NA in all loaded runs (the Rmds always drop them). The "0 = keep all NA" hint under the absolute mode is not honoured; the minimum is forced to 1.
- **Protein volcano** repeats each protein once per site (results are site-keyed); counts are inflated. Inherited from the Rmds.
- **Export filenames and titles** use `_vs_rest` even in *vs control* mode, and the heatmap export takes the top consistent sites regardless of significance.
- **`vs control` mode** fits the design on all loaded lines (variance pooled across them) while testing only target vs control; the reproducibility filter uses target and control only.

**Failures on non-standard data**
- A `Protein.Group` accession absent from the FASTA stops Filtering with a "subscript out of bounds" error (also in the Rmds). Use the search FASTA; keep cRAP removal on, or add custom contaminant prefixes.
- Missing replicate or injection order is not validated at *Load*; duplicates surface as a cryptic row-names error.
- The *Proteome run map CSV* field is currently not used.

**Scope**
- PSP annotation keeps only `ORGANISM == human`; the ubiquitin-gene filter uses human symbols (case-sensitive). For other organisms, disable *Drop ubiquitin genes* and treat PSP results as not applicable.
- PSP matching takes the accession from the first entry of `Protein.Group` rather than from `site_id`; sites mapped to a later accession in a multi-protein group can be reported as novel.
- Gene-level Stouffer statistics assume independent sites; sites of one protein share evidence, so gene-level FDRs are optimistic. Proteome integration uses `treat` P-values while the .rnk files use t-statistics.
- Sites on a doubly modified peptide share precursor evidence; `ambiguous_pos` and `is_prot_cterm` flags are carried in the per-layer DE tables but not in the combined results table.

---

## 7. Changelog

**v1 → v2 fixes** (carried over from the previous README)

| Module | Fix |
|--------|-----|
| `mod_config.R` | Duplicated file content removed |
| `mod_normalization.R` | `targets` used before assignment in `renderPlot` |
| `mod_data.R` | `readDIANN` with `q.cutoffs = Inf`; q-value filter moved to Filtering so annotation columns stay intact |
| `mod_filtering.R` | Precursor q-value filter added as an explicit step; per-pattern `stringr::fixed()` matching |
| `mod_export.R` | Spurious `contrasts =` argument removed from the `dpcDE` export call |
| `mod_plots.R` | SE-weighted MDS removed; `make_volcano` defined at top level for reuse by export |
| `hybrid_path_ui` | `icon()` replaced by an emoji for `shinyDirButton` |
