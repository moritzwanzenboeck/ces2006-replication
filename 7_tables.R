# 7_tables.R
# Generates all LaTeX tables for the CES 2006 replication paper.
# Outputs: output/tables/*.tex (booktabs format)
# Tables requiring dcc_fits.RData are skipped with a warning if missing.
# Run after: 2_returns_compute, 3_stage1_garch_fit, 4_nonparametric_tests,
#            5_bootstrap, and (for Tables 6b/7/8) 6_dcc_estimation.

rm(list = ls())
tryCatch(
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)),
  error = function(e) {
    setwd("C:\\Users\\morit\\OneDrive\\Uni\\Universität Zürich\\FS 2026\\Topics in Time Series Analysis\\Seminar Paper\\implementation")
  }
)

packages <- c("dplyr", "tidyr", "moments", "rugarch", "readxl")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# LOAD DATA
# -------------------------------------------------------------------------

load("data/weekly_returns_usd.RData")    # returns_mat, returns_df, asset_meta
load("data/stage1_fits.RData")           # stage1_best, stage1_selected, e_std, h_mat
load("data/nonparam_tests.RData")        # test_results, table4_summary

model_names <- c("DCC", "ADCC", "GDCC", "AGDCC")
fits <- setNames(vector("list", length(model_names)), model_names)
for (m in model_names) {
  cache <- sprintf("data/dcc_fit_%s_3stage_tmp.rds", m)
  if (file.exists(cache)) {
    fits[[m]] <- readRDS(cache)
    cat(sprintf("Loaded %s from 3-stage cache.\n", m))
  }
}
has_dcc <- !is.null(fits[["AGDCC"]])
if (!has_dcc)
  cat("NOTE: 3-stage VCV cache not found â€” Tables 6b, 7, 8 will be skipped.\n\n")

has_boot <- file.exists("data/bootstrap_corr.RData")
if (has_boot) {
  load("data/bootstrap_corr.RData")      # bootstrap_corr (obs_medians, p_values)
} else {
  cat("NOTE: data/bootstrap_corr.RData not found â€” bootstrap tests unavailable.\n\n")
}

T_obs   <- nrow(returns_mat)
k       <- ncol(returns_mat)
nms     <- colnames(returns_mat)
eq_nms  <- nms[startsWith(nms, "eq_")]
bd_nms  <- nms[startsWith(nms, "bd_")]
eq_idx  <- which(startsWith(nms, "eq_"))
bd_idx  <- which(startsWith(nms, "bd_"))

# Regional membership (update if asset set changes)
reg <- list(
  eq_australasia = c("eq_australia","eq_hong_kong","eq_japan",
                     "eq_new_zealand","eq_singapore"),
  eq_europe      = c("eq_austria","eq_belgium","eq_denmark","eq_france",
                     "eq_germany","eq_ireland","eq_italy","eq_netherlands",
                     "eq_norway","eq_spain","eq_sweden","eq_switzerland",
                     "eq_united_kingdom"),
  eq_namerica    = c("eq_canada","eq_mexico","eq_usa"),
  bd_australasia = c("bd_australia","bd_japan","bd_new_zealand"),
  bd_europe      = c("bd_austria","bd_belgium","bd_denmark","bd_france",
                     "bd_germany","bd_ireland","bd_italy","bd_netherlands",
                     "bd_norway","bd_spain","bd_sweden","bd_switzerland",
                     "bd_united_kingdom"),
  bd_namerica    = c("bd_canada","bd_united_states")
)
reg <- lapply(reg, intersect, nms)

# -------------------------------------------------------------------------
# HELPER: write booktabs LaTeX table to file
# -------------------------------------------------------------------------

save_tex <- function(lines, filename) {
  path <- file.path("output/tables", filename)
  writeLines(lines, path)
  cat("Saved:", path, "\n")
}

fmt  <- function(x, digits = 4) formatC(x, format = "f", digits = digits)
fmtp <- function(x) ifelse(x < 0.001, "<0.001", fmt(x, 3))
ch   <- function(x) paste0("\\multicolumn{1}{c}{", x, "}")

# Wrap multi-word labels in \makecell{Word1 \\ Word2} to save column width
make_col_label <- function(lbl) {
  words <- strsplit(lbl, " ")[[1]]
  if (length(words) == 1) return(lbl)
  paste0("\\makecell{", paste(words, collapse = " \\\\ "), "}")
}

booktabs_table <- function(header_line, data_lines, caption, label,
                           col_spec = NULL, notes = NULL) {

  nc    <- length(strsplit(header_line, " & ")[[1]])
  cspec <- if (!is.null(col_spec)) col_spec else paste0("l", strrep("S", nc - 1))

  cleaned_data <- sapply(data_lines, function(line) {
    trimmed <- trimws(line)
    if (trimmed == "" || grepl("\\\\\\\\$", trimmed) ||
        grepl("\\\\hline", trimmed) || grepl("\\\\midrule", trimmed) ||
        grepl("\\\\toprule", trimmed) || grepl("\\\\bottomrule", trimmed))
      return(line)
    paste0(line, " \\\\")
  }, USE.NAMES = FALSE)

  tabular <- c(
    paste0("\\begin{tabular}{", cspec, "}"),
    "\\toprule",
    paste0(header_line, " \\\\"),
    "\\midrule",
    cleaned_data,
    "\\bottomrule",
    "\\end{tabular}"
  )

  inner <- if (!is.null(notes))
    c("\\begin{threeparttable}", tabular,
      "\\begin{tablenotes}[flushleft]", "\\small",
      paste0("\\item ", notes),
      "\\end{tablenotes}", "\\end{threeparttable}")
  else tabular

  c(
    "\\begin{table}[H]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    inner,
    "\\end{table}"
  )
}

# Pretty series name: eq_united_states â†’ United States (Equity)
pretty_name <- function(nm) {
  cls  <- ifelse(startsWith(nm, "eq_"), "Equity", "Bond")
  base <- sub("^(eq|bd)_", "", nm)
  base <- gsub("_", " ", base)
  base <- tools::toTitleCase(base)
  paste0(base, " (", cls, ")")
}
bare_name <- function(nm) {
  b <- tools::toTitleCase(gsub("_", " ", sub("^(eq|bd)_", "", nm)))
  if (b == "Usa") "United States" else b
}

# -------------------------------------------------------------------------
# TABLE 1 â€” Descriptive Statistics
# -------------------------------------------------------------------------

cat("--- Table 1: Descriptive Statistics ---\n")

# Columns 1-2: annualised moments of raw returns
ann_mean <- colMeans(returns_mat) * 52 * 100
ann_sd   <- apply(returns_mat, 2, sd) * sqrt(52) * 100

# Columns 3-4: skewness and kurtosis of raw returns (raw kurtosis, normal = 3)
sk_raw <- apply(returns_mat, 2, moments::skewness)
ku_raw <- apply(returns_mat, 2, moments::kurtosis)

# Columns 5-6: skewness and kurtosis of GARCH standardised residuals e_it / sqrt(h_it)
sk_std <- apply(e_std, 2, moments::skewness)
ku_std <- apply(e_std, 2, moments::kurtosis)   # raw kurtosis (normal = 3)

# Jarque-Bera on standardised residuals: JB = T/6 * (S^2 + (K-3)^2/4) ~ chi^2(2)
jb_stat <- T_obs / 6 * (sk_std^2 + (ku_std - 3)^2 / 4)
jb_pval <- pchisq(jb_stat, df = 2, lower.tail = FALSE)

# Stars on kurtosis of standardised residuals (reverse convention: stars = appears normal)
# * = insignificant at 5% (p > 0.05); ** = insignificant at 1% (0.01 < p <= 0.05)
jb_star <- function(p) {
  if (is.na(p) || p <= 0.01) ""
  else if (p <= 0.05) "\\text{**}"
  else "\\text{*}"
}

make_t1_rows <- function(idxs) {
  sapply(nms[idxs], function(nm) {
    paste(
      paste0("\\hspace{1em}", pretty_name(nm)),
      fmt(ann_mean[nm], 2),
      fmt(ann_sd[nm],   2),
      fmt(sk_raw[nm],   2),
      fmt(ku_raw[nm],   2),
      fmt(sk_std[nm],   2),
      paste0(fmt(ku_std[nm], 2), jb_star(jb_pval[nm])),
      sep = " & "
    )
  })
}

hdr <- paste("Asset", ch("Mean (\\%)"), ch("SD (\\%)"), ch("Skew"), ch("Kurt"),
             ch("$\\hat{\\varepsilon}$ Skew"), ch("$\\hat{\\varepsilon}$ Kurt"), sep = " & ")
rows <- c(
  "\\multicolumn{7}{l}{\\textit{Equities (21 series)}} \\\\",
  make_t1_rows(eq_idx),
  "\\midrule",
  "\\multicolumn{7}{l}{\\textit{Bonds (18 series)}} \\\\",
  make_t1_rows(bd_idx)
)

save_tex(
  booktabs_table(hdr, rows,
    caption  = "Descriptive statistics",
    label    = "tab:descriptive",
    col_spec = "lSSSSSSS",
    notes    = paste0(
      "Weekly USD log returns, 2000--2025 ($T=1{,}355$, $k=39$). ",
      "Mean and SD are annualised ($\\times 52$ and $\\times\\sqrt{52}$, $\\times 100$). ",
      "Skewness and kurtosis are for raw returns; kurtosis is raw (normal $=3$). ",
      "$\\hat{\\varepsilon}$ columns: moments of Stage-1 BIC-selected GARCH standardised residuals. ",
      "Jarque--Bera: $JB=(T/6)(S^2+(K-3)^2/4)\\sim\\chi^2(2)$ under normality. ",
      "* and ** denote failure to reject normality at 5\\% and 1\\%, respectively.")),
  "table1_descriptive.tex"
)

# -------------------------------------------------------------------------
# TABLE A1 â€” Data sources and index identifiers
# -------------------------------------------------------------------------
# Built from the data_sources.xlsx provenance sheets (Datastream_indeces = bonds,
# MSCI = equities). Two-panel layout matching Table 1 (Equities / Bonds).
# Emitted as a longtable so it can break across pages; the enclosing landscape
# environment lives in main.tex.

cat("--- Table A1: Data sources ---\n")

ds_xlsx     <- "data_sources.xlsx"
has_sources <- file.exists(ds_xlsx)

if (!has_sources) {
  cat("NOTE: data_sources.xlsx missing â€” Table A1 skipped.\n\n")
} else {
  ds_bd <- as.data.frame(readxl::read_excel(ds_xlsx, sheet = "Datastream_indeces"))
  ds_eq <- as.data.frame(readxl::read_excel(ds_xlsx, sheet = "MSCI"))

  # Map xlsx country labels to the bd_*/eq_* asset names used in returns_mat
  bd_key <- function(country) {
    country <- gsub("U\\.K\\.", "United Kingdom", country)
    country <- gsub("U\\.S\\.", "United States", country)
    paste0("bd_", gsub("\\s+", "_", tolower(country)))
  }
  eq_key <- function(country) {
    country <- gsub("United States", "usa", country)   # dict uses eq_usa
    paste0("eq_", gsub("\\s+", "_", tolower(country)))
  }
  ds_bd$key <- bd_key(ds_bd$country)
  ds_eq$key <- eq_key(ds_eq$country)

  ds_row <- function(nm) {
    if (startsWith(nm, "bd_")) {
      r       <- ds_bd[match(nm, ds_bd$key), ]
      idxname <- "TR 10-Year Gov. Benchmark"
      mnem    <- r[["mnemonic_BM105Y"]]
      code    <- r[["investing.com"]]
      url     <- r[["investing.com_url"]]
    } else {
      r       <- ds_eq[match(nm, ds_eq$key), ]
      idxname <- r[["msci_imi_index_name"]]
      mnem    <- ""
      code    <- as.character(r[["msci_index_code"]])
      url     <- r[["msci.com_url"]]
    }
    # Clickable link; displayed text drops the https://www. prefix. Whether a
    # web link opens in a new browser tab is decided by the PDF viewer, not the
    # document: PDF /NewWindow is only valid for file-open actions, not URI (web)
    # links, so no hyperref option forces it. NOTE: on rotated pdflscape pages the
    # clickable rectangle can appear offset in pdf.js viewers (VS Code preview);
    # verify in a spec-compliant reader.
    url_disp <- sub("^https://www\\.", "", url)
    paste(paste0("\\hspace{1em}", bare_name(nm)), idxname, mnem, code,
          paste0("\\href{", url, "}{", url_disp, "}"), sep = " & ")
  }

  eq_rows_ds <- vapply(nms[eq_idx], ds_row, character(1))
  bd_rows_ds <- vapply(nms[bd_idx], ds_row, character(1))

  hdr_ds <- paste("Country", "Index name", "DS mnemonic", "Code",
                  "Source URL", sep = " & ")

  ds_notes <- paste0(
    "\\textit{Notes:} Index identifiers and provenance for the $k=39$ series. ",
    "Equities: MSCI IMI country indices (gross); Code is the MSCI index code. ",
    "Bonds: Thomson Reuters 10-Year Government Benchmark total-return indices; ",
    "Mnemonic is the Datastream code (\\texttt{BM\\textit{cc}10Y}) and Code the ",
    "investing.com identifier.")

  # 39 series is too tall for one landscape page, so use a longtable that breaks
  # across pages (repeating header). The enclosing landscape environment lives in
  # main.tex; this file stays orientation-agnostic.
  notes_row <- paste0(
    "\\multicolumn{5}{@{}p{\\linewidth}@{}}{\\footnotesize ", ds_notes, "} \\\\")

  ds_lines <- c(
    "{\\footnotesize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\setlength{\\LTcapwidth}{\\linewidth}",   # left-align caption (no indent)
    "\\begin{longtable}{l l l l l}",
    "\\caption{Data sources and index identifiers}\\label{tab:data_sources} \\\\",
    "\\toprule",
    paste0(hdr_ds, " \\\\"),
    "\\midrule",
    "\\endfirsthead",
    "\\multicolumn{5}{c}{\\tablename~\\thetable{} -- continued from previous page} \\\\",
    "\\toprule",
    paste0(hdr_ds, " \\\\"),
    "\\midrule",
    "\\endhead",
    "\\midrule",
    "\\multicolumn{5}{r}{\\footnotesize Continued on next page} \\\\",
    "\\endfoot",
    "\\bottomrule",
    notes_row,
    "\\endlastfoot",
    sprintf("\\multicolumn{5}{l}{\\textit{Equities (%d series)}} \\\\", length(eq_idx)),
    paste0(eq_rows_ds, " \\\\"),
    # Force the Bonds panel onto a fresh page
    "\\pagebreak",
    sprintf("\\multicolumn{5}{l}{\\textit{Bonds (%d series)}} \\\\", length(bd_idx)),
    paste0(bd_rows_ds, " \\\\"),
    "\\end{longtable}",
    "}"
  )

  save_tex(ds_lines, "table_a1_data_sources.tex")
}

# -------------------------------------------------------------------------
# TABLE 2 â€” Average Correlations
# -------------------------------------------------------------------------

cat("--- Table 2: Average Correlations ---\n")

C_mat <- cor(returns_mat)

# Unique pairwise correlations: upper triangle for symmetric (same) series groups
corr_pairs2 <- function(r, c) {
  r <- intersect(r, nms); c <- intersect(c, nms)
  if (length(r) == 0 || length(c) == 0) return(numeric(0))
  sub <- C_mat[r, c, drop = FALSE]
  if (identical(r, c)) sub[lower.tri(sub, diag = TRUE)] <- NA
  vals <- as.numeric(sub); vals[!is.na(vals)]
}

fmtc2 <- function(x) if (length(x) == 0 || is.na(x)) "{}" else fmt(x, 4)

reg_lbl2 <- c(
  eq_australasia = "Australasia", eq_europe = "Europe",   eq_namerica = "N.~America",
  bd_australasia = "Australasia", bd_europe = "Europe",   bd_namerica = "N.~America"
)

# 3x3 regional matrix rows (lower-tri for within-class, full for cross-class)
make_mat2 <- function(row_regs, col_regs, lower_tri = TRUE) {
  rn <- names(row_regs); cn <- names(col_regs)
  sapply(seq_along(rn), function(i) {
    vals <- sapply(seq_along(cn), function(j) {
      if (lower_tri && j > i) return("")
      fmtc2(mean(corr_pairs2(row_regs[[i]], col_regs[[j]])))
    })
    paste(c(paste0("\\hspace{1em}", reg_lbl2[rn[i]]), vals), collapse = " & ")
  })
}

eq_reg_list2 <- list(eq_australasia = reg$eq_australasia,
                     eq_europe      = reg$eq_europe,
                     eq_namerica    = reg$eq_namerica)
bd_reg_list2 <- list(bd_australasia = reg$bd_australasia,
                     bd_europe      = reg$bd_europe,
                     bd_namerica    = reg$bd_namerica)

# Within-class panel: aggregate row + 3x3 lower-triangular matrix
panel_within2 <- function(title, all_r, all_c, reg_list) {
  v   <- corr_pairs2(all_r, all_c)
  mr  <- make_mat2(reg_list, reg_list, lower_tri = TRUE)
  cn  <- names(reg_list)
  chdr <- paste(c("", sapply(cn, function(x) ch(reg_lbl2[x]))), collapse = " & ")
  c(
    paste0("\\multicolumn{4}{l}{\\textit{", title, "}} \\\\"),
    paste0(paste("", ch("Mean"), ch("Minimum"), ch("Maximum"), sep = " & "), " \\\\"),
    "\\cmidrule(l){2-4}",
    paste0(paste("\\hspace{1em}All pairs",
                 fmtc2(mean(v)), fmtc2(min(v)), fmtc2(max(v)), sep = " & "), " \\\\"),
    "\\midrule",
    paste0(chdr, " \\\\"),
    "\\cmidrule(l){2-4}",
    paste0(mr, " \\\\")
  )
}

# Cross-asset panel: equities in rows, bonds in columns
panel_cross2 <- function() {
  v  <- corr_pairs2(eq_nms, bd_nms)
  mr <- make_mat2(eq_reg_list2, bd_reg_list2, lower_tri = FALSE)
  bd_names <- c("bd_australasia", "bd_europe", "bd_namerica")
  c(
    "\\multicolumn{4}{l}{\\textit{Bond and equity indices}} \\\\",
    paste0(paste("", ch("Mean"), ch("Minimum"), ch("Maximum"), sep = " & "), " \\\\"),
    "\\cmidrule(l){2-4}",
    paste0(paste("\\hspace{1em}All pairs",
                 fmtc2(mean(v)), fmtc2(min(v)), fmtc2(max(v)), sep = " & "), " \\\\"),
    "\\midrule",
    paste0(" & \\multicolumn{3}{c}{Bonds} \\\\"),
    "\\cmidrule(lr){2-4}",
    paste0(paste(c("Equities",
                   sapply(bd_names, function(x) ch(reg_lbl2[x]))),
                 collapse = " & "), " \\\\"),
    "\\cmidrule(l){2-4}",
    paste0(mr, " \\\\")
  )
}

t2_lines <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\small",
  "\\caption{Average unconditional correlations}",
  "\\label{tab:avg_corr}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lSSS}",
  "\\toprule",
  panel_within2("Equity indices", eq_nms, eq_nms, eq_reg_list2),
  "\\midrule",
  panel_within2("Bond indices",   bd_nms, bd_nms, bd_reg_list2),
  "\\midrule",
  panel_cross2(),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  "\\small",
  paste0("\\item Weekly USD log returns, 2000--2025. ",
         "\\textit{All pairs}: mean, minimum, and maximum across all unique ",
         "pairwise correlations in the group. ",
         "Regional matrices: average correlation within each region (diagonal) ",
         "and across regions (off-diagonal); lower-triangular for within-class panels. ",
         "Australasia equity: AU, HK, JP, NZ, SG. Australasia bond: AU, JP, NZ."),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)

save_tex(t2_lines, "table2_avg_correlations.tex")

# -------------------------------------------------------------------------
# TABLES 3a / 3b / 3c â€” Correlation Matrices
# -------------------------------------------------------------------------

cat("--- Tables 3a/b/c: Correlation Matrices ---\n")

write_corr_matrix <- function(row_nms, col_nms, filename, caption, label,
                              notes = NULL) {
  row_nms   <- intersect(row_nms, nms)
  col_nms   <- intersect(col_nms, nms)
  sub       <- C_mat[row_nms, col_nms, drop = FALSE]
  symmetric <- identical(row_nms, col_nms)

  # For lower-triangular symmetric tables, drop the first row (all blank)
  # and last column (all blank) to avoid empty margins.
  if (symmetric) {
    row_nms <- row_nms[-1]
    col_nms <- col_nms[-length(col_nms)]
    sub     <- sub[row_nms, col_nms, drop = FALSE]
  }

  col_labels <- sub("^(eq|bd)_", "", col_nms)
  col_labels <- tools::toTitleCase(gsub("_", " ", col_labels))
  col_labels <- sapply(col_labels, make_col_label)
  hdr_line   <- paste(c("", sapply(col_labels, ch)), collapse = " & ")
  nc         <- length(col_nms)
  cspec      <- paste0("l", strrep("S", nc))

  data_lines <- character(length(row_nms))
  for (ii in seq_along(row_nms)) {
    vals <- fmt(sub[ii, ], 3)
    if (symmetric && ii < nc) vals[(ii + 1):nc] <- "{}"  # lower triangular
    row_label       <- sub("^(eq|bd)_", "", row_nms[ii])
    row_label       <- tools::toTitleCase(gsub("_", " ", row_label))
    data_lines[ii]  <- paste(c(row_label, vals), collapse = " & ")
  }

  tabular <- c(
    paste0("\\begin{tabular}{", cspec, "}"),
    "\\toprule",
    paste0(hdr_line, " \\\\"),
    "\\midrule",
    paste0(data_lines, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}"
  )
  inner <- if (!is.null(notes))
    c("\\begin{threeparttable}", tabular,
      "\\begin{tablenotes}[flushleft]", "\\small",
      paste0("\\item ", notes),
      "\\end{tablenotes}", "\\end{threeparttable}")
  else tabular

  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\small",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{adjustbox}{max width=1.5\\textheight}",
    inner,
    "\\end{adjustbox}",
    "\\end{table}"
  )
  save_tex(lines, filename)
}

t3_note <- "Unconditional correlations from weekly USD log returns, 2000--2025 ($T=1{,}355$)."
write_corr_matrix(eq_nms, eq_nms,
  "table3a_corr_equity.tex",
  "Equity return correlations",
  "tab:corr_equity",
  notes = paste(t3_note, "Lower triangular; diagonal omitted."))

write_corr_matrix(eq_nms, bd_nms,
  "table3b_corr_eq_bond.tex",
  "Equity--bond return correlations",
  "tab:corr_eq_bond",
  notes = t3_note)

write_corr_matrix(bd_nms, bd_nms,
  "table3c_corr_bond.tex",
  "Bond return correlations",
  "tab:corr_bond",
  notes = paste(t3_note, "Lower triangular; diagonal omitted."))

# -------------------------------------------------------------------------
# TABLE 4 â€” Conditional Partial Covariances
# -------------------------------------------------------------------------

cat("--- Table 4: Conditional Partial Covariances ---\n")

# Build series â†’ region lookup
series_region <- setNames(rep(NA_character_, length(nms)), nms)
for (rn in names(reg)) series_region[intersect(reg[[rn]], nms)] <- rn

tr4       <- test_results
tr4$reg_i <- series_region[tr4$series_i]
tr4$reg_j <- series_region[tr4$series_j]

get_pairs4 <- function(r1, r2) {
  tr4[(tr4$reg_i == r1 & tr4$reg_j == r2) |
      (tr4$reg_i == r2 & tr4$reg_j == r1), ]
}

prop4 <- function(df, alpha) {
  if (is.null(df) || nrow(df) == 0) return(NA_real_)
  mean(df$p_reg < alpha, na.rm = TRUE)
}
fmtp4 <- function(x) if (is.na(x)) "{--}" else fmt(x, 3)

all_pairs4   <- tr4
eq_pairs4    <- tr4[tr4$type == "intra_equity", ]
bd_pairs4    <- tr4[tr4$type == "intra_bond",   ]
inter_pairs4 <- tr4[tr4$type == "inter",         ]

eq_regs4 <- c("eq_australasia", "eq_europe", "eq_namerica")
bd_regs4 <- c("bd_australasia", "bd_europe", "bd_namerica")
rlabel4  <- c(
  eq_australasia = "Australasia", eq_europe = "Europe",   eq_namerica = "N.~America",
  bd_australasia = "Australasia", bd_europe = "Europe",   bd_namerica = "N.~America"
)

# Matrix rows for a regional section; proportions at 10%; 5-col tabular (pad 4th val to "")
make_mat4_rows <- function(row_regs, col_regs, lower_tri = TRUE, bond_label_first_col = FALSE) {
  sapply(seq_along(row_regs), function(i) {
    vals <- sapply(seq_along(col_regs), function(j) {
      if (lower_tri && j > i) return("{}")
      fmtp4(prop4(get_pairs4(row_regs[i], col_regs[j]), 0.10))
    })
    rlbl <- paste0("\\multicolumn{1}{l}{\\hspace{1em}", rlabel4[row_regs[i]], "}")
    if (bond_label_first_col && i == 2) {
      paste(c("Bonds", rlbl, vals), collapse = " & ")
    } else {
      paste(c("", rlbl, vals), collapse = " & ")
    }
  })
}

# Column header row for a matrix section (3 region labels + empty 4th)
mat_col_hdr4 <- function(col_regs) {
  paste(c("", sapply(col_regs, function(r) ch(rlabel4[r]))), collapse = " & ")
}

# Within-class section (Stocks only / Bonds only): lower-triangular 3x3
make_section4 <- function(title, row_regs, col_regs) {
  c(
    "\\cmidrule(l){2-5}",
    paste0("& \\multicolumn{4}{l}{\\textit{", title, "}} \\\\"),
    "\\cmidrule(l){2-5}",
    paste0("& ", mat_col_hdr4(col_regs), " \\\\"),
    "\\cmidrule(l){3-5}",
    paste0(make_mat4_rows(row_regs, col_regs, lower_tri = TRUE), " \\\\")
  )
}

# Cross-asset section: bonds in rows, stocks in columns
make_cross4 <- function() {
  bonds_sub  <- paste0(c("& & ", paste(sapply(eq_regs4, function(r) ch(rlabel4[r])), collapse = " & "), "\\\\"), " ")

  c(
    "\\cmidrule(l){2-5}",
    "& \\multicolumn{1}{l}{\\textit{Across stocks -- bonds}} & \\multicolumn{3}{c}{Stocks} \\\\",
    "\\cmidrule(l){2-5}",
    bonds_sub,
    "\\cmidrule(l){3-5}",
    paste0(make_mat4_rows(bd_regs4, eq_regs4, lower_tri = FALSE, bond_label_first_col = TRUE), " \\\\")
  )
}

t4_lines <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\small",
  "\\caption{Covariance asymmetry tests}",
  "\\label{tab:partial_cov}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lSSSS}",
  "\\toprule",
  paste0(paste("", ch("Overall"), ch("Intrastock"), ch("Intrabond"), ch("Interstock-bond"),
               sep = " & "), " \\\\"),
  "\\midrule",
  paste0(paste("Significant at 10\\%",
               fmtp4(prop4(all_pairs4, 0.10)), fmtp4(prop4(eq_pairs4, 0.10)),
               fmtp4(prop4(bd_pairs4, 0.10)), fmtp4(prop4(inter_pairs4, 0.10)),
               sep = " & "), " \\\\"),
  paste0(paste("Significant at 20\\%",
               fmtp4(prop4(all_pairs4, 0.20)), fmtp4(prop4(eq_pairs4, 0.20)),
               fmtp4(prop4(bd_pairs4, 0.20)), fmtp4(prop4(inter_pairs4, 0.20)),
               sep = " & "), " \\\\"),
  make_section4("Stocks only", eq_regs4, eq_regs4),
  make_section4("Bonds only",  bd_regs4, bd_regs4),
  make_cross4(),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  "\\small",
  paste0("\\item Proportion of the $C(39,2)=741$ asset pairs rejecting $H_0\\colon\\beta_{--}=\\beta_{++}$ ",
         "in the regression $\\varepsilon_{it}\\varepsilon_{jt}=\\alpha+\\beta_{--}I_{--,t}+\\beta_{++}I_{++,t}+u_t$, ",
         "where $I_{--,t}$ ($I_{++,t}$) indicates joint negative (positive) shocks. ",
         "Newey--West HAC standard errors. ",
         "Top panel: proportions by group at 10\\% and 20\\% levels. ",
         "Regional panels: 10\\% level; lower-triangular for within-class, full matrix for cross-class."),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)

save_tex(t4_lines, "table4_partial_cov.tex")

# -------------------------------------------------------------------------
# TABLE 5 â€” Univariate GARCH Parameters
# -------------------------------------------------------------------------

cat("--- Table 5: Univariate GARCH Parameters ---\n")

model_label <- c(
  "sGARCH"  = "GARCH",   "gjrGARCH" = "GJR-GARCH",
  "eGARCH"  = "EGARCH",  "apARCH"   = "APARCH",
  "fGARCH"  = "fGARCH"   # submodel appended below
)
fgarch_sub <- c(
  "AVGARCH" = "AVGARCH", "NGARCH" = "NARCH",
  "TGARCH"  = "ZARCH",   "NAGARCH" = "NAGARCH"
)

get_garch_row <- function(nm) {
  f <- stage1_best[[nm]]
  if (is.null(f)) return(paste(
    c(paste0("\\hspace{1em}", bare_name(nm)), "", rep("{}", 4)),
    collapse = " & "
  ))

  vmod <- f@model$modeldesc$vmodel
  vsub <- f@model$modeldesc$vsubmodel
  cf   <- coef(f)

  # Safe coefficient extraction (NA if absent)
  p <- function(n) { v <- cf[n]; if (length(v) == 0 || is.na(v)) NA_real_ else as.numeric(v) }
  f5 <- function(x) if (is.null(x) || is.na(x)) "{}" else fmt(x, 4)

  om   <- p("omega");  a1 <- p("alpha1"); b1 <- p("beta1")
  gam1 <- p("gamma1"); e11 <- p("eta11"); e21 <- p("eta21"); delt <- p("delta")

  if (vmod == "sGARCH") {
    mlab <- "GARCH"
    om_d <- om;  a_d <- a1;  g_d <- NA_real_;  b_d <- b1

  } else if (vmod == "fGARCH" && vsub == "AVGARCH") {
    # AVGARCH: no asymmetry (eta11=0 forced); Î±_cap = Î±_ru
    mlab <- "AVGARCH"
    om_d <- om;  a_d <- a1;  g_d <- NA_real_;  b_d <- b1

  } else if (vmod == "fGARCH" && vsub == "NGARCH") {
    # NARCH: no asymmetry; Î´ (power) shown in Î³/Î´ column
    mlab <- "NARCH"
    om_d <- om;  a_d <- a1;  g_d <- delt;  b_d <- b1

  } else if (vmod == "fGARCH" && vsub == "TGARCH") {
    # ZARCH: Î±_cap = Î±_ru*(1âˆ’Î·11), Î³_cap = 2*Î±_ru*Î·11
    mlab <- "ZARCH"
    if (is.na(e11)) e11 <- 0
    om_d <- om;  a_d <- a1 * (1 - e11);  g_d <- 2 * a1 * e11;  b_d <- b1

  } else if (vmod == "fGARCH" && vsub == "NAGARCH") {
    # NAGARCH: Î³_cap = âˆ’Î·21
    mlab <- "NAGARCH"
    om_d <- om;  a_d <- a1;  g_d <- -e21;  b_d <- b1

  } else if (vmod == "gjrGARCH") {
    # GJR-GARCH: trivially identical
    mlab <- "GJR-GARCH"
    om_d <- om;  a_d <- a1;  g_d <- gam1;  b_d <- b1

  } else if (vmod == "eGARCH") {
    # EGARCH: Î± and Î³ swapped vs Cappiello; Ï‰_cap = Ï‰_ru âˆ’ Î³_ruÂ·E|z|
    mlab <- "EGARCH"
    om_d <- om - gam1 * sqrt(2 / pi)  # E|z| = sqrt(2/pi) for N(0,1)
    a_d  <- gam1   # Î±_cap = Î³_ru (symmetric absolute term)
    g_d  <- a1     # Î³_cap = Î±_ru (asymmetric signed term)
    b_d  <- b1

  } else if (vmod == "apARCH") {
    # APARCH: Î±_cap = Î±_ru*(1âˆ’Î³_ru)^Î´, Î³_cap = Î±_ru*(1+Î³_ru)^Î´ âˆ’ Î±_cap
    mlab <- "APARCH"
    if (is.na(delt)) delt <- 2
    a_cap <- a1 * (1 - gam1)^delt
    g_cap <- a1 * (1 + gam1)^delt - a_cap
    om_d  <- om;  a_d <- a_cap;  g_d <- g_cap;  b_d <- b1

  } else {
    mlab <- if (!is.na(vsub) && vsub != "") paste0(vmod, "(", vsub, ")") else vmod
    om_d <- om;  a_d <- a1;  g_d <- NA_real_;  b_d <- b1
  }

  paste(paste0("\\hspace{1em}", bare_name(nm)), mlab, f5(om_d), f5(a_d), f5(g_d), f5(b_d), sep = " & ")
}

`%||%` <- function(a, b) if (!is.na(a) && !is.null(a)) a else b

make_t5_rows <- function(idxs) sapply(nms[idxs], get_garch_row)

hdr5 <- paste("Asset", "Model", ch("$\\omega$"), ch("$\\alpha$"), ch("$\\lambda$ or $\\gamma$"), ch("$\\beta$"), sep = " & ")
rows5 <- c(
  "\\multicolumn{6}{l}{\\textit{Equities:}} \\\\",
  make_t5_rows(eq_idx),
  "\\midrule",
  "\\multicolumn{6}{l}{\\textit{Bonds:}} \\\\",
  make_t5_rows(bd_idx)
)

save_tex(
  booktabs_table(hdr5, rows5,
    caption  = "Univariate GARCH estimates",
    label    = "tab:garch_params",
    col_spec = "llSSSS",
    notes    = paste0(
      "Stage-1 BIC-selected GARCH model and parameter estimates for weekly USD log returns. ",
      "$\\lambda$ or $\\gamma$: power or asymmetry parameter (model-dependent); blank if absent. ",
      "Parameters mapped to CES~(2006) notation: ",
      "EGARCH --- $\\alpha$ and $\\gamma$ exchanged relative to rugarch; ",
      "$\\omega_{\\text{cap}}=\\omega_{\\text{ru}}-\\gamma_{\\text{ru}}\\sqrt{2/\\pi}$. ",
      "ZARCH --- $\\alpha_{\\text{cap}}=\\alpha(1-\\eta_{11})$, $\\gamma_{\\text{cap}}=2\\alpha\\eta_{11}$. ",
      "APARCH --- $\\alpha_{\\text{cap}}=\\alpha(1-\\gamma)^\\delta$, ",
      "$\\gamma_{\\text{cap}}=\\alpha(1+\\gamma)^\\delta-\\alpha_{\\text{cap}}$.")),
  "table5_garch_params.tex"
)

# -------------------------------------------------------------------------
# TABLE 6a â€” Diagonal DCC Parameters (requires dcc_fits.RData)
# -------------------------------------------------------------------------

if (!has_dcc) {
  cat("Skipping Table 6a (dcc_fits.RData not found).\n")
} else {
  cat("--- Table 6a: Diagonal DCC Parameters ---\n")

  gdcc  <- fits[["GDCC"]]
  agdcc <- fits[["AGDCC"]]
  dcc   <- fits[["DCC"]]
  adcc  <- fits[["ADCC"]]

  # Asset label: indented "Australia stocks" / "Austria bonds"
  t6a_name <- function(nm) {
    paste0("\\hspace{1em}", bare_name(nm), if (startsWith(nm, "eq_")) " stocks" else " bonds")
  }

  # Display est^2; * if p > 0.05 (insignificant at 5%)
  fmt6a <- function(fit, pname) {
    idx <- match(pname, names(fit$psi_D))
    if (is.na(idx)) return("{--}")
    est  <- fit$psi_D[[idx]]
    pv   <- fit$pval[[idx]]
    star <- if (!is.na(pv) && pv > 0.05) "\\text{*}" else ""
    paste0(fmt(est^2, 4), star)
  }

  get_t6a_row <- function(nm) {
    pa <- paste0("[", nm, "].a1"); pg <- paste0("[", nm, "].g1"); pb <- paste0("[", nm, "].b1")
    paste(
      t6a_name(nm),
      fmt6a(gdcc,  pa),
      fmt6a(gdcc,  pb),
      fmt6a(agdcc, pa),
      fmt6a(agdcc, pg),
      fmt6a(agdcc, pb),
      sep = " & "
    )
  }

  scalar_row6a <- paste(
    "Scalar model",
    fmt6a(dcc,  "a"),
    fmt6a(dcc,  "b"),
    fmt6a(adcc, "a"),
    fmt6a(adcc, "g"),
    fmt6a(adcc, "b"),
    sep = " & "
  )

  t6a_lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\small",
    "\\caption{Diagonal DCC parameter estimates}",
    "\\label{tab:dcc_params}",
    "\\begin{threeparttable}",
    "\\begin{tabular}{lSSSSS}",
    "\\toprule",
    paste0(paste("", "\\multicolumn{2}{c}{Symmetric model}",
                     "\\multicolumn{3}{c}{Asymmetric model}", sep = " & "), " \\\\"),
    "\\cmidrule(lr){2-3}\\cmidrule(lr){4-6}",
    paste0(paste("", ch("$a_i^2$"), ch("$b_i^2$"),
                     ch("$a_i^2$"), ch("$g_i^2$"), ch("$b_i^2$"), sep = " & "), " \\\\"),
    "\\midrule",
    "\\multicolumn{6}{l}{\\textit{Equities:}} \\\\",
    paste0(sapply(nms[eq_idx], get_t6a_row), " \\\\"),
    "\\midrule",
    "\\multicolumn{6}{l}{\\textit{Bonds:}} \\\\",
    paste0(sapply(nms[bd_idx], get_t6a_row), " \\\\"),
    "\\midrule",
    paste0(scalar_row6a, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}[flushleft]",
    "\\small",
    paste0("\\item Stage-3 diagonal DCC and ADCC parameter estimates. ",
           "$a_i^2$, $g_i^2$, $b_i^2$: squared values entering the $Q_t$ recursion. ",
           "Standard errors from the full three-stage sandwich (Engle \\& Sheppard, 2001). ",
           "* denotes insignificance at the 5\\% level ($|z|<1.96$). ",
           "Bottom row: scalar model estimates."),
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{table}"
  )

  save_tex(t6a_lines, "table6a_dcc_params.tex")
}

# -------------------------------------------------------------------------
# TABLE 6b â€” DCC Model Comparison (requires dcc_fits.RData)
# -------------------------------------------------------------------------

if (!has_dcc) {
  cat("Skipping Table 6b (dcc_fits.RData not found).\n")
} else {
  cat("--- Table 6b: DCC Model Comparison ---\n")

  np_map <- c(DCC = 2L, ADCC = 3L, GDCC = 2L * k, AGDCC = 3L * k)

  t6_rows <- character(length(fits))
  for (i in seq_along(fits)) {
    m   <- names(fits)[i]
    f   <- fits[[m]]
    if (is.null(f)) { t6_rows[i] <- paste(c(m, rep("{--}", 3)), collapse=" & "); next }
    np  <- np_map[m]
    ll  <- f$ll
    bic <- -2 * ll + np * log(T_obs)
    t6_rows[i] <- paste(m, np, fmt(ll, 1), fmt(bic, 3), sep=" & ")
  }

  save_tex(
    booktabs_table(
      paste("Model", ch("$np$"), ch("LL"), ch("BIC"), sep = " & "),
      t6_rows,
      caption = "DCC model comparison",
      label   = "tab:dcc_comparison",
      notes   = paste0("LL: Stage-3 quasi-log-likelihood. ",
                       "BIC: $np\\cdot\\ln T - 2\\,\\text{LL}$, $T=1{,}355$. ",
                       "All LR tests vs.\\ DCC reject at $p<0.001$.")
    ),
    "table6b_dcc_comparison.tex"
  )
}

# -------------------------------------------------------------------------
# TABLE 7 â€” Equity Varianceâ€“Correlation Relationship
# -------------------------------------------------------------------------

if (!has_dcc) {
  cat("Skipping Table 7 (dcc_fits.RData not found).\n")
} else {
  cat("--- Table 7: Equity Variance-Correlation Relationship ---\n")

  Rt  <- fits[["AGDCC"]]$Rt
  Hu  <- fits[["AGDCC"]]$H_univ   # T Ã— k

  t7_rows <- sapply(eq_idx, function(i) {
    vi    <- sqrt(Hu[, i])
    j_set <- setdiff(eq_idx, i)
    col1  <- mean(sapply(j_set, function(j) cor(vi, Rt[i, j, ])), na.rm = TRUE)
    col2  <- cor(vi, colMeans(Rt[i, j_set, ]), use = "complete.obs")
    paste(bare_name(nms[i]), fmt(col1, 4), fmt(col2, 4), sep = " & ")
  })

  save_tex(
    booktabs_table(
      paste("Asset", ch("Mean pairwise corr."), ch("Corr.\\ w/ avg."), sep = " & "),
      t7_rows,
      caption = "Equity variance--correlation relationship",
      label   = "tab:equity_var_corr",
      notes   = paste0("Correlation between each equity series' conditional standard deviation ",
                       "$\\sqrt{\\hat{h}_{it}}$ and AG-DCC conditional correlations $\\hat{R}_{ij,t}$. ",
                       "Col.~1: mean of $\\text{cor}(\\sqrt{\\hat{h}_{it}},\\hat{R}_{ij,t})$ over all $j\\neq i$. ",
                       "Col.~2: correlation with the equally-weighted average correlation.")
    ),
    "table7_equity_var_corr.tex"
  )
}

# -------------------------------------------------------------------------
# TABLE 8 â€” Bond Volatilityâ€“Correlation Relationship
# -------------------------------------------------------------------------

if (!has_dcc) {
  cat("Skipping Table 8 (dcc_fits.RData not found).\n")
} else {
  cat("--- Table 8: Bond Volatility-Correlation Relationship ---\n")

  Rt  <- fits[["AGDCC"]]$Rt
  Hu  <- fits[["AGDCC"]]$H_univ

  t8_rows <- sapply(bd_idx, function(i) {
    vi    <- sqrt(Hu[, i])
    j_set <- setdiff(bd_idx, i)
    col1  <- mean(sapply(j_set, function(j) cor(vi, Rt[i, j, ])), na.rm = TRUE)
    col2  <- cor(vi, colMeans(Rt[i, j_set, ]), use = "complete.obs")
    paste(bare_name(nms[i]), fmt(col1, 4), fmt(col2, 4), sep = " & ")
  })

  save_tex(
    booktabs_table(
      paste("Asset", ch("Mean pairwise corr."), ch("Corr.\\ w/ avg."), sep = " & "),
      t8_rows,
      caption = "Bond variance--correlation relationship",
      label   = "tab:bond_var_corr",
      notes   = paste0("Correlation between each bond series' conditional standard deviation ",
                       "$\\sqrt{\\hat{h}_{it}}$ and AG-DCC conditional correlations $\\hat{R}_{ij,t}$. ",
                       "Col.~1: mean of $\\text{cor}(\\sqrt{\\hat{h}_{it}},\\hat{R}_{ij,t})$ over all $j\\neq i$. ",
                       "Col.~2: correlation with the equally-weighted average correlation.")
    ),
    "table8_bond_var_corr.tex"
  )
}

cat("\nDone. Tables written to output/tables/\n")
