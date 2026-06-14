# 8_charts.R
# Generates all figures for the CES 2006 replication paper.
# Outputs: output/charts/*.pdf
# Figures requiring dcc_fits.RData are skipped with a warning if missing.
# Run after: 3_stage1_garch_fit and (for Figs 2/4-10) 6_dcc_estimation.

rm(list = ls())
tryCatch(
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)),
  error = function(e) {
    setwd("C:\\Users\\morit\\OneDrive\\Uni\\Universität Zürich\\FS 2026\\Topics in Time Series Analysis\\Seminar Paper\\implementation")
  }
)

# plyr is listed before dplyr so dplyr's verbs (mutate/summarise) mask plyr's,
# not the other way round; plyr is only used qualified as plyr::round_any().
packages <- c("plyr", "ggplot2", "patchwork", "scales", "dplyr", "tidyr", "rugarch", "gridExtra", "purrr", "zoo")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/charts", recursive = TRUE, showWarnings = FALSE)

# Set TRUE to also write .tex TikZ files alongside each PDF.
# Requires the tikzDevice package and a working LaTeX installation.
USE_TIKZ <- TRUE
HAS_TIKZ <- USE_TIKZ && requireNamespace("tikzDevice", quietly = TRUE)
if (USE_TIKZ && !HAS_TIKZ) {
  message("tikzDevice not installed — install.packages('tikzDevice') to enable TikZ export.")
}
if (HAS_TIKZ) {
  options(tikzMetricsDictionary = NULL)  # discard cached font metrics; force full rerender
}

# -------------------------------------------------------------------------
# LOAD DATA
# -------------------------------------------------------------------------

load("data/weekly_returns_usd.RData")   # returns_mat, returns_df, asset_meta
load("data/stage1_fits.RData")          # stage1_best, h_mat, e_std

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
  cat("NOTE: 3-stage VCV cache not found — Figures 2, 3, 6-10 will be skipped.\n\n")

dates   <- as.Date(rownames(returns_mat))
T_obs   <- nrow(returns_mat)
nms     <- colnames(returns_mat)

# Fixed x-axis breaks: 2000, 2005, 2010, … through end of sample
x_breaks <- seq(as.Date("2000-01-01"), max(dates) + 365, by = "5 years")
eq_idx  <- which(startsWith(nms, "eq_"))
bd_idx  <- which(startsWith(nms, "bd_"))

reg <- list(
  eq_australasia = intersect(c("eq_australia","eq_hong_kong","eq_japan",
                                "eq_new_zealand","eq_singapore"), nms),
  eq_europe      = intersect(c("eq_austria","eq_belgium","eq_denmark","eq_france",
                                "eq_germany","eq_ireland","eq_italy","eq_netherlands",
                                "eq_norway","eq_spain","eq_sweden","eq_switzerland",
                                "eq_united_kingdom"), nms),
  eq_namerica    = intersect(c("eq_canada","eq_mexico","eq_usa"), nms),
  bd_emu         = intersect(c("bd_austria","bd_belgium","bd_france","bd_germany",
                                "bd_ireland","bd_italy","bd_netherlands","bd_spain"), nms),
  bd_europe_nemu = intersect(c("bd_denmark","bd_norway","bd_sweden",
                                "bd_switzerland","bd_united_kingdom"), nms),
  bd_namerica    = intersect(c("bd_canada","bd_united_states"), nms),
  bd_asia        = intersect(c("bd_australia","bd_japan","bd_new_zealand"), nms),
  eq_emu         = intersect(c("eq_austria","eq_belgium","eq_france","eq_germany",
                                "eq_ireland","eq_italy","eq_netherlands","eq_spain"), nms),
  eq_europe_nemu = intersect(c("eq_denmark","eq_norway","eq_sweden",
                                "eq_switzerland","eq_united_kingdom"), nms)
)

# Major events for overlays
events <- data.frame(
  date  = as.Date(c("2001-09-14","2008-09-15","2011-07-22","2020-03-19","2022-02-24")),
  label = c("9/11", "GFC", "Euro\nDebt", "COVID", "Ukraine")
)
events <- events[events$date >= min(dates) & events$date <= max(dates), ]

save_pdf <- function(p, filename, width = 8, height = 5, p_tikz = NULL) {
  path <- file.path("output/charts", filename)
  if (inherits(p, "gg") || inherits(p, "patchwork")) {
    ggsave(path, plot = p, width = width, height = height, device = "pdf")
  } else {
    pdf(path, width = width, height = height)
    p()
    dev.off()
  }
  cat("Saved:", path, "\n")

  if (HAS_TIKZ) {
    p_out <- if (!is.null(p_tikz)) p_tikz else p
    tex_path <- file.path("output/charts", sub("\\.pdf$", ".tex", filename))
    tw <- min(width, 6); th <- height * tw / width
    tryCatch({
      tikzDevice::tikz(tex_path, width = tw, height = th, standAlone = FALSE)
      if (inherits(p_out, "gg") || inherits(p_out, "patchwork")) print(p_out) else p_out()
      dev.off()
      cat("Saved TikZ:", tex_path, "\n")
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      message("TikZ export failed for ", filename, ": ", e$message)
    })
  }
}

theme_ces <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position  = "bottom")

# Figs 4-10: titles size 7, all other text base_size 9
theme_ces_panel <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position  = "bottom",
        plot.title       = element_text(size = 7),
        axis.title.y     = element_blank())

# Pretty label helper
plabel <- function(nm) {
  base <- sub("^(eq|bd)_", "", nm)
  base <- tools::toTitleCase(gsub("_", " ", base))
  if (base == "Usa") "United States" else base
}

# -------------------------------------------------------------------------
# FIGURE 1 — Volatility News Impact Curves
# -------------------------------------------------------------------------

cat("--- Figure 1: Volatility NIC ---\n")

nic_series <- list(
  equity = c("eq_usa", "eq_france", "eq_sweden"),
  bond   = c("bd_canada", "bd_switzerland")
)

build_nic_df <- function(snames) {
  purrr::map_dfr(snames, function(nm) {
    f  <- stage1_best[[nm]]
    if (is.null(f)) return(NULL)
    ni <- tryCatch(newsimpact(f), error = function(e) NULL)
    if (is.null(ni)) return(NULL)
    data.frame(series = plabel(nm), epsilon = ni$zx, h = ni$zy,
               stringsAsFactors = FALSE)
  })
}

nic_eq  <- build_nic_df(intersect(nic_series$equity, nms))
nic_bd  <- build_nic_df(intersect(nic_series$bond,   nms))

nic_pal  <- c("#D73027", "#4575B4", "#1A9641", "#FC8D59", "#762A83")
nic_ltys <- c("solid", "dashed", "dotted", "longdash", "twodash")

make_nic_panel <- function(df, title) {
  snames <- unique(df$series)
  cols <- setNames(nic_pal[seq_along(snames)],  snames)
  ltys <- setNames(nic_ltys[seq_along(snames)], snames)
  ggplot(df, aes(epsilon, h, colour = series, linetype = series)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = cols) +
    scale_linetype_manual(values = ltys) +
    labs(x = expression(epsilon[t-1]), y = expression(hat(h)[t]),
         title = title, colour = NULL, linetype = NULL) +
    theme_ces
}

p_nic <- if (nrow(nic_eq) > 0 && nrow(nic_bd) > 0) {
  make_nic_panel(nic_eq, "Equity") / make_nic_panel(nic_bd, "Bond")
} else if (nrow(nic_eq) > 0) {
  make_nic_panel(nic_eq, "Equity")
} else {
  make_nic_panel(nic_bd, "Bond")
}

# TikZ version: two separate ggplots via gridExtra::grid.arrange.
# Avoids patchwork tikzDevice rendering issues; each panel gets its own $\hat{h}_t$ y-label.
make_nic_panel_tikz <- function(df, title, show_xlab = TRUE) {
  snames <- unique(df$series)
  cols <- setNames(nic_pal[seq_along(snames)],  snames)
  ltys <- setNames(nic_ltys[seq_along(snames)], snames)
  p <- ggplot(df, aes(epsilon, h, colour = series, linetype = series)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = cols) +
    scale_linetype_manual(values = ltys) +
    labs(x = if (show_xlab) "$\\varepsilon_{t-1}$" else "",
         y = "$\\hat{h}_t$", title = title, colour = NULL, linetype = NULL) +
    theme_ces
  if (!show_xlab)
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                   axis.title.x = element_blank())
  p
}
p_nic_tikz <- if (HAS_TIKZ && nrow(nic_eq) > 0 && nrow(nic_bd) > 0) {
  pt1 <- make_nic_panel_tikz(nic_eq, "Equity", show_xlab = FALSE)
  pt2 <- make_nic_panel_tikz(nic_bd, "Bond",   show_xlab = TRUE)
  function() gridExtra::grid.arrange(pt1, pt2, ncol = 1)
} else NULL

save_pdf(p_nic, "fig1_volatility_nic.pdf", height = 6, p_tikz = p_nic_tikz)

# -------------------------------------------------------------------------
# FIGURES 2a / 2b & 3 — Correlation and Covariance NIC Surfaces
# -------------------------------------------------------------------------

if (!has_dcc) {
  cat("Skipping Figures 2a, 2b, 3 (dcc_fits.RData not found).\n")
} else {
  cat("--- Figures 2a/2b/3: NIC Surfaces ---\n")

  agdcc   <- fits[["AGDCC"]]
  psi     <- agdcc$psi_D
  P_mat   <- agdcc$P
  N_mat   <- agdcc$N

  # Choose pair: Germany and US equity
  i_nm <- "eq_usa"; j_nm <- "eq_germany"
  if (!all(c(i_nm, j_nm) %in% nms)) {
    i_nm <- eq_nms[1]; j_nm <- eq_nms[2]
  }
  ii <- which(nms == i_nm); jj <- which(nms == j_nm)

  a_i <- psi[paste0("[", i_nm, "].a1")]; g_i <- psi[paste0("[", i_nm, "].g1")]; b_i <- psi[paste0("[", i_nm, "].b1")]
  a_j <- psi[paste0("[", j_nm, "].a1")]; g_j <- psi[paste0("[", j_nm, "].g1")]; b_j <- psi[paste0("[", j_nm, "].b1")]

  P_ij  <- P_mat[ii, jj]; P_ii <- P_mat[ii, ii]; P_jj <- P_mat[jj, jj]
  N_ij  <- N_mat[ii, jj]; N_ii <- N_mat[ii, ii]; N_jj <- N_mat[jj, jj]

  # Q NIC functions: Qt given Q_{t-1} = P (long-run level).
  # Full update: Qt_ij = P_ij*(1 - a_i*a_j + b_i*b_j) - g_i*g_j*N_ij + shock terms.
  Q_ij_fn <- function(ei, ej) {
    P_ij * (1 - a_i * a_j + b_i * b_j) - g_i * g_j * N_ij +
      (a_i * a_j + g_i * g_j * (ei < 0) * (ej < 0)) * ei * ej
  }
  Q_ii_fn <- function(ei) {
    P_ii * (1 - a_i^2 + b_i^2) - g_i^2 * N_ii +
      (a_i^2 + g_i^2 * (ei < 0)) * ei^2
  }
  Q_jj_fn <- function(ej) {
    P_jj * (1 - a_j^2 + b_j^2) - g_j^2 * N_jj +
      (a_j^2 + g_j^2 * (ej < 0)) * ej^2
  }
  R_ij_fn <- function(ei, ej) {
    Q_ij_fn(ei, ej) / sqrt(pmax(Q_ii_fn(ei) * Q_jj_fn(ej), 1e-10))
  }

  eg_x <- seq(-3, 3, by = 0.15)
  eg_y <- seq(-3, 3, by = 0.15)
  eg <- seq(-3, 3, by = 0.15)
  z_corr   <- outer(eg_x, eg_y, Vectorize(R_ij_fn))
  R0       <- R_ij_fn(0, 0)
  z_impact <- z_corr - R0
  iclim    <- c(
    plyr::round_any(min(z_impact, na.rm = TRUE), 0.05, f = floor),
    plyr::round_any(max(z_impact, na.rm = TRUE), 0.05, f = ceiling)
  )
  

  # Figure 2a: 3D correlation NIC — pgfplots .tex + persp() PDF preview
  {
    skip  <- 2
    s_eg_x <- eg_x[seq(1, length(eg_x), by = skip)]
    s_eg_y <- eg_y[seq(1, length(eg_y), by = skip)]
    s_z    <- z_impact[seq(1, nrow(z_impact), by = skip),
                       seq(1, ncol(z_impact), by = skip)]

    # Export full grid as .dat for pgfplots (x varies fastest per scanline)
    n_x <- length(eg_x); n_y <- length(eg_y)
    dat_path <- file.path("output/charts", "fig2a_data.dat")
    write.table(
      data.frame(x = rep(eg_x, times = n_y),
                 y = rep(eg_y, each  = n_x),
                 z = as.vector(z_impact)),
      dat_path, row.names = FALSE, quote = FALSE, sep = " "
    )
    cat("Saved grid data:", dat_path, "\n")

    # pgfplots .tex — style matching nis_plot_style.tex
    tex_lines <- c(
      "% Auto-generated by 8_charts.R -- do not edit",
      "\\begin{tikzpicture}",
      "\\begin{axis}[",
      "  view={340}{5},",
      "  width=12cm, height=8cm,",
      paste0("  xlabel={\\footnotesize Shock: ", plabel(i_nm), "},"),
      paste0("  ylabel={\\footnotesize Shock: ", plabel(j_nm), "},"),
      "  zlabel={\\footnotesize Corr.\\ Impact},",
      "  tick label style={font=\\scriptsize},",
      "  label style={font=\\footnotesize},",
      "  grid=major,",
      "  grid style={densely dotted, gray!60},",
      paste0("  zmin=", iclim[1], ", zmax=", iclim[2], ","),
      "  colormap/viridis,",
      "  axis lines*=left,",
      "  axis line style={black},",
      "  tick align=outside,",
      "  xticklabel style={anchor=north, yshift=-2pt},",
      "  yticklabel style={anchor=north, yshift=-2pt},",
      "  xtick distance=1,",
      "  ytick distance=1,",
      "  ztick distance=0.05,",
      "  scaled z ticks=false,",
      "  z tick label style={",
      "    /pgf/number format/fixed,",
      "    /pgf/number format/precision=2",
      "  },",
      "]",
      paste0("\\addplot3[surf, shader=interp, mesh/cols=", n_x, "]"),
      "  table[x=x, y=y, z=z] {fig2a_data.dat};",
      paste0("\\addplot3[surf, mesh/cols=", n_x, ", draw=black, opacity=0.1]"),
      "  table[x=x, y=y, z=z] {fig2a_data.dat};",
      "\\end{axis}",
      "\\end{tikzpicture}"
    )
    writeLines(tex_lines, file.path("output/charts", "fig2a_corr_nic_3d.tex"))
    cat("Saved pgfplots TeX: output/charts/fig2a_corr_nic_3d.tex\n")

    # PDF preview (persp, subsampled)
    pdf(file.path("output/charts", "fig2a_corr_nic_3d.pdf"), width = 8, height = 5)
    op <- par(cex.axis = 0.7, cex.lab = 0.8)
    persp(s_eg_x, s_eg_y, s_z,
          xlab = paste("Shock:", plabel(i_nm)),
          ylab = paste("Shock:", plabel(j_nm)),
          zlab = "Impact", main = NULL,
          theta = 340, phi = 2, col = "lightblue", shade = .75,
          r = 100, d = 100, border = "black", ltheta = -135, lphi = 0,
          ticktype = "detailed", zlim = iclim)
    par(op)
    dev.off()
    cat("Saved PDF preview: output/charts/fig2a_corr_nic_3d.pdf\n")
  }
  
  # Figure 2b: 4-panel contour views
  save_pdf(function() {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    for (theta in c(210, 330, 30, 150)) {
      persp(eg, eg, z_impact,
            xlab = plabel(i_nm), ylab = plabel(j_nm),
            zlab = "Impact", main = sprintf("theta = %d deg", theta),
            theta = theta, phi = 20, col = "lightblue",
            shade = 0.4, ticktype = "simple", zlim = iclim)
    }
  }, "fig2b_corr_nic_views.pdf", width = 10, height = 8)

  # Figure 3: Covariance NIC surface (return units, matching CES Figure 3)
  ni_i <- tryCatch(newsimpact(stage1_best[[i_nm]]), error = function(e) NULL)
  ni_j <- tryCatch(newsimpact(stage1_best[[j_nm]]), error = function(e) NULL)
  if (!is.null(ni_i) && !is.null(ni_j)) {
    sigma_i <- sqrt(mean(h_mat[, ii], na.rm = TRUE))
    sigma_j <- sqrt(mean(h_mat[, jj], na.rm = TRUE))
    h_i_fn  <- approxfun(ni_i$zx, ni_i$zy, rule = 2)
    h_j_fn  <- approxfun(ni_j$zx, ni_j$zy, rule = 2)
    # Baseline: covariance at zero shock (eps=0), consistent with Fig 2a impact baseline
    cov_lr  <- sqrt(pmax(h_i_fn(0), 0)) * sqrt(pmax(h_j_fn(0), 0)) * R_ij_fn(0, 0)

    eg_ret <- seq(-0.05, 0.05, by = .1/100)
    z_cov3 <- outer(eg_ret, eg_ret, function(ri, rj) {
      # h_i_fn expects raw innovation units (same as ni$zx); R_ij_fn expects standardized residuals
      hi <- pmax(h_i_fn(ri), 0)
      hj <- pmax(h_j_fn(rj), 0)
      sqrt(hi) * sqrt(hj) * R_ij_fn(ri / sigma_i, rj / sigma_j) - cov_lr
    })
    clim3 <- quantile(z_cov3, c(0.01, 0.99), na.rm = TRUE)

    {
      n_ret <- length(eg_ret)
      dat3_path <- file.path("output/charts", "fig3_data.dat")
      write.table(
        data.frame(x = rep(eg_ret, times = n_ret),
                   y = rep(eg_ret, each  = n_ret),
                   z = as.vector(z_cov3)),
        dat3_path, row.names = FALSE, quote = FALSE, sep = " "
      )
      cat("Saved grid data:", dat3_path, "\n")

      # Compute a round z tick distance from clim3
      z3ticks    <- pretty(as.numeric(clim3), n = 5)
      z3tick_dist <- diff(z3ticks)[1]
      z3prec     <- max(4L, as.integer(ceiling(-log10(abs(z3tick_dist)))) + 1L)

      tex3_lines <- c(
        "% Auto-generated by 8_charts.R -- do not edit",
        "\\begin{tikzpicture}",
        "\\begin{axis}[",
        "  view={340}{5},",
        "  width=12cm, height=8cm,",
        paste0("  xlabel={\\footnotesize Return: ", plabel(i_nm), "},"),
        paste0("  ylabel={\\footnotesize Return: ", plabel(j_nm), "},"),
        "  zlabel={\\footnotesize Cov.\\ Change},",
        "  tick label style={font=\\scriptsize},",
        "  label style={font=\\footnotesize},",
        "  grid=major,",
        "  grid style={densely dotted, gray!60},",
        paste0("  zmin=", clim3[1], ", zmax=", clim3[2], ","),
        "  colormap/viridis,",
        "  axis lines*=left,",
        "  axis line style={black},",
        "  tick align=outside,",
        "  xticklabel style={anchor=north, yshift=-2pt},",
        "  yticklabel style={anchor=north, yshift=-2pt},",
        "  xtick distance=0.01,",
        "  ytick distance=0.01,",
        paste0("  ztick distance=", format(z3tick_dist, scientific = FALSE), ","),
        "  scaled x ticks=false,",
        "  scaled y ticks=false,",
        "  scaled z ticks=false,",
        "  x tick label style={/pgf/number format/fixed, /pgf/number format/precision=2},",
        "  y tick label style={/pgf/number format/fixed, /pgf/number format/precision=2},",
        "  z tick label style={",
        "    /pgf/number format/fixed,",
        paste0("    /pgf/number format/precision=", z3prec),
        "  },",
        "]",
        paste0("\\addplot3[surf, shader=interp, mesh/cols=", n_ret, "]"),
        "  table[x=x, y=y, z=z] {output/charts/fig3_data.dat};",
        paste0("\\addplot3[surf, mesh/cols=", n_ret, ", draw=black, opacity=0.1]"),
        "  table[x=x, y=y, z=z] {output/charts/fig3_data.dat};",
        "\\end{axis}",
        "\\end{tikzpicture}"
      )
      writeLines(tex3_lines, file.path("output/charts", "fig3_covariance_nic_3d.tex"))
      cat("Saved pgfplots TeX: output/charts/fig3_covariance_nic_3d.tex\n")

      # PDF preview
      pdf(file.path("output/charts", "fig3_covariance_nic_3d.pdf"), width = 8, height = 5)
      persp(eg_ret, eg_ret, z_cov3,
            xlab = paste("Return:", plabel(i_nm)),
            ylab = paste("Return:", plabel(j_nm)),
            zlab = "Covariance Change", main = NULL,
            theta = 340, phi = 2, col = "lightyellow", shade = 0.5,
            ticktype = "detailed", zlim = clim3)
      dev.off()
      cat("Saved PDF preview: output/charts/fig3_covariance_nic_3d.pdf\n")
    }
  }
}

# -------------------------------------------------------------------------
# FIGURE 4 — Conditional Equity Volatility (regional averages)
# -------------------------------------------------------------------------

cat("--- Figure 4: Conditional Equity Volatility ---\n")

vol_ann <- sqrt(h_mat * 52) * 100   # annualised %

make_vol_panel <- function(idx, region_label) {
  idx <- intersect(idx, seq_along(nms))
  if (length(idx) == 0) return(NULL)
  avg_vol <- rowMeans(vol_ann[, idx, drop = FALSE], na.rm = TRUE)
  df <- data.frame(date = dates, vol = avg_vol)
  ggplot(df, aes(date, vol)) +
    geom_line(colour = "steelblue", linewidth = 0.5) +
    geom_vline(data = events, aes(xintercept = date),
               linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    labs(x = NULL, y = "Ann. vol.", title = region_label) +
    scale_x_date(breaks = x_breaks, date_labels = "%Y") +
    theme_ces_panel
}

panels4 <- list(
  make_vol_panel(which(nms %in% reg$eq_emu),         "EMU Equity"),
  make_vol_panel(which(nms %in% reg$eq_europe_nemu), "Non-EMU Europe Equity"),
  make_vol_panel(which(nms %in% reg$eq_namerica),    "Americas Equity"),
  make_vol_panel(which(nms %in% reg$eq_australasia), "Australasia Equity")
)
panels4 <- Filter(Negate(is.null), panels4)

if (length(panels4) > 0)
  save_pdf(Reduce(`/`, panels4), "fig4_equity_volatility.pdf",
           width = 8, height = 2.5 * length(panels4))

# -------------------------------------------------------------------------
# FIGURE 5 — Conditional Bond Volatility (regional averages)
# -------------------------------------------------------------------------

cat("--- Figure 5: Conditional Bond Volatility ---\n")

vol_ann_bond <- sqrt(h_mat * 52) * 100

bond_panels_spec <- list(
  list(reg$bd_emu,      "EMU Bond Volatility"),
  list(reg$bd_namerica, "N. America Bond Volatility"),
  list(reg$bd_asia,     "Australasia Bond Volatility")
)

panels5 <- lapply(bond_panels_spec, function(x) {
  idx <- which(nms %in% x[[1]])
  if (length(idx) == 0) return(NULL)
  avg_vol <- rowMeans(vol_ann_bond[, idx, drop = FALSE], na.rm = TRUE)
  ggplot(data.frame(date = dates, vol = avg_vol), aes(date, vol)) +
    geom_line(colour = "tomato3", linewidth = 0.5) +
    geom_vline(data = events, aes(xintercept = date),
               linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    labs(x = NULL, y = "Ann. vol.", title = x[[2]]) +
    scale_x_date(breaks = x_breaks, date_labels = "%Y") +
    theme_ces_panel
})
panels5 <- Filter(Negate(is.null), panels5)

if (length(panels5) > 0)
  save_pdf(Reduce(`/`, panels5), "fig5_bond_volatility.pdf",
           width = 8, height = 2.5 * length(panels5))

# -------------------------------------------------------------------------
# FIGURE 6 — Conditional Equity Correlations (selected pairs)
# -------------------------------------------------------------------------

if (!has_dcc) {
  cat("Skipping Figures 6-10 (dcc_fits.RData not found).\n")
} else {
  cat("--- Figure 6: Conditional Equity Correlations ---\n")

  Rt  <- fits[["AGDCC"]]$Rt

  # Fig 6: upper-triangular 3×3 patchwork among France, Germany, Italy, UK
  make_panel6 <- function(p) {
    ii <- which(nms == p[1]); jj <- which(nms == p[2])
    df <- data.frame(date = dates, corr = Rt[ii, jj, ])
    ggplot(df, aes(date, corr)) +
      geom_line(colour = "steelblue", linewidth = 0.5) +
      geom_vline(data = events, aes(xintercept = date),
                 linetype = "dashed", colour = "grey50", linewidth = 0.4) +
      labs(x = NULL, y = "Cond. Corr.",
           title = paste(plabel(p[1]), "—", plabel(p[2]))) +
      scale_x_date(breaks = x_breaks, date_labels = "%Y") +
      theme_ces_panel
  }

  pairs6 <- list(
    c("eq_france",  "eq_germany"),
    c("eq_france",  "eq_italy"),
    c("eq_france",  "eq_united_kingdom"),
    c("eq_germany", "eq_italy"),
    c("eq_germany", "eq_united_kingdom"),
    c("eq_italy",   "eq_united_kingdom")
  )
  pairs6 <- Filter(function(p) all(p %in% nms), pairs6)

  if (length(pairs6) == 6) {
    ylim6 <- range(unlist(lapply(pairs6, function(p) {
      ii <- which(nms == p[1]); jj <- which(nms == p[2]); Rt[ii, jj, ]
    })), na.rm = TRUE)
    p6l <- lapply(pairs6, make_panel6)
    p6 <- patchwork::wrap_plots(
      p6l[[1]], p6l[[2]], p6l[[3]],
      patchwork::plot_spacer(), p6l[[4]], p6l[[5]],
      patchwork::plot_spacer(), patchwork::plot_spacer(), p6l[[6]],
      ncol = 3, widths = c(1, 1, 1), heights = c(1, 1, 1)
    ) & coord_cartesian(ylim = ylim6)
    save_pdf(p6, "fig6_equity_correlations.pdf", width = 12, height = 8)
  }

  # -----------------------------------------------------------------------
  # FIGURE 7 — Regional Average Equity Correlations
  # -----------------------------------------------------------------------

  cat("--- Figure 7: Regional Average Equity Correlations ---\n")

  avg_reg_corr <- function(nms_a, nms_b) {
    ia <- which(nms %in% nms_a); ib <- which(nms %in% nms_b)
    pairs <- if (identical(ia, ib)) {
      combn(ia, 2, simplify = FALSE)
    } else {
      as.list(as.data.frame(t(expand.grid(ia, ib))))
    }
    pairs <- Filter(function(p) p[1] != p[2], pairs)
    if (length(pairs) == 0) return(rep(NA_real_, T_obs))
    mat <- sapply(pairs, function(p) Rt[p[1], p[2], ])
    rowMeans(mat, na.rm = TRUE)
  }

  # All 6 upper-triangular pairs among 4 equity regions
  reg_pairs7 <- list(
    "EMU–non-EMU Europe"   = list(reg$eq_emu,         reg$eq_europe_nemu),
    "EMU–Americas"         = list(reg$eq_emu,         reg$eq_namerica),
    "EMU–Australasia"      = list(reg$eq_emu,         reg$eq_australasia),
    "Non-EMU–Americas"     = list(reg$eq_europe_nemu, reg$eq_namerica),
    "Non-EMU–Australasia"  = list(reg$eq_europe_nemu, reg$eq_australasia),
    "Americas–Australasia" = list(reg$eq_namerica,    reg$eq_australasia)
  )

  make_panel7 <- function(nms_a, nms_b, label) {
    corr_v <- avg_reg_corr(nms_a, nms_b)
    corr_sm <- zoo::rollmean(corr_v, k = 26, fill = NA, align = "right")
    df <- data.frame(date = dates, corr = corr_sm)
    df <- df[!is.na(df$corr), ]
    ggplot(df, aes(date, corr)) +
      geom_line(colour = "steelblue", linewidth = 0.6) +
      geom_vline(data = events, aes(xintercept = date),
                 linetype = "dashed", colour = "grey50", linewidth = 0.4) +
      labs(x = NULL, y = "Avg. Corr.", title = label) +
      scale_x_date(breaks = x_breaks, date_labels = "%Y") +
      theme_ces_panel
  }

  reg_args7 <- list(
    list(reg$eq_emu,         reg$eq_europe_nemu, "EMU — non-EMU Europe"),
    list(reg$eq_emu,         reg$eq_namerica,    "EMU — Americas"),
    list(reg$eq_emu,         reg$eq_australasia, "EMU — Australasia"),
    list(reg$eq_europe_nemu, reg$eq_namerica,    "Non-EMU — Americas"),
    list(reg$eq_europe_nemu, reg$eq_australasia, "Non-EMU — Australasia"),
    list(reg$eq_namerica,    reg$eq_australasia, "Americas — Australasia")
  )
  ylim7 <- range(unlist(lapply(reg_args7, function(x) avg_reg_corr(x[[1]], x[[2]]))),
                 na.rm = TRUE)
  p7l <- lapply(reg_args7, function(x) make_panel7(x[[1]], x[[2]], x[[3]]))
  p7 <- patchwork::wrap_plots(
    p7l[[1]], p7l[[2]], p7l[[3]],
    patchwork::plot_spacer(), p7l[[4]], p7l[[5]],
    patchwork::plot_spacer(), patchwork::plot_spacer(), p7l[[6]],
    ncol = 3, widths = c(1, 1, 1), heights = c(1, 1, 1)
  ) & coord_cartesian(ylim = ylim7)
  save_pdf(p7, "fig7_regional_equity_corr.pdf", width = 12, height = 8)

  # -----------------------------------------------------------------------
  # FIGURE 8 — Bond Return Correlations (regional averages)
  # -----------------------------------------------------------------------

  cat("--- Figure 8: Bond Return Correlations ---\n")

  # Upper-triangular 3-region bond correlation matrix: EMU, non-EMU Europe, Americas
  reg_pairs8 <- list(
    "EMU vs non-EMU Europe"      = list(reg$bd_emu,         reg$bd_europe_nemu),
    "EMU vs Americas"            = list(reg$bd_emu,         reg$bd_namerica),
    "non-EMU Europe vs Americas" = list(reg$bd_europe_nemu, reg$bd_namerica)
  )

  make_panel8 <- function(nms_a, nms_b, label) {
    corr_v <- avg_reg_corr(nms_a, nms_b)
    corr_sm <- zoo::rollmean(corr_v, k = 26, fill = NA, align = "right")
    df <- data.frame(date = dates, corr = corr_sm)
    df <- df[!is.na(df$corr), ]
    ggplot(df, aes(date, corr)) +
      geom_line(colour = "steelblue", linewidth = 0.6) +
      geom_vline(data = events, aes(xintercept = date),
                 linetype = "dashed", colour = "grey50", linewidth = 0.4) +
      labs(x = NULL, y = "Avg. Corr.", title = label) +
      scale_x_date(breaks = x_breaks, date_labels = "%Y") +
      theme_ces_panel
  }

  ylim8 <- range(unlist(lapply(list(
    list(reg$bd_emu, reg$bd_europe_nemu),
    list(reg$bd_emu, reg$bd_namerica),
    list(reg$bd_europe_nemu, reg$bd_namerica)
  ), function(x) avg_reg_corr(x[[1]], x[[2]]))), na.rm = TRUE)
  p8_tl <- make_panel8(reg$bd_emu,         reg$bd_europe_nemu, "EMU — non-EMU Europe")
  p8_tr <- make_panel8(reg$bd_emu,         reg$bd_namerica,    "EMU — Americas")
  p8_br <- make_panel8(reg$bd_europe_nemu, reg$bd_namerica,    "non-EMU Europe — Americas")
  p8 <- ((p8_tl | p8_tr) / (patchwork::plot_spacer() | p8_br)) &
        coord_cartesian(ylim = ylim8)
  save_pdf(p8, "fig8_bond_correlations.pdf", width = 10, height = 6)

  # -----------------------------------------------------------------------
  # FIGURE 9 — Selected Bond Pair Correlations
  # -----------------------------------------------------------------------

  cat("--- Figure 9: Selected Bond Correlations ---\n")

  # Upper-triangular: Germany vs Japan (top-left), Germany vs US (top-right),
  # Japan vs US (bottom-right); bottom-left is empty
  sel_bond_pairs <- list(
    c("bd_germany", "bd_japan"),
    c("bd_germany", "bd_united_states"),
    c("bd_japan",   "bd_united_states")
  )
  sel_bond_pairs <- Filter(function(p) all(p %in% nms), sel_bond_pairs)

  make_bond_panel9 <- function(p) {
    ii <- which(nms == p[1]); jj <- which(nms == p[2])
    df <- data.frame(date = dates, corr = Rt[ii, jj, ])
    ggplot(df, aes(date, corr)) +
      geom_line(colour = "steelblue", linewidth = 0.5) +
      geom_vline(data = events, aes(xintercept = date),
                 linetype = "dashed", colour = "grey50", linewidth = 0.4) +
      labs(x = NULL, y = "Conditional Correlation",
           title = paste(plabel(p[1]), "—", plabel(p[2]))) +
      scale_x_date(breaks = x_breaks, date_labels = "%Y") +
      theme_ces_panel
  }

  ylim9 <- range(unlist(lapply(sel_bond_pairs, function(p) {
    ii <- which(nms == p[1]); jj <- which(nms == p[2]); Rt[ii, jj, ]
  })), na.rm = TRUE)
  if (length(sel_bond_pairs) == 3) {
    p9_tl <- make_bond_panel9(sel_bond_pairs[[1]])
    p9_tr <- make_bond_panel9(sel_bond_pairs[[2]])
    p9_br <- make_bond_panel9(sel_bond_pairs[[3]])
    p9 <- ((p9_tl | p9_tr) / (patchwork::plot_spacer() | p9_br)) &
          coord_cartesian(ylim = ylim9)
  } else {
    panels9 <- lapply(sel_bond_pairs, make_bond_panel9)
    p9 <- Reduce(`|`, panels9) & coord_cartesian(ylim = ylim9)
  }
  save_pdf(p9, "fig9_bond_corr_selected.pdf", width = 10, height = 6)

  # -----------------------------------------------------------------------
  # FIGURE 10 — Equity-Bond Cross-Asset Correlations
  # -----------------------------------------------------------------------

  cat("--- Figure 10: Equity-Bond Cross-Asset Correlations ---\n")

  emu_bond_idx <- which(nms %in% reg$bd_emu)

  cross_corr <- function(bd_idx_v, eq_idx_v) {
    pairs <- as.list(as.data.frame(t(expand.grid(bd_idx_v, eq_idx_v))))
    if (length(pairs) == 0) return(rep(NA_real_, T_obs))
    mat <- sapply(pairs, function(p) Rt[p[1], p[2], ])
    rowMeans(mat, na.rm = TRUE)
  }

  reg_pairs10 <- list(
    "EMU Bonds vs EMU Equity"         = list(emu_bond_idx, which(nms %in% reg$eq_emu)),
    "EMU Bonds vs Americas Equity"    = list(emu_bond_idx, which(nms %in% reg$eq_namerica)),
    "EMU Bonds vs Australasia Equity" = list(emu_bond_idx, which(nms %in% reg$eq_australasia))
  )

  df10 <- purrr::map_dfr(names(reg_pairs10), function(lb) {
    p <- reg_pairs10[[lb]]
    data.frame(date = dates, corr = cross_corr(p[[1]], p[[2]]), region = lb)
  })
  df10 <- df10[!is.na(df10$corr), ]
  df10 <- df10 %>%
    group_by(region) %>% arrange(date) %>%
    mutate(corr_sm = zoo::rollmean(corr, 26, fill = NA, align = "right")) %>%
    ungroup()

  df10$region <- factor(df10$region, levels = names(reg_pairs10))

  ylim10 <- range(df10$corr_sm, na.rm = TRUE)
  panels10 <- lapply(levels(df10$region), function(lb) {
    sub10 <- df10[df10$region == lb, ]
    ggplot(sub10, aes(date, corr_sm)) +
      geom_line(colour = "steelblue", linewidth = 0.6) +
      geom_vline(data = events, aes(xintercept = date),
                 linetype = "dashed", colour = "grey50", linewidth = 0.4) +
      labs(x = NULL, y = "Avg. Corr. (26-wk smooth)", title = lb) +
      scale_x_date(breaks = x_breaks, date_labels = "%Y") +
      theme_ces_panel
  })

  p10 <- Reduce(`/`, panels10) & coord_cartesian(ylim = ylim10)
  save_pdf(p10, "fig10_equity_bond_corr.pdf", width = 8, height = 8)
}

cat("\nDone. Charts written to output/charts/\n")
