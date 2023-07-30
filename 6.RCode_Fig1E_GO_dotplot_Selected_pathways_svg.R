setwd("/Users/sohyunmoon/Desktop/Collaboration_Work/NJ/Project4_Abdel_RNAseq_NG_S961/Code_Manuscript")
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Cairo))

CairoFonts(
  regular="Arial:style=Regular",
  bold="Arial:style=Bold",
  italic="Arial:style=Italic",
  bolditalic="Arial:style=Bold Italic,BoldItalic",
  symbol="Symbol"
)

# ===========================================

log10_floor <- function(x) {
  10^(floor(log10(x)))
}

stopifnot(log10_floor(.004) == .001)
stopifnot(log10_floor(.04) == .01)
stopifnot(log10_floor(.4) == .1)
stopifnot(log10_floor(4) == 1)
stopifnot(log10_floor(44) == 10)
stopifnot(log10_floor(444) == 100)
stopifnot(log10_floor(4444) == 1000)

scooch_out <- function(x) {
  if (x == 0) return(0)
  lbound <- log10_floor(x)
  floor((x + lbound) / lbound) * lbound
}
stopifnot(scooch_out(.0044) == 0.005)
stopifnot(scooch_out(.004) == 0.005)
stopifnot(scooch_out(.044) == 0.05)
stopifnot(scooch_out(.04) == 0.05)
stopifnot(scooch_out(.44) == 0.5)
stopifnot(scooch_out(.4) == 0.5)
stopifnot(scooch_out(0) == 0)
stopifnot(scooch_out(4) == 5)
stopifnot(scooch_out(40) == 50)
stopifnot(scooch_out(44) == 50)
stopifnot(scooch_out(400) == 500)
stopifnot(scooch_out(444) == 500)
stopifnot(scooch_out(4000) == 5000)
stopifnot(scooch_out(4444) == 5000)

make_graph <- function(csv_filename, svg_filename, plot_title, plot_cb = NULL) {
  print(csv_filename)

  GO <- read.csv(csv_filename, sep = ",", header = TRUE, row.names = NULL)

   if (any(grepl("REAC", csv_filename))) {
    GO$term_name <- lapply(strwrap(GO$term_name, width = 35, simplify = FALSE), paste, collapse="\n")
  }
    else if (any(grepl("process", GO$term_name))) {
    GO$term_name <- gsub(".process.*","", GO$term_name)
  }

  GO_top20 <- head(GO[order(GO$adjusted_p_value),], n=20)
  GO_top20$term_name <- factor(
    GO_top20$term_name,
    levels=GO_top20$term_name[rev(order(GO_top20$adjusted_p_value))]
  )

  x_scooch <- scooch_out(max(GO_top20$adjusted_p_value))

  size_breaks <- scales::breaks_extended()(GO_top20$intersection_size)
  
  size_breaks[1] <- floor(min(GO_top20$intersection_size))
  size_breaks[length(size_breaks)] <- ceiling(max(GO_top20$intersection_size))

  CairoSVG(svg_filename, width = 5, height = 8, unit = "px")

    ggplot(GO_top20, aes(x = adjusted_p_value, y = term_name,
                         size = intersection_size, color = Gene_ratio)) + 
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() + 
    scale_size(
      limits = c(size_breaks[1], size_breaks[length(size_breaks)]),
      breaks = size_breaks,
    ) +
    scale_x_continuous(
      limits = c(-0, x_scooch)
    ) +
    ylab(NULL) + # hide y-axis label
    xlab("adjusted_p_value") +
    labs(
      title = plot_title,
      x = "adjusted_p_value",
      color="Gene_ratio"
    ) +
    theme(
      axis.text.y = element_text(size = 0, color = "black"),
      legend.title = element_text(size = 7.5, color = "black")
    )+
    guides(
      size = guide_legend(order = 1)
      )
}


make_graph("NG_Pathway_Selelcted_FDR005_Name_GO_KEGG_REAC.csv", "NG_Pathway_Selected.svg", "GO_NG_enrichment analysis") 

