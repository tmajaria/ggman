

#######

#' Plot a manhattan plot using ggplot2.
#'
#' @param data
#' @param SNP
#' @param chr
#' @param bp
#' @param P
#' @param logP
#' @param build
#' @param theme_base
#' @param scale_color
#' @export
ggmanhattanTM <- function(data, SNP = "SNP", chr = "CHR", bp = "BP", P = "P", group = "group_id", P_char = NULL, logP = TRUE, build = 'hg19',
                        nominal = c(1.0e-5),
                        significance = c(5.0e-8), ylim = NULL,
                        lead_snp = NULL, annotate_snp = NULL,
                        theme_base = theme_publication(),
                        scale_color = scale_colour_dichromatic(),
                        highlight = NULL, highlight_col = c("mediumblue", "deeppink"),
                        plot.grid = FALSE,
                        expand.x = c(0.03, 0.03), expand.y = c(0.03, 0.03)) {
  requireNamespace('ggplot2')

  idx = match(c(SNP, chr, bp, P, group), colnames(data))
  if (sum(is.na(idx)) > 0) {
    stop(paste(c(SNP, chr, bp, P, group)[which(is.na(idx))], "not found.", sep = " "))
  }
  colnames(data)[idx] = c("SNP", "CHR", "BP", "P", "GROUP")

  if (is.null(data$BP)) {
    stop("NULL BP")
  }
  if (is.null(data$P)) {
    stop("NULL P")
  }
  if (is.null(data$SNP) & (!is.null(lead_snp) | !is.null(annotate_snp))) {
    stop("NULL SNP")
  }
  if (is.null(data$GROUP)) {
    stop("NULL GROUP")
  }
  if (is.function(theme_base)) {
    theme_base = theme_base()
  }
  if (is.function(scale_color)) {
    scale_color = scale_color()
  }

  conv = .convert2posX(data$CHR, data$BP, build)
  data$x = conv$posX
  data$y = if (logP) -log10(data$P) else data$P
  
  grey_vals = rep(c("grey70", "grey80"), 13)
  # col_vals = rep(c("#FF0000", "#008B00", "#0000FF", "#ffff00", "#EE00EE", "#009ACD", "#EE7600"), 5)
  col_vals <- rep("blue", 30)
  # all.cols <- c("#FF0000" = "#FF0000", "#008B00" = "#008B00", "#0000FF" = "#0000FF", "#ffff00" = "#ffff00", "#EE00EE" = "#EE00EE", "#009ACD" = "#009ACD", "#EE7600" = "#EE7600","grey40" = "grey40", "grey70" = "grey70")
  all.cols <- c("blue" = "blue", "grey70" = "grey70", "grey80" = "grey80")
  all.alpha <- c("0.3" = 0.3, "1" = 1)
  data$color <- as.factor(data$CHR)
  data$color <- grey_vals[data$CHR]

  # get groups with P < nominal
  # nom_groups <- unique(data[data$P < nominal, ]$GROUP)
  # nom_groups <- nom_groups[nom_groups != "NA"]

  # if (length(nom_groups) > 0){
  #   data[data$GROUP %in% nom_groups, "color"] <- col_vals[data[data$GROUP %in% nom_groups, "CHR"]]
  # }

  data$color <- grey_vals[data$CHR]
  data[data$P < nominal & data$GROUP != "NA", "color"] <- col_vals[data[data$P < nominal & data$GROUP != "NA", "CHR"]]

  data[data$GROUP %in% data[data$color %in% col_vals,"GROUP"], "color"] <- col_vals[1]

  data$color <- as.factor(data$color)
  data$alpha <- "1"
  data[data$color %in% grey_vals, "alpha"] <- "0.3"

  data$GROUP <- factor( data$GROUP, levels =  unique(data$GROUP[order(data$GROUP)]))

  if (nrow(data[data$color %in% col_vals,]) > 0){
    plt = ggplot(data) + geom_point(data = base::subset(data, color %in% c("grey70","grey80")), aes(x, y, color = color), size = .5) +
          geom_point(data = base::subset(data, color == "blue"), aes(x, y, color = color), size = .5) +
          geom_hline(yintercept = -log10(nominal), linetype = "dashed", color = "red") +
          geom_hline(yintercept = -log10(significance), linetype = "dashed", color = "black") +
          scale_x_continuous(breaks = conv$breaks, labels = conv$labels, expand = expand.x) +
          scale_y_continuous(expand = expand.y) +
          theme_base +
          scale_color_manual(values = all.cols) +
          scale_alpha_manual(values = all.alpha) + 
          theme(axis.text.x = element_text(size = rel(0.5)), legend.position = "none") +
          theme(axis.text.y = element_text(size = rel(0.5)), legend.position = "none") +
          theme(axis.title.x = element_text(size = rel(0.5)), legend.position = "none") +
          theme(axis.title.y = element_text(size = rel(0.5)), legend.position = "none") +
          xlab(conv$xlabel) + ylab(expression(-log[10](italic(P))))
  } else {
    plt = ggplot(data) + geom_point(data = base::subset(data, color %in% c("grey70","grey80")), aes(x, y, color = color), size = .5) +
          geom_hline(yintercept = -log10(nominal), linetype = "dashed", color = "red") +
          geom_hline(yintercept = -log10(significance), linetype = "dashed", color = "black") +
          scale_x_continuous(breaks = conv$breaks, labels = conv$labels, expand = expand.x) +
          scale_y_continuous(expand = expand.y) +
          theme_base +
          scale_color_manual(values = all.cols) +
          scale_alpha_manual(values = all.alpha) + 
          theme(axis.text.x = element_text(size = rel(0.5)), legend.position = "none") +
          theme(axis.text.y = element_text(size = rel(0.5)), legend.position = "none") +
          theme(axis.title.x = element_text(size = rel(0.5)), legend.position = "none") +
          theme(axis.title.y = element_text(size = rel(0.5)), legend.position = "none") +
          xlab(conv$xlabel) + ylab(expression(-log[10](italic(P))))
  }
  

  # if (length(nom_groups) > 0){
  # plt = ggplot(data) + geom_point(data = base::subset(data, color %in% c("grey70", "grey80")), aes(x, y, color = color)) +
  #         geom_point(data = base::subset(data, color == "blue"), aes(x, y, color = color)) +
  #         geom_hline(yintercept = -log10(nominal), linetype = "dashed", color = "red") +
  #         geom_hline(yintercept = -log10(significance), linetype = "dashed", color = "black") +
  #         scale_x_continuous(breaks = conv$breaks, labels = conv$labels, expand = expand.x) +
  #         scale_y_continuous(expand = expand.y) +
  #         theme_base +
  #         scale_color_manual(values = all.cols) +
  #         scale_alpha_manual(values = all.alpha) + 
  #         theme(axis.text.x = element_text(size = rel(0.75)), legend.position = "none") +
  #         xlab(conv$xlabel) + ylab(expression(-log[10](italic(P))))
  # } else {
  #   plt = ggplot(data) + geom_point(data = data, aes(x, y, color = color)) +
  #         geom_hline(yintercept = -log10(nominal), linetype = "dashed", color = "red") +
  #         geom_hline(yintercept = -log10(significance), linetype = "dashed", color = "black") +
  #         scale_x_continuous(breaks = conv$breaks, labels = conv$labels, expand = expand.x) +
  #         scale_y_continuous(expand = expand.y) +
  #         theme_base +
  #         scale_color_manual(values = all.cols) +
  #         scale_alpha_manual(values = all.alpha) + 
  #         theme(axis.text.x = element_text(size = rel(0.75)), legend.position = "none") +
  #         xlab(conv$xlabel) + ylab(expression(-log[10](italic(P))))
  # }

  if (!is.null(lead_snp)) {
    lead_snp = subset(data, data$SNP %in% lead_snp)
  }

  if (!is.null(annotate_snp)) {
    lsnp = if(!is.null(ylim)) subset(lead_snp, y < ylim[2]) else lead_snp
    plt = plt + geom_text(data = lsnp, aes(label = data$SNP), color = "black", angle = 90, hjust = -0.1)
  }

  if (!is.null(highlight)) {
    if (!is.list(highlight)) highlight = list(highlight)
    for (i in 1:length(highlight)) {
      plt = plt + geom_point(data = subset(data, data$SNP %in% highlight[[i]]), color = highlight_col[i])
      if (!is.null(lead_snp)) {
        lead_snp$color[lead_snp$SNP %in% highlight[[i]]] = i
      }
    }
    lead_snp$color = factor(lead_snp$color, levels = 1:length(highlight))
  }


  if (!is.null(ylim)) {
    plt = plt + coord_cartesian(ylim = ylim)
    print(lead_snp)
    print(ylim)
    if (!is.null(lead_snp)) {
      lead_snp = subset(lead_snp, y >= ylim[2])
      scale_color_highlight = if(!is.null(highlight)) scale_color_manual(values = highlight_col) else scale_color
      if (is.null(P_char)) {
        lead_snp$label = if(logP) stringr::str_c(sprintf("\u25BA"), sprintf("%.1e",lead_snp$P), sep = "") else sprintf("\u25BA")
      } else if (P_char != FALSE) {
        lead_snp$label = stringr::str_c(sprintf("\u25BA"), lead_snp[[P_char]], sep = "")
      } else {
        lead_snp$label = sprintf("\u25BA")
      }
      plt_annot = ggplot(lead_snp) +
                    geom_text(aes(x, 0, color = color, label = label), size = 3, vjust = 0.5, angle = 90) +
                    scale_x_continuous(limits = c(0, max(data$x)), expand = expand.x) +
                    scale_color_highlight +
                    theme_void() + theme(plot.margin = margin(6,6,-2,6), legend.position = "none")
      g1 = ggplotGrob(plt_annot)
      g2 = ggplotGrob(plt + theme(plot.margin = margin(-2,6,6,6)))
      maxWidth = grid::unit.pmax(g1$widths[2:5], g2$widths[2:5])
      g1$widths[2:5] <- as.list(maxWidth)
      g2$widths[2:5] <- as.list(maxWidth)
      g = gridExtra::arrangeGrob(g1, g2, heights = c(1, 11))
      if (plot.grid) {
        grid::grid.newpage()
        grid::grid.draw(g)
      }
      return(g)
    }
  } else if (!is.null(annotate_snp)) {
    plt = plt + coord_cartesian(ylim = c(0, ceiling(max(data$y) + 1)))
  }
  return(plt)
}