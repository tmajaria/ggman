# hg19 chromosome sizes
.hg19 <<- c(0,           # chr1
          249250621,
          492449994,
          690472424,
          881626700,
          1062541960,
          1233657027,
          1392795690,
          1539159712,
          1680373143,  # chr10
          1815907890,
          1950914406,
          2084766301,
          2199936179,
          2307285719,
          2409817111,
          2500171864,
          2581367074,
          2659444322,
          2722469842,  # chr20
          2781598825,
          2829728720,  # chr22
          2881033286,  # chrX (w/o PAR)
          3033334811,  # chrY (w/o PAR)
          3089739342,  # chrXY (PAR1 & PAR2)
          3095677412)

.hg19.breaks <<- c(124625310,
                 370850307,
                 591461209,
                 786049562,
                 972084330,
                 1148099493,
                 1313226358,
                 1465977701,
                 1609766427,
                 1748140516,
                 1883411148,
                 2017840353,
                 2142351240,
                 2253610949,
                 2358551415,
                 2454994487,
                 2540769469,
                 2620405698,
                 2690957082,
                 2752034333,
                 2805663772,
                 2855381003,
                 2957184048,
                 3061537076,
                 3092708377,
                 3095685697)

.hg38 <<- c(0,          #chr1
          248956422,
          491149951,
          689445510,
          879660065,
          1061198324,
          1232004303,
          1391350276,
          1536488912,
          1674883629,   #chr10
          1808681051,
          1943767673,
          2077042982,
          2191407310,
          2298451028,
          2400442217,
          2490780562,
          2574038003,
          2654411288,
          2713028904,    #chr20
          2777473071,
          2824183054,     # chr22
          2875001522,     
          3031042417,
          3088269832,
          3095677412)

.hg38.breaks <<- c(124478211,
                  370053186,
                  590297730,
                  784552788,
                  970429194,
                  1146601314,
                  1311677290,
                  1463919594,
                  1605686270,
                  1741782340,
                  1876224362,
                  2010405328,
                  2134225146,
                  2244929169,
                  2349446622,
                  2445611390,
                  2532409282,
                  2614224646,
                  2683720096,
                  2745250988,
                  2800828062,
                  2849592288,
                  2953021970,
                  3059656124,
                  3092708377,
                  3095685697)

.labels <<- c(as.character(seq(22)), "X", "Y", "XY", "MT")

scale_colour_dichromatic = function(values = c("grey30", "grey60")) {
  requireNamespace('ggplot2')
  values = rep(values, 13)
  scale_color_manual(values = values)
}

scale_colour_traditional = function() {
  requireNamespace('ggplot2')
  traditional = rep(c("#000000", "#FF0000", "#008B00", "#0000FF", "#454545", "#EE00EE", "#009ACD", "#EE7600"), 4)
  scale_color_manual(values = traditional)
}

theme_publication = function(base_size = 16, base_family = "Helvetica", ...) {
  requireNamespace('ggplot2')
  thm = theme_bw(base_size = base_size, base_family = base_family, ...)
  thm + theme(axis.line.x = element_blank(),
              panel.background  = element_blank(),
              panel.grid.major  = element_blank(),
              panel.grid.minor  = element_blank())
}

..convert2posXR <- function(chr, bp, chr_pos){
  n_var = length(bp)
  pos <- c()
  for (i in 1:n_var){
    pos <- c(pos, bp[i] + chr_pos[chr[i]])
  }
  return(pos)
}

.convert2posXR <- function(chr, bp, build) {
  if (length(chr) != length(bp)) {
    stop("SIZE DIFFER.");
  }

  n_chr = length(unique(chr))
  if (n_chr == 1) {
    return (list(posX = bp,
                 breaks = ggplot2::waiver(),
                 labels = ggplot2::waiver(),
                 xlabel = paste("Chromosome", chr[1])))
  }

  max_chr = max(chr)
  if (build == "hg19") {
    return (list(posX = ..convert2posXR(chr, bp, .hg19),
                 breaks = .hg19.breaks[1:max_chr],
                 labels = .labels[1:max_chr],
                 xlabel = "Chromosome"))
  } else {
    return (list(posX = ..convert2posXR(chr, bp, .hg38),
                 breaks = .hg38.breaks[1:max_chr],
                 labels = .labels[1:max_chr],
                 xlabel = "Chromosome"))
  }
}



ggmanhattanGrouped <- function(data, SNP = "SNP", chr = "CHR", bp = "BP", P = "P", group = "group_id", color = "color", P_char = NULL, logP = TRUE, build = 'hg19',
                        nominal = c(1.0e-5),
                        significance = c(5.0e-8), ylim = NULL,
                        lead_snp = NULL, annotate_snp = NULL,
                        theme_base = theme_publication(),
                        scale_color = scale_colour_dichromatic(),
                        highlight = NULL, highlight_col = c("mediumblue", "deeppink"),
                        plot.grid = FALSE,
                        expand.x = c(0.03, 0.03), expand.y = c(0.03, 0.03)) {
  requireNamespace('ggplot2')

  # initial checks on data
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

  # get absolute positions for variants
  conv = .convert2posXR(data$CHR, data$BP, build)
  data$x = conv$posX

  # convert to log scale
  data$y = if (logP) -log10(data$P) else data$P
  
  # define colors for each chromosome

  # grey scale for non-highlighted groups
  grey_vals = rep(c("gray30", "gray50"), 20)

  # set colors based on chromosomes
  data$COLOR <- as.factor(data$CHR)
  data$COLOR <- grey_vals[data$CHR]

  # add colors to groups if we have any
  data$GROUP <- factor( data$GROUP, levels =  unique(data$GROUP[order(data$GROUP)]))
  data[!is.na(data[[color]]), "COLOR"] <- as.character(data[!is.na(data[[color]]), color])
  data$COLOR <- as.factor(data$COLOR)

  # make list of all colors for later in ggplot
  all.cols <- levels(data$COLOR)
  names(all.cols) <- all.cols

  # add alpha level (this may be removed later)
  data$alpha <- "1"
  data[data$COLOR %in% grey_vals, "alpha"] <- "0.3"
  data$alpha <- as.factor(data$alpha)

  # make list of all colors for later in ggplot
  all.alpha <- levels(data$alpha)
  names(all.alpha) <- all.alpha

  ### Plots

  if (!all(levels(data$COLOR) %in% grey_vals)){
    print("Multicolor plot")
    plt = ggplot(data) + geom_point(data = base::subset(data, COLOR %in% grey_vals), aes(x, y, color = COLOR), size = .5) +
          geom_point(data = base::subset(data, !(COLOR %in% grey_vals)), aes(x, y, color = COLOR), size = 1) +
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
    print("Monochromatic plot")
    plt = ggplot(data) + geom_point(data = base::subset(data, COLOR %in% grey_vals), aes(x, y, color = COLOR), size = .5) +
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
  
  return(plt)
}