#!/usr/bin/env Rscript
library(optparse)
library(ape)

option.list <- list(
  make_option(c("-c", "--colour"),
              help = "Colour the branches using a provided colour file,
                 where first column is the label and
                 second column is the colour."),
  make_option(c("-l", "--legend"),
                help = "Generate a legend of line colour used, where
                 first column is the label and
                 second column is the colour."),
  make_option(c("-L", "--label"), action = "store_true",
              help = "Add label to tree tip (Default: Do not add label).",
              default = FALSE),
  make_option(c("-e", "--edge_colourisation"),
              help = "Algorithm for colouring edges (Default: complete).
                 Choices are: \"complete\", \"tip_only\"",
              default = "complete"),
  make_option(c("-O", "--outtree"),
              help = "Output file tree name (Default: no output)."),
  make_option(c("-o", "--output"),
              help = "Output file name (Default: input file name + pdf)."),
  make_option(c("-t", "--type"),
              help = "Phylogenetic tree type (Default: unrooted).
                 Choices are: \"phylogram\", \"cladogram\",
                 \"fan\", \"unrooted\", \"radial\"",
              default = "unrooted"),
  make_option(c("-m", "--method"),
              help = "Clustering method (Default: Neighbor joining).
                 Choices are: \"nj\", \"upgma\"",
              default = "nj"),
  make_option("--cex",
              type = "double",
              help = "Font size for leaf label (Default: 0.7)",
              default = 0.7)
)
parser <- OptionParser(usage = "%prog [options] distance",
                       option_list = option.list)
args <- parse_args(parser, positional_arguments = 1)

distance <- as.dist(read.delim(args$args, check.names=F))
if (args$options$method == "upgma") {
  phylo <- as.phylo(hclust(distance, method = "average"))
} else {
  phylo <- nj(distance)
}
tree <- reorder.phylo(phylo, order = "postorder")

# second column of the edge matrix is the tip.label index
if (is.null(args$options$colour)) {
  tip.colours <- "#191919"
  edge.colours <- "#191919"
} else {
  annotations <- read.delim(args$options$colour)
  tip.colours <- as.vector(annotations$COLOUR)
  # sanity check
  if (Ntip(tree) != length(tip.colours)) {
    warning("Number of sample (", Ntip(tree), ") does not match with ",
            "number of colours (", length(tip.colours), ")")
  }
  edge.colours <- tip.colours[tree$edge[, 2]]
  if (args$options$edge_colourisation == "tip_only") {
    edge.colours[is.na(edge.colours)] <- "#E5E4E2"
  } else if (args$options$edge_colourisation == "complete") {
      edge.colours[is.na(edge.colours)] <- "####"
      # FIXME: how to traverse previous and adjacent edge in R idiom
      for (edge in 1:nrow(tree$edge)) {
        if (edge.colours[edge] == "####") {
            nodes <- which(tree$edge[, 1] == tree$edge[edge, 2])
            if (edge.colours[nodes[1]] == edge.colours[nodes[2]]) {
              edge.colours[edge] <- edge.colours[nodes[1]]
            } else {
              edge.colours[edge] <- "#E5E4E2"
            }
        }
      }
  } else {
      stop("Unknown parameters supplied: ", args$options$edge_colourisation)
  }
}

if (is.null(args$options$output)) {
  output.name <- basename(tools::file_path_sans_ext(args$args))
  if (!is.null(args$options$colour)) {
    output.name <- paste(output.name, 
        basename(tools::file_path_sans_ext(args$options$colour)), sep = "_")
    if (!is.null(args$options$legend)) {
      output.name <- paste(output.name,
          basename(tools::file_path_sans_ext(args$options$legend)), sep = "_")
    }
    output.name <- paste(output.name, args$options$edge_colourisation,
                         sep = "_")
  }
  if (args$options$label) {
    output.name <- paste(output.name, "use", "label", sep = "_")
  } else {
      output.name <- paste(output.name, "no", "label", sep = "_")
  }
  output.name <- paste(output.name, args$options$type,
                       args$options$method, sep = "_")
  output.name <- paste(output.name, "pdf", sep = ".")
} else {
  output.name <- args$options$output
}
#output.file <- paste(output.name, "pdf", sep = ".")
output.file <- output.name
pdf(file = output.file, title = output.name, width = 16.5, height = 9.27)
if (!is.null(args$options$outtree)) {
  outtree.file <- args$options$outtree
  write.tree(tree, file = outtree.file)
}

phylo.type <- args$option$type
use.label <- args$options$label
plot.phylo(tree, type = phylo.type, show.tip.label = use.label,
           tip.color = tip.colours, lab4ut = "axial",
           edge.color = edge.colours, edge.width = 0.5,
           cex = args$options$cex, rotate.tree = 90, no.margin = TRUE)

if (!is.null(args$options$colour) && !is.null(args$options$legend)) {
  legends <- read.delim(args$options$legend)
  legends.group <- as.vector(legends$GROUP)
  legends.colour <- as.vector(legends$COLOUR)
  legend("topleft", legend = legends.group, col = legends.colour,
         lwd = 3)
}
