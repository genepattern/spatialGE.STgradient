# Load necessary libraries
if (!require("optparse"))  stop("The optparse package is required but is not installed.")
if (!require("spatialGE")) stop("The spatialGE package is required but is not installed.")

library('optparse')
library('spatialGE')
library('tibble')
library('dplyr')

# Parse the command line options
parser = OptionParser()
parser <- add_option(parser, c("-f", "--file"), type = 'character', help = "input.file")
parser <- add_option(parser, c("-v", "--ngenes"),  type = 'integer', help = "num.variable.genes")
parser <- add_option(parser, c("-r", "--robust"),  type = "character", help = "robust.regression")
parser <- add_option(parser, c("-i", "--ignore"),  type = "character", help = "ignore.outliers")
parser <- add_option(parser, c("-l", "--limit"),  type = "double", help = "correlation.limit")
parser <- add_option(parser, c("-n", "--neighbors"),  type = 'integer', help = "min.neighbors")
parser <- add_option(parser, c("-d", "--distance"),  type = "character", help = "distance.summary")
parser <- add_option(parser, c("-s", "--samples"),  type = "character", help = "samples")
parser <- add_option(parser, c("-a", "--annotation"),  type = "character", help = "annotation.to.test")
parser <- add_option(parser, c("-c", "--cluster"),  type = "character", help = "reference.cluster")
parser <- add_option(parser, c("-e", "--exclude"),  type = "character", help = "exclude.clusters")
args <- parse_args(parser)

# Load the RDS file
data <- readRDS(args$file)

# Process arguments
robust_regression <- tolower(args$robust) == "true"
ignore_outliers <- tolower(args$ignore) == "true"
correlation_limit <- if (args$limit != "") args$limit else NULL
samples <- if (args$samples != "") unlist(strsplit(args$samples, ",")) else NULL
exclude_clusters <- if (args$exclude != "") unlist(strsplit(args$exclude, ",")) else NULL
annotations = if (args$annotation != "") list(args$annotation) else NULL
reference_cluster = if (args$cluster != "") list(args$cluster) else NULL

# If all annotations, assemble the list
if (is.null(annotations) {
    annotations <- c()
    for (t in data@spatial_meta[]) annotations <- unique(c(annotations, unlist(colnames(t)))
    annotations <- annotations[!annotations %in% c('libname', 'xpos', 'ypos', 'total_counts', 'total_genes')]
    # colnames(x@spatial_meta[[i]])
}

# Call STgradient
for (a in annotations) {
    print(paste('Beginning loop for annotation: ', a))
    if (is.null(reference_cluster)) reference_cluster <- 1:as.integer(sub(".*_k(\\d+)$", "\\1", a))

    for (c in reference_cluster) {
        print(paste('Beginning loop for reference cluster: ', c))
        mod_stlist <- STgradient(
            x = data,
            samples = samples,
            topgenes = args$ngenes,
            annot = a,
            ref = c,
            exclude = exclude_clusters,
            out_rm = ignore_outliers,
            limit = correlation_limit,
            distsumm = args$distance,
            min_nb = args$neighbors,
            robust = robust_regression)
            # cores=1)  # Multiple cores may not work on macOS

        # Lazily create subdirectory for result files
        if (!dir.exists(a)) dir.create(a, recursive = TRUE)

        lapply(names(mod_stlist), function(i) {
            file_name <- paste0('stgradient_', i, '_ref_', c, '.csv')
            file_path <- file.path(a, file_name)
            print(file_path)
            write.csv(mod_stlist[[i]], file=file_path, row.names=F, quote=F)
        })
    }
}
