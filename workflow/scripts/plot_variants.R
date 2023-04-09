# plot_variants.R
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
output_file <- args[2]

# Read VCF file and extract quality scores
extract_quality_scores <- function(vcf_file) {
  qual_scores <- c()
  con <- file(vcf_file, "r")
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
      break
    }
    if (substring(line, 1, 1) != "#") {
      fields <- strsplit(line, "\t")[[1]]
      qual_scores <- c(qual_scores, as.numeric(fields[6]))
    }
  }
  close(con)
  return(qual_scores)
}

qual_scores <- extract_quality_scores(vcf_file)

# Create data frame for ggplot
qual_scores_df <- data.frame(QualityScore = qual_scores)

# Generate scatterplot using ggplot2
scatterplot <- ggplot(qual_scores_df, aes(x = seq_along(QualityScore), y = QualityScore)) +
  geom_point() +
  labs(title = "Scatterplot of Variant Calling Quality Scores",
       x = "Variant Index",
       y = "Quality Score") +
  theme_minimal()

# Save scatterplot to file
ggsave(filename = output_file, plot = scatterplot, width = 10, height = 6, dpi = 300)
