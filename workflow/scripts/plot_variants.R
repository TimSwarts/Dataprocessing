# Load the required R packages
library(ggplot2)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
output_file <- args[2]

# Read the VCF file
data <- read.table(vcf_file, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)

# Set columnames to standard VCF format
colnames(data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

# Create a scatter plot of the quality scores of the variant calls
plot <- ggplot(data, aes(x = 1:length(QUAL), y = QUAL)) +
  geom_point() +
  labs(title = "Quality Scores Scatterplot",
       x = "Variant Index",
       y = "Quality Score")
ggsave(output_file, plot)
