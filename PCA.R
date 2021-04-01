################################################################################
# Install and load libraries
################################################################################

# Uncomment the lines below to install any missing packages
#install.packages("tidyverse")
#install.packages("ggthemes")
#install.packages("ggrepel")
#install.packages("here")

# Load
library(tidyverse)
library(ggthemes)
library(ggrepel)
library(here)

################################################################################
# Read in and process our data
################################################################################

# PCA does not supposed missing values. You must either impute them before hand
# or remove any samples with missing values. This has already been done for our
# data.
data <- read_tsv(here("example data", "proteomics.tsv"))

# Our annotations describing the conditions and replicates are stored in a
# separate file. We'll change colors and shapes based on this information.
sample_annotations <- read_csv(here("example data", "Experimental Design Annotations.csv"))

# Our data contains a column for gene names and protein IDs
# The PCA function requires a numerical matrix so we must remove them
data_numerical <- data %>% select(-`Gene names`, -`Protein IDs`)

# Each sample should be a separate row
# Each feature should be a separate column
#     (e.g. the abunance of a particular protein)
# Our data is in the opposite orientation, so we transpose it
data_transposed <- t(data_numerical)

################################################################################
# Principal Component Analysis (PCA)
################################################################################

# Compute the PCA
pca <- prcomp(data_transposed, center=T, scale=F)

# Compute the variance explained by each principal component
eigen <- pca$sdev^2
variance <- (eigen/sum(eigen)) * 100

# Round to 2 digits after the decimal for better aesthetics
variance <- format(round(variance, 2), nsmall = 2)


################################################################################
# Make the plot
################################################################################

# 1) Append the annotation columns to our PCA results. They are the same size
# and in the same order so we can just concatenate them.
# 2) Create a new Label column that will equal the condition column, but only
# for the first replicate. The other samples won't have a label so we don't
# crowd our figure with text.
plotting.data <- bind_cols(as_tibble(pca$x), sample_annotations) %>%
  mutate(Label = ifelse(Replicate == 1, Condition, ""))


# Pick components to plot
pc_x_axis = 1
pc_y_axis = 2

# Determine column names for the corresponding principal components
pc_x_axis_column = str_c("PC", pc_x_axis)
pc_y_axis_column = str_c("PC", pc_y_axis)

# Create figure panel
panel <- (ggplot(plotting.data, aes(x=get(pc_x_axis_column),
                                    y=get(pc_y_axis_column),
                                    shape=as.factor(Infected),
                                    color=Time,
                                    label=Label))
          # Geometries
          + geom_point(size=2, alpha=0.8)
          + geom_label_repel(size=2, min.segment.length = Inf, fill=alpha(c("white"), 0.25))

          # Scales
          + scale_x_continuous(expand = c(0.1, 0))
          + scale_y_continuous(expand = c(0.1, 0))
          + scale_shape_discrete(name = "Condition")
          + scale_color_discrete(guide = 'none')

          # Labels
          + xlab(str_c(pc_x_axis_column, " (", variance[pc_x_axis], "% explained variance)", sep=""))
          + ylab(str_c(pc_y_axis_column, " (", variance[pc_y_axis], "% explained variance)", sep=""))
          + ggtitle(label = "Principal Component Analysis (PCA)",
                    subtitle = str_c("n =",
                                     prettyNum(nrow(data_numerical), big.mark=","),
                                     "proteins",
                                     sep=" "))

          # Themes
          + theme_clean(base_size = 7)
          + theme(panel.grid.major.y = element_blank(),
                  legend.background = element_rect(fill=grey(0.95), color=NA),
                  legend.title = element_text(size=6),
                  legend.text = element_text(size=6),
                  legend.key = element_blank(),
                  legend.key.size = unit(3.0, 'mm'),
                  legend.position = c(0.88, 0.15),
                  legend.margin = margin(3,3,3,3),
                  plot.title = element_text(hjust=0.5),
                  plot.subtitle = element_text(hjust=0.5, face="italic"))
)

print(panel)

ggsave(here("PCA.pdf"), width=2.5, height=2.5, compress=F, panel)
