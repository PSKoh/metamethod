#### START OF CODE #####

### Set Up --------------

# Set working directory to that of script's current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# R version 4.5.0

# Install packages (if not already installed)
install.packages("metafor")
install.packages("lmerTest")
install.packages("dplyr")

# Load packages
library(metafor)  # version 4.8-0
library(lmerTest) # version 3.1-3
library(dplyr)    # version 1.1.4

# Display settings (to disable scientific notation)
options(scipen = 9999, digits = 4)

# Read in data drawn from Lua et al. (2023)
# Original paper: https://doi.org/10.1007/s11031-023-10047-w
mlmmeta_raw = read.csv("NFCWB.csv")

### Prepare Data --------------

# Clean data file (reverse correlation for negative well-being)
mlmmeta_new = mlmmeta_raw %>% 
  mutate(corr_nfcwb = ifelse(wellbeing_category == "Negative well-being", -corr_nfcwb, corr_nfcwb)) 

# Compute effect sizes for each study
mlmmeta = escalc(
  # Type of effect size measure
  measure = "ZCOR",
  
  # Column for raw correlation coefficients
  ri = corr_nfcwb,
  
  # Column for sample sizes
  ni = sample_size,
  
  # Specify data.frame that the information will be extracted from
  data = mlmmeta_new
  )

# Convert Publication to a factor with specified levels
mlmmeta$publication_type = factor(
  mlmmeta$publication_type,
  levels = c("Journal article", "Conference", "Panel data", "Thesis/dissertation", "Unpublished data")
)

### Compute Overall Effect Size -------------------  
# Effect size estimates
mlmmetaresults = rma.mv(
  # Effect size estimates
  yi = yi,
  # Sampling variances
  V = vi,
  # Include random effects for grouping variable (i.e., sample)
  random = ~ 1 | sample_id/meta_id,
  # Specify where to get the data from
  data = mlmmeta
)

# summary function used to provide detailed results of the meta-analysis
summary(mlmmetaresults)

### Forest Plot --------------

# Save the forest plot as a PDF file
# Name the pdf file of the forest plot
pdf(file = "NFCWBforestplot.pdf", width = 15, height = 40)

# Start creating the forest plot itself
# Specify dataset
forest(
  mlmmetaresults, 
  
  # Arrangement of studies
  order = "obs",
  
  # Add y-axis limits
  ylim = c(-3, 111),
  
  # Add sample size information for presence and absence of smartphones group
  # Values indicate the x-axis position of the sample size columns  
  ilab = sample_size,
  ilab.xpos = -3,
  
  # Label studies on the forest plot
  slab = paste(author, year, sep = ", "),
  
  # Add x-axis limits
  xlim = c(-5, 3),
  
  # Add confidence interval limits
  # Adjust intervals based on the number of steps
  alim = c(-2, 2),
  steps = 9,
  
  # Change size of effect size polygons
  efac = 0.3,
  
  # Show (TRUE) or hide (FALSE) default headers
  # Hide when we want to manually specify our own headers
  header = FALSE,
  
  # Add label for confidence interval, in this case, "Pearson's corr"
  xlab = "Pearson's corr"
)

# For the following lines of code,
# Use text function to manually include text within the plot

# Add "Author(s) Year" header
text(x = -4.6, y = 110, "Author(s) Year", font = 2)

# Add "Sample Size" header
text(x = -3, y = 110, "Sample Size", font = 2)

# Add "r [95% CI]" header
text(x = 2.7, y = 110, "r [95% CI]", font = 2) 

# Close the forest plot and finalise it as a saved file
dev.off()

### Tests for Publication Bias --------------------

# Funnel Plot

# Save the funnel plot as a PDF file
# Name the pdf file of the funnel plot
# Adjust the width and height of the pdf file
pdf(file = "mlmfunnelplot.pdf", width = 8, height = 5)
# funnel function to create the funnel plot, specify the data to create the plot
funnel(mlmmetaresults, legend = TRUE, xlab = "Pearson's corr")
# Close the funnel plot and finalise it as a saved file
dev.off()

# Rank Correlation Test

ranktest(mlmmetaresults)


# Eggers' Test

# Calculate standard error (SE)
mlmmeta$sei = with(
  mlmmeta, 
  (sqrt(vi) / sqrt(sample_size)))
lmer(
  # g weighted by SE is predicted by intercept and inverse SE
  # with random intercept by sample
  I(yi / sei) ~ 1 + I(1 / sei) + (1 | sample_id),
  data = mlmmeta
) |>
  # Estimate of interest is the slope
  summary(correlation = FALSE)

## Test of Moderators --------------
# Panel data was not included as there was only one study 

# Categorical Variable (i.e., publication type)
rma.mv(
  yi = yi,
  V = vi,
  random = ~ 1 | sample_id/meta_id,
  # Specify categorical moderator (i.e., Journal Article)
  subset = (publication_type == "Journal article"),
  data = mlmmeta,
  # To address convergence issues (if it exists)
  control=list(rel.tol=1e-8) 
)

rma.mv(
  yi = yi,
  V = vi,
  random = ~ 1 | sample_id/meta_id,
  # Specify categorical moderator (i.e., Conference)
  subset = (publication_type == "Conference"),
  data = mlmmeta,
  # To address convergence issues (if it exists)
  control=list(rel.tol=1e-8) 
)

rma.mv(
  yi = yi,
  V = vi,
  random = ~ 1 | sample_id/meta_id,
  # Specify categorical moderator (i.e., Thesis/dissertation)
  subset = (publication_type == "Thesis/dissertation"),
  data = mlmmeta,
  # To address convergence issues (if it exists)
  control=list(rel.tol=1e-8) 
)

rma.mv(
  yi = yi,
  V = vi,
  random = ~ 1 | sample_id/meta_id,
  # Specify categorical moderator (i.e., Unpublished data)
  subset = (publication_type == "Unpublished data"),
  data = mlmmeta,
  # To address convergence issues (if it exists)
  control=list(rel.tol=1e-8) 
)

# Continuous variable (i.e., female proportion)
rma.mv(
  yi = yi,
  V = vi,
  random = ~ 1 | sample_id/meta_id,
  # Specify continuous moderator (i.e., female proportion)
  mods = ~female_proportion,
  method = "REML",
  data = mlmmeta
) |>
  summary()

### Forest Plot of Moderators --------------

# Save the forest plot as a PDF file

### Forest Plot of Moderators --------------

# Save the forest plot as a PDF file
# Name the pdf file of the forest plot
# Adjust the width and height of the pdf file
pdf(file = "NFCWBforestplotwithmod.pdf", width = 15, height = 45) 
forest(
  mlmmetaresults,
  
  # Manually arrange effect sizes by creativity measure type
  # - Journal article: Rows 108 to 51
  # - Conference: Rows 50 to 48
  # - Panel data: Rows 47
  # - Thesis/dissertations: Rows 46 to 25
  # - Unpublished data: Rows 25 to 2
  # The arrangement must consider spacing and must end at row 2
  rows = c(108:51, 50:48,, 47, 46:25, 25:2),
  
  # Add y-axis limits
  ylim = c(-3, 147),
  
  # Add sample size information for presence and absence of smartphones group
  # Values indicate the x-axis position of the sample size columns  
  ilab = sample_size,
  ilab.xpos = -4.2, 
  
  # Label studies on the forest plot
  slab = paste(author, year, sep = ", "),
  
  # Add x-axis limits
  xlim = c(-7, 4),
  
  # Add confidence interval limits
  # Adjust intervals based on the number of steps
  alim = c(-1.5, 1.5),
  steps = 11,
  
  # Change size of effect size polygons
  efac = 0.3,
  
  # Remove headers (if any), for manual input
  header = FALSE,
  
  # Add label for confidence interval, in this case, "Pearson's corr"
  xlab = "Pearson's corr"
)

# For the following lines of code,
# Use text function to manually include text within the plot

# Add text labels for moderator (type of publication)
# x values to indicate the horizontal arrangement of the text
# Labels for different creativity task types (Moderator Analysis)
# y values indicate the vertical arrangement of the text
# - "Journal article" at y = 7
# - "Thesis/Dissertations" at y = 15
# - "Non-Conference" at y = 21
# `pos = 4` ensures left alignment of the text
text(
  x = -7,
  pos = 4,
  y = c(4, 76, 144),
  c("Conference", "Thesis/dissertation", "Journal article"),
  font = 2
)

# Moderation analysis
res.j = rma(
  yi,
  vi,
  random = ~ 1 | sample_id/meta_id,
  subset = (publication_type == "Journal article"),
  data = mlmmeta,
  # To address convergence issues (if it exists)
  control=list(rel.tol=1e-8) 
)
res.t = rma(
  yi,
  vi,
  random = ~ 1 | sample_id/meta_id,
  subset = (publication_type == "Thesis/dissertation"),
  data = mlmmeta
)
res.c = rma(
  yi, 
  vi, 
  random = ~ 1 | sample_id/meta_id,
  subset = (publication_type == "Conference"), 
  data = mlmmeta
)

# Add summary effect sizes for each of the moderators
addpoly(res.c, row = 1) # summary effect for "Conference" group
addpoly(res.t, row = 6) # summary effect for "Thesis/Dissertation" group
addpoly(res.j, row = 78) # summary effect for "Journal article" group

# Add"Author(s) Year" header
text(x = -6.5, y = 146, "Author(s) Year", font = 2) 

# Add “Sample Size” header
text(x = -3.9, y = 146, "Sample Size", font = 2)


# Add "g [95% CI]" header
text(x = 3.6, y = 146, "g [95% CI]", font = 2) 

# Close the forest plot and finalise it as a saved file
dev.off()

#### END OF CODE ####