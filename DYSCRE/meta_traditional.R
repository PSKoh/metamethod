#### START OF CODE ####

### Set Up --------------

# R version 4.5.0

# Set working directory to that of script's current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Install packages (if not already installed)
install.packages("metafor")

# Load packages
library(metafor) # version 4.8-0

# Display settings (to disable scientific notation)
options(scipen = 9999, digits = 4)

# Read in data drawn from Majeed et al. 2021
# Original paper: https://doi.org/10.1002/dys.1677
tradmeta_raw = read.csv("DYSCRE.csv")

### Prepare Data --------------

# Compute effect sizes for each study
tradmeta = escalc(
  # Type of effect size measure
  measure = "SMD",

  # Columns for sample size of each group
  n1i = n_dys,
  n2i = n_control,

  # Columns for means of each group
  m1i = Mean_CRE_dys,
  m2i = Mean_CRE_control,

  # Columns for standard deviation of each group
  sd1i = SD_CRE_dys,
  sd2i = SD_CRE_control,

  # Specify data.frame that the information will be extracted from
  data = tradmeta_raw
)

# Convert Creativity.Measure_type to a factor with specified levels
tradmeta$Creativity.Measure_type = factor(
  tradmeta$Creativity.Measure_type,
  levels = c("Verbal", "Mixed", "Non-verbal")
)

# Order the data frame by Creativity.Measure_type and effect sizes (yi)
tradmeta = tradmeta[order(tradmeta$Creativity.Measure_type, tradmeta$yi), ]

### Compute overall effect size --------------

# Estimate the overall effect size using the rma() function
tradmetaresults = rma(
  # Effect size estimates
  yi = yi,
  # Sampling variance
  vi = vi,
  # Specify method to estimate heterogeneity
  method = "REML",
  # Specify where to get the data from
  data = tradmeta
)

# summary function used to provide detailed results of the meta-analysis
summary(tradmetaresults)

### Forest Plot --------------

# Save the forest plot as a PDF file
# Name the pdf file of the forest plot
pdf(file = "tradforestplot.pdf", width = 11, height = 8)

# Start creating the forest plot itself
# Specify dataset
forest(
  tradmetaresults,

  # Arrangement of studies
  order = "obs",

  # Add y-axis limits
  ylim = c(-2, 16),

  # Add sample size information for dyslexia (n_dys) and control (n_control) group
  # Values indicate the x-axis position of the sample size columns  
  # -4.2 for Dslx (Dyslexia Group)
  # -3.5 for Ctrl (Control Group)
  ilab = cbind(n_dys, n_control),
  ilab.xpos = c(-4.2, -3.5),

  # Label studies on the forest plot
  slab = paste(Paper, paste("Study", Study), sep = ", "),

  # Add x-axis limits
  xlim = c(-8, 4),

  # Add confidence interval limits
  # Adjust intervals based on the number of steps
  alim = c(-2.5, 2.5),
  steps = 11,

  # Show (TRUE) or hide (FALSE) default headers
  # Hide when we want to manually specify our own headers
  header = FALSE,
  
  # Add label for confidence interval, in this case, "Hedge's g"
  xlab = "Hedge's g"
)

# For the following lines of code,
# Use text function to manually include text within the plot

# Add "Author(s) Year" header
text(x = -7.2, y = 14.5, "Author(s) Year", font = 2)

# Add "Sample Size" header
text(x = -3.85, y = 15, "Sample Size", font = 2)

# Add specific sample size column headers for dyslexia and control groups
# x values indicate the horizontal arrangement of the columns
# x = -4.2 for Dslx (Dyslexia Group)
# x = -3.5 for Ctrl (Control Group)
# y values indicate the vertical arrangement of the columns
# y = 14.5 for both
text(c(x = -4.2, x = -3.5), y = 14.5, c("Dslx", "Ctrl"), font = 2)

# Add "g [95% CI]" header
text(x = 3.5, y = 14.5, "g [95% CI]", font = 2)

# Close the forest plot and finalise it as a saved file
dev.off()

### Tests for Publication Bias --------------

# Funnel Plot # 

# Save the funnel plot as a PDF file
# Name the pdf file of the funnel plot
# Adjust the width and height of the pdf file
pdf(file = "tradfunnelplot.pdf", width = 8, height = 5)
# funnel function to create the funnel plot, specify the data to create the plot
funnel(tradmetaresults, legend = TRUE, xlab = "Hedge's g")
# Close the funnel plot and finalise it as a saved file
dev.off()

# Rank Correlation Test #
ranktest(tradmetaresults)

# Egger's Test #

# Calculate standard error (SE)
tradmeta$sei_corrected = with(
  tradmeta,
  sqrt((n_dys + n_control) / (n_dys * n_control))
)
rma(
  # Effect size estimates
  yi = yi,
  # Sampling variance
  vi = vi,
  # Indicate moderator which is SE/sei_corrected
  mods = ~sei_corrected,
  # Indicate weight which is inverse SE^2
  weights = 1 / sei_corrected^2,
  # Specify dataset
  data = tradmeta
) |>
  # Estimate of interest is the intercept
  summary()

### Moderation Analysis --------------

# Continuous variable (i.e., female proportion)
rma(
  yi = yi,
  vi = vi,
  # Specify continuous moderator (i.e., sex)
  mods = ~Proportion.of.female,
  method = "REML",
  data = tradmeta
)

# Categorical variable (i.e., type of creativity measure)
rma(
  yi = yi,
  vi = vi,
  # Specify categorical moderator (i.e., verbal)
  subset = (Creativity.Measure_type == "Verbal"),
  method = "REML",
  data = tradmeta
)
rma(
  yi = yi,
  vi = vi,
  # Specify categorical moderator (i.e., mixed)
  subset = (Creativity.Measure_type == "Mixed"),
  method = "REML",
  data = tradmeta
)
rma(
  yi = yi,
  vi = vi,
  # Specify categorical moderator (i.e., non-verbal)
  subset = (Creativity.Measure_type == "Non-verbal"),
  method = "REML",
  data = tradmeta
)

### Forest Plot of Moderators --------------

# Save the forest plot as a PDF file
# Name the pdf file of the forest plot
# Adjust the width and height of the pdf file
pdf(file = "tradforestplotwithmoderators.pdf", width = 11, height = 9)

# Start creating the forest plot itself
# Specify dataset
forest(
  tradmetaresults,

  # Manually arrange effect sizes by creativity measure type
  # - Verbal: Rows 20 to 18
  # - Mixed: Rows 14 to 10
  # - Non-verbal: Rows 6 to 2
  # The arrangement must consider spacing and must end at row 2
  rows = c(20:18, 14:10, 6:2),

  # Add y-axis limits
  ylim = c(-2, 24),

  # Add sample size information for dyslexia (n_dys) and control (n_control) group
  # Values indicate the x-axis position of the sample size columns
  # -3.8 for Dslx (Dyslexia Group)
  # -3.3 for Ctrl (Control Group)
  ilab = cbind(n_dys, n_control),
  ilab.xpos = c(-3.8, -3.3),

  # Label studies on the forest plot
  slab = paste(Paper, paste("Study", Study), sep = ", "),

  # Add x-axis limits
  xlim = c(-7, 4),

  # Add confidence interval limits
  # Adjust intervals based on the number of steps
  alim = c(-1.5, 1.5),
  steps = 7,

  # Remove headers (if any), for manual input
  header = FALSE,

  # Add label for confidence interval, in this case, "Hedge's g" 
  xlab = "Hedge's g"
)
# For the following lines of code,
# Use text function to manually include text within the plot

# Add text labels for moderator (type of creativity task)
# x values to indicate the horizontal arrangement of the text
# Labels for different creativity task types (Moderator Analysis)
# y values indicate the vertical arrangement of the text
# - "Non-verbal" at y = 7
# - "Mixed" at y = 15
# - "Verbal" at y = 21
# `pos = 4` ensures left alignment of the text
text(
  x = -7,
  y = c(7, 15, 21),
  pos = 4,
  c("Non-verbal", "Mixed", "Verbal"),
  font = 2
)

# Moderation analysis
res.v = rma(
  yi,
  vi,
  subset = (Creativity.Measure_type == "Verbal"),
  data = tradmeta
)
res.n = rma(
  yi,
  vi,
  subset = (Creativity.Measure_type == "Non-verbal"),
  data = tradmeta
)
res.m = rma(
  yi,
  vi,
  subset = (Creativity.Measure_type == "Mixed"),
  data = tradmeta
)

# Add summary effect sizes for each of the moderators
addpoly(res.n, row = 1) # summary effect for "non-verbal" group
addpoly(res.m, row = 9) # summary effect for "mixed" group
addpoly(res.v, row = 17) # summary effect for "verbal" group

# Add "Author(s) Year" header
text(x = -6.3, y = 23, "Author(s) Year", font = 2) 
# Add “Sample Size” header
text(x = -3.6, y = 23.7, "Sample Size", font = 2)

# Add specific sample size column headers for dyslexia and control groups
# x values indicate the horizontal arrangement of the columns
# x = -3.8 for Dslx (dyslexia Group)
# x = -3.3 for Ctrl (control Group)
# y values indicate the vertical arrangement of the columns
# y = 23 for both
text(c(x = -3.8, x = -3.3), y = 23, c("Dslx", "Ctrl"), font = 2)

# Add "g [95% CI]" header
text(x = 3.5, y = 23, "g [95% CI]", font = 2) 
# Close the forest plot and finalise it as a saved file
dev.off()

#### END OF CODE ####
