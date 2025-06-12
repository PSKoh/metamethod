#### START OF CODE #####

### Set Up --------------

# Set working directory to that of script's current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# R version 4.5.0

# Install packages (if not already installed)
install.packages("metafor")
install.packages("lmerTest")
install.packages("psych")

# Load packages
library(metafor)  # version 4.8-0
library(lmerTest) # version 3.1-3
library(psych)    # version 2.5.3

# Display settings (to disable scientific notation)
options(scipen = 9999, digits = 4)

# Read in data drawn from Lua et al. (2023)
# Original paper: https://doi.org/10.1007/s11031-023-10047-w
mlmmeta_raw = read.csv("NFCWB.csv")

### Prepare Data --------------

# Clean data file (reverse correlation for negative well-being)
mlmmeta_raw$corr_nfcwb = with(mlmmeta_raw, ifelse(wellbeing_category == "Negative well-being", -corr_nfcwb, corr_nfcwb))

# Compute effect sizes for each study
mlmmeta = escalc(
  # Type of effect size measure
  measure = "ZCOR",
  
  # Column for raw correlation coefficients
  ri = corr_nfcwb,
  
  # Column for sample sizes
  ni = sample_size,
  
  # Specify data.frame that the information will be extracted from
  data = mlmmeta_raw
  )

# Categorise publication type into "published" and "unpublished"
# Published: Journal articles
# Unpublished: Conference, Panel Data, Thesis/dissertation, Unpublished data
mlmmeta$publication_type = ifelse(
  mlmmeta$publication_type == "Journal article",
  "Published",
  "Unpublished")

# Order the data frame based on publication type and effect sizes (yi)
mlmmeta = mlmmeta[order(mlmmeta$publication_type, mlmmeta$yi), ]

### Compute Overall Effect Size ------------------- 

# Estimate the overall effect size using the rma.mv() function
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

# Convert from Fisher's Z to Pearson r
mlmmetaresults$b |> fisherz2r() 

### Forest Plot --------------

# Save the forest plot as a PDF file
# Name the pdf file of the forest plot
# cairo_pdf function used for font compatibility
# Adjust the width and height of the pdf file
cairo_pdf(file = "NFCWBforestplot.pdf", width = 14, height = 35)

# Start creating the forest plot itself
# Specify dataset
forest(
  mlmmetaresults, 
  
  # Arrangement of studies
  order = "obs",
  
  # Add y-axis limits
  ylim = c(-3, 111),
  
  # Add sample size information for need for cognition and well-being group
  # Values indicate the x-axis position of the sample size columns  
  ilab = sample_size,
  ilab.xpos = -3,
  
  # Label studies on the forest plot
  slab = paste(author, year, sep = ", "),
  
  # Add x-axis limits
  xlim = c(-5, 3),
  
  # Add confidence interval limits
  # Adjust intervals based on the number of steps
  alim = c(-1.5, 1.5),
  steps = 7,
  
  # Change size of effect size polygons
  efac = 0.3,
  
  # Show (TRUE) or hide (FALSE) default headers
  # Hide when we want to manually specify our own headers
  header = FALSE,
  
  # Add label for confidence interval, in this case, "Fisher's Z"
  xlab = "Fisher's Z"
)

# For the following lines of code,
# Use text function to manually include text within the plot

# Add "Author(s) Year" header
text(x = -4.6, y = 110, "Author(s) Year", font = 2)

# Add "Sample Size" header
text(x = -3, y = 110, "Sample Size", font = 2)

# Add "r [95% CI]" header
text(x = 2.7, y = 110, "Z [95% CI]", font = 2) 

# Close the forest plot and finalise it as a saved file
dev.off()

### Tests for Publication Bias --------------------

# Funnel Plot #

# Save the funnel plot as a PDF file
# Name the pdf file of the funnel plot
# Adjust the width and height of the pdf file
pdf(file = "NFCWBfunnelplot.pdf", width = 8, height = 5)
# par function to adjust margins of the funnel plot
# mar = c(bottom, left, top, right)
par(mar = c(4, 4, 0.3, 1))
# funnel function to create the funnel plot, specify the data to create the plot
funnel(mlmmetaresults, legend = TRUE, xlab = "Fisher's Z")
# Close the funnel plot and finalise it as a saved file
dev.off()

# Rank Correlation Test #
ranktest(mlmmetaresults)

# Eggers' Test # 
lmer(
  # g weighted by SE is predicted by intercept and inverse SE
  # with random intercept by sample
  I(yi / vi) ~ 1 + I(1 / vi) + (1 | sample_id),
  data = mlmmeta
) |>
  # Estimate of interest is the intercept
  summary(correlation = TRUE)

### Moderation Analysis --------------

# Categorical Variable (i.e., publication type)
rma.mv(
  yi = yi,
  V = vi,
  random = ~ 1 | sample_id/meta_id,
  # Specify categorical moderator (i.e., Published articles)
  subset = (publication_type == "Published"),
  data = mlmmeta,
  # To address convergence issues (if it exists)
  control=list(rel.tol=1e-8) 
)

rma.mv(
  yi = yi,
  V = vi,
  random = ~ 1 | sample_id/meta_id,
  # Specify categorical moderator (i.e., Unpublished articles)
  subset = (publication_type == "Unpublished"),
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
       mods = ~ female_proportion,
       method = "REML",
       data = mlmmeta
) |>
       summary()

### Forest Plot of Moderators --------------

# Save the forest plot as a PDF file
# Name the pdf file of the forest plot
# cairo_pdf function used for font compatibility
# Adjust the width and height of the pdf file
cairo_pdf(file = "NFCWBforestplotwithmod.pdf", width = 13, height = 35) 

# Start creating the forest plot itself
# Specify dataset
forest(
  mlmmetaresults,
  # Manually arrange effect sizes by publication type
  # - Unpublished: Rows 108 to 51
  # - Published: Rows 50 to 48
  # The arrangement must consider spacing and must end at row 2
  rows = c(112:40, 36:2),
  
  # Add y-axis limits
  ylim = c(-3, 116),
  
  # Add sample size information for need for cognition and well-being group
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
  steps = 7,
  
  # Change size of effect size polygons
  efac = 0.3,
  
  # Remove headers (if any), for manual input
  header = FALSE,
  
  # Add label for confidence interval, in this case, "Fisher's Z"
  xlab = "Fisher's Z"
)

# For the following lines of code,
# Use text function to manually include text within the plot

# Add text labels for moderator (type of publication)
# x values to indicate the horizontal arrangement of the text
# Labels for different publication types (Moderator Analysis)
# y values indicate the vertical arrangement of the text
# - "Unpublished" (Unpublished data, Panel Data, Thesis/Dissertations) at y = 37
# - "Published" (Journal Articles, Conference) at y = 113
# `pos = 4` ensures left alignment of the text
text(
  x = -7,
  y = c(37, 113),
  pos = 4,
  c("Unpublished", "Published"),
  font = 2
)

# Moderation analysis
res.p = rma.mv(
  yi,
  vi,
  random = ~ 1 | sample_id/meta_id,
  subset = (publication_type == "Published"),
  data = mlmmeta,
  # To address convergence issues (if it exists)
  control=list(rel.tol=1e-8) 
)
res.u = rma.mv(
  yi,
  vi,
  random = ~ 1 | sample_id/meta_id,
  subset = (publication_type == "Unpublished"),
  data = mlmmeta,
  # To address convergence issues (if it exists)
  control = list(rel.tol=1e-8) 
)

# Add summary effect sizes for each of the moderators
addpoly(res.u, row = 1) # summary effect for "Unpublished" group
addpoly(res.p, row = 39) # summary effect for "Published" group

# Add"Author(s) Year" header
text(x = -6.4, y = 115,"Author(s) Year", font = 2) 

# Add “Sample Size” header
text(x = -4.0, y = 115,"Sample Size", font = 2)

# Add "r [95% CI]" header
text(x = 3.6, y = 115, "Z [95% CI]", font = 2) 

# Close the forest plot and finalise it as a saved file
dev.off()

#### END OF CODE ####
