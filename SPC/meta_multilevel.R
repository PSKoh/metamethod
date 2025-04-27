#### START OF CODE ####

### Set Up --------------

# Set working directory to that of script's current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# R version 4.5.0

# Install packages (if not already installed)
install.packages("metafor")
install.packages("lmerTest")

# Load packages
library(metafor) # version 4.8-0
library(lmerTest) # version 3.1-3

# Display settings (to disable scientific notation)
options(scipen = 9999, digits = 4)

# Read in data drawn from Hartanto et al. 2024
# Original paper: https://doi.org/10.1037/tmb0000123
multilevelmeta_raw = read.csv("SPC.csv")

# Create new column with unique IDs
multilevelmeta_raw$ID = 1:nrow(multilevelmeta_raw)

### Prepare Data --------------

multilevelmeta = escalc(
       # Type of effect size measure
       measure = "SMD",

       # Columns for sample size of each group
       n1i = n_p,
       n2i = n_a,

       # Columns for means of each group
       m1i = cog_M_p,
       m2i = cog_M_a,

       # Columns for standard deviation of each group
       sd1i = cog_SD_p,
       sd2i = cog_SD_a,

       # Specify data.frame that the information will be extracted from
       data = multilevelmeta_raw
)

# Convert Publication to a factor with specified levels
multilevelmeta$publication = factor(
       multilevelmeta$publication,
       levels = c("Journal article", "Thesis/dissertation", "Conference")
)
# Order the data frame based on publication and effect sizes (yi)
multilevelmeta = multilevelmeta[order(multilevelmeta$publication, multilevelmeta$yi), ]

### Compute Overall Effect size --------------

mlmmetaresults = rma.mv(
       # Effect size estimates
       yi = yi,
       # Sampling variances
       V = vi,
       # Include random effects for grouping variable (i.e., lab)
       random = ~ 1 | lab_id / ID,
       # Specify where to get the data from
       data = multilevelmeta
)

# summary function used to provide detailed results of the meta-analysis
summary(mlmmetaresults)

### Forest Plot --------------

# Save the forest plot as a PDF file
# Name the pdf file of the forest plot
pdf(file = "mlmforestplot.pdf", width = 15, height = 40)

# Start creating the forest plot itself
# Specify dataset
forest(
       mlmmetaresults, 

       # Arrangement of studies
       order = "obs",

       # Add y-axis limits
       ylim = c(-3, 140),

       # Add sample size information for presence and absence of smartphones group
       # Values indicate the x-axis position of the sample size columns  
       # -3 for presence of smartphones (n_p)
       # -2.55 for absence of smartphones (n_a)
       ilab = cbind(n_p, n_a),
       ilab.xpos = c(-3, -2.55),

       # Label studies on the forest plot
       slab = paste(author, year_published, sep = ", "),

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

       # Add label for confidence interval, in this case, "Hedge's g"
       xlab = "Hedge's g"
)

# For the following lines of code,
# Use text function to manually include text within the plot

# Add "Author(s) Year" header
text(x = -4.6, y = 139, "Author(s) Year", font = 2)

# Add "Sample Size" header
text(x = -2.8, y = 139.6, "Sample Size", font = 2)

# Add specific sample size column headers, “Presence” and “Absence”
# x values indicate the horizontal arrangement of the columns
# x = -3 for for presence of smartphones
# x = -2.55 for absence of smartphones
# y values indicate the vertical arrangement of the columns
# y = 139 for both
text(c(x = -3, x = -2.55), y = 139, c("Presence", "Absence"), font = 2)

# Add "g [95% CI]" header
text(x = 2.7, y = 139, "g [95% CI]", font = 2) 

# Close the forest plot and finalise it as a saved file
dev.off()

### Tests for Publication Bias --------------

# Funnel Plot #

# Save the funnel plot as a PDF file
# Name the pdf file of the funnel plot
# Adjust the width and height of the pdf file
pdf(file = "mlmfunnelplot.pdf", width = 8, height = 5)
# funnel function to create the funnel plot, specify the data to create the plot
funnel(mlmmetaresults, legend = TRUE, xlab = "Hedge's g")
# Close the funnel plot and finalise it as a saved file
dev.off()

# Egger's Test #

# Calculate standard error (SE)
multilevelmeta$sei_corrected = with(
       multilevelmeta, 
       sqrt((n_p + n_a) / (n_p * n_a)))

# metafor::rma.mv does not have a weights argument
# metafor::regtest does not support rma.mv objects
# For three (or more) level meta-analysis, use lmerTest::lmer instead
lmer(
       # g weighted by SE is predicted by intercept and inverse SE
       # with random intercept by sample
       I(yi / sei_corrected) ~ 1 + I(1 / sei_corrected) + (1 | lab_id),
       data = multilevelmeta
) |>
       # Estimate of interest is the intercept
       summary(correlation = FALSE)

### Moderation Analysis --------------

# Categorical variable (i.e., publication)
rma.mv(
       yi = yi,
       V = vi,
       random = ~ 1 | lab_id / ID,
       # Specify categorical moderator (i.e., Journal Article)
       subset = (publication == "Journal article"),
       data = multilevelmeta,
       # To address convergence issues (if it exists)
       control=list(rel.tol=1e-8) 
)
rma.mv(
       yi = yi,
       V = vi,
       random = ~ 1 | lab_id / ID,
       # Specify categorical moderator (i.e., Thesis/dissertation)
       subset = (publication == "Thesis/dissertation"),
       data = multilevelmeta
)
rma.mv(
       yi = yi,
       V = vi,
       random = ~ 1 | lab_id / ID,
       # Specify categorical moderator (i.e., Conference)
       subset = (publication == "Conference"),
       data = multilevelmeta
)

# Continuous variable (i.e., female proportion)
rma.mv(
       yi = yi,
       V = vi,
       random = ~ 1 | lab_id / ID,
       # Specify continuous moderator (i.e., female proportion)
       mods = ~female_proportion,
       method = "REML",
       data = multilevelmeta
) |>
       summary()

### Forest Plot of Moderators --------------

# Save the forest plot as a PDF file
# Name the pdf file of the forest plot
# Adjust the width and height of the pdf file
pdf(file = "mlmforestplotwithmoderators.pdf", width = 15, height = 45) 
forest(
       mlmmetaresults,

       # Manually arrange effect sizes by creativity measure type
       # - Journal article: Rows 143 to 79
       # - Thesis/Dissertations: Rows 75 to 7
       # - Conference: Rows 3 to 2
       # The arrangement must consider spacing and must end at row 2
       rows = c(143:79, 75:7, 3:2),

       # Add y-axis limits
       ylim = c(-3, 147),

       # Add sample size information for presence and absence of smartphones group
       # Values indicate the x-axis position of the sample size columns  
       # -4.2 for presence of smartphones (n_p)
       # -3.6 for absence of smartphones (n_a)
       ilab = cbind(n_p, n_a),
       ilab.xpos = c(-4.2, -3.6),

       # Label studies on the forest plot
       slab = paste(author, year_published, sep = ", "),

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

       # Add label for confidence interval, in this case, "Hedge's g"
       xlab = "Hedge's g"
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
       y = c(4, 76, 144),
       pos = 4,
       c("Conference", "Thesis/dissertation", "Journal article"),
       font = 2
)

# Moderation analysis
res.j = rma.mv(
       yi,
       vi,
       random = ~ 1 | lab_id / ID,
       subset = (publication == "Journal article"),
       data = multilevelmeta,
       # To address convergence issues (if it exists)
       control=list(rel.tol=1e-8) 
)
res.t = rma.mv(
       yi,
       vi,
       random = ~ 1 | lab_id / ID,
       subset = (publication == "Thesis/dissertation"),
       data = multilevelmeta
)
res.c = rma.mv(
       yi, 
       vi, 
       random = ~ 1 | lab_id / ID,
       subset = (publication == "Conference"), 
       data = multilevelmeta
)

# Add summary effect sizes for each of the moderators
addpoly(res.c, row = 1) # summary effect for "Conference" group
addpoly(res.t, row = 6) # summary effect for "Thesis/Dissertation" group
addpoly(res.j, row = 78) # summary effect for "Journal article" group

# Add"Author(s) Year" header
text(x = -6.5, y = 146, "Author(s) Year", font = 2) 

# Add “Sample Size” header
text(x = -3.9, y = 146.7, "Sample Size", font = 2)

# Add specific sample size column headers, “Presence” and “Absence”
# x values indicate the horizontal arrangement of the columns
# x = -4.2 for for presence of smartphones
# x = -3.6 for absence of smartphones
# y values indicate the vertical arrangement of the columns
# y = 146 for both
text(c(x = -4.2, x = -3.6), y = 146, c("Presence", "Absence"), font = 2)

# Add "g [95% CI]" header
text(x = 3.6, y = 146, "g [95% CI]", font = 2) 

# Close the forest plot and finalise it as a saved file
dev.off()

#### END OF CODE ####
