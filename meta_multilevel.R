#### START OF CODE ####

### Set Up --------------

# Set working directory to that of script's current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# R version 4.4.3

# Install packages (if not already installed)
install.packages("metafor")

# Load packages
library("metafor") # version 4.8-0
library("lmerTest") # version 3.1-3

# Display settings (to disable scientific notation)
options(scipen = 9999, digits = 4)

# Read in data drawn from Hartanto et al. 2024
# Original paper: https://doi.org/10.1037/tmb0000123
mlmmeta = read.csv("SPC.csv")

# Create new column with unique IDs
mlmmeta$ID = 1:nrow(mlmmeta)

### Prepare Data --------------

# escalc function is used to compute effect sizes
multimmeta = escalc(
       # Type of effect size measure
       # See ?escalc for more information
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
       data = mlmmeta
)

# Fix row order of data frame
# For easier plotting later
# Convert Publication to a factor with specified levels
mlmmeta$publication = factor(
       mlmmeta$publication,
       levels = c("Journal article", "Thesis/dissertation", "Conference")
)
# Order data frame based on publication
mlmmeta = mlmmeta[order(mlmmeta$publication), ]

### Calculate overall effect size, using the effect size (yi) and sampling variances (vi)  --------------
# Using traditional meta code (refer to meta_traditional.R), which is not recommended for multilevel meta-analysis
usingsimplemetacode = rma(
       yi = yi,
       vi = vi,
       data = multimmeta
)
summary(usingsimplemetacode)

# Using multilevel code
mlmmetaresults = rma.mv(
       # Effect size estimates
       yi = yi,
       # Sampling variances
       V = vi,
       # Include random effects for grouping variable (i.e., lab)
       random = ~ 1 | lab_id / ID,
       # Specify where to get the data from
       data = multimmeta
)

# summary function used to provide detailed results of the meta-analysis
summary(mlmmetaresults)

###### Forest Plot ###### ### Forest Plot --------------

# pdf function starts the graphics device driver to create PDF files
# Name the file of the forest plot
# Adjust the width and height of the pdf file
# pdf function starts the graphics device driver to create PDF files
# Name the file of the forest plot
# Adjust the width and height of the pdf file
pdf(file = "mlmforestplot.pdf", width = 15, height = 40)

# Start creating the forest plot itself
forest(
       mlmmetaresults, # Specify dataset

       # Arrangement of studies
       # "obs" to arrange by effect sizes
       order = "obs",

       # Add y-axis limits
       ylim = c(-1, 140),

       # Add sample size information for presence of smartphones (n_p) and absence of smartphones(n_a) group into forest plot
       # Adjust positioning of sample size information with x (horizontal) function
       # x values represents the x-coordinates where the sample size values of the experimental and control groups will be placed
       # cbind function combines the columns indicating the sample size of the groups (n_p and n_a)
       # ilab.xpos specifies the horizontal arrangement of the columns
       ilab = cbind(n_p, n_a),
       ilab.xpos = c(x = -3, x = -2.5),

       # Label studies on the forest plot
       # Extracts info from the "author" and "year_published" column of data
       # slab argument used to label each effect size with its respective study
       # paste function creates the label
       # sep function used to separate the columns apart with ','
       slab = paste(author, year_published, sep = ", "),

       # Add x-axis limits
       xlim = c(-5, 3),

       # Add confidence interval limits
       # Adjust intervals based on the number of steps
       alim = c(-2, 2),
       steps = 9,

       # Change size of polygons
       efac = 0.5,

       # Show (TRUE) or hide (FALSE) default headers
       # Hide when we want to manually specify our own headers
       header = FALSE
)

# For the following lines of code,
# use text function to manually include text within the plot

# Add "Author(s) Year" header
# Include desired text of header within the double prime symbol ""
# Adjust the position of the header with the x (horizontal) and y (vertical) function
# Adjust font size of header with the font function
text(x = -4.6, y = 139, "Author(s) Year", font = 2)

# Add "Sample Size" header
# Include desired text of header within the double prime symbol ""
# Adjust the position of the header with the x (horizontal) and y (vertical) function
# Adjust font size of header with the font function
text(-2.8, y = 140, "Sample Size", font = 2)

# Add specific sample size column headers, “Presence” and “Absence”
# Include desired text of header within the double prime symbol ""
# Adjust the position of the header with the x (horizontal) and y (vertical) function
# x values represents the x-coordinates of where the "Presence" and "Absence" headers will be placed
# Adjust font size of header with the font function
text(c(-3, -2.5), y = 139, c("Presence", "Absence"), font = 2)

# Add "g [95% CI]" header
# Include desired text of header within the double prime symbol ""
# Adjust the position of the header with the x (horizontal) and y (vertical) function
# Adjust font size of header with the font function
text(x = 2.7, y = 139, "g [95% CI]", font = 2) # text function includes text within the plot

# Close the forest plot and finalise it as a saved file
dev.off()

### Checking for Publication Bias --------------

# Funnel Plot

# pdf function starts the graphics device driver to create PDF files
# Name the file of the funnel
# Adjust the width and height of the pdf file
pdf(file = "mlmfunnelplot.pdf", width = 8, height = 5)
# funnel argument to create the funnel plot, specify the data to create the plot
funnel(mlmmetaresults, legend = TRUE)
#Close the funnel plot and finalise it as a saved file
dev.off()

# Rank Correlation Test

# ranktest argument to compute kendall tau value
ranktest(mlmmetaresults)

# Egger' Test

# Calculate standard error (SE)
multimmeta$sei_corrected = with(multimmeta, sqrt((n_p + n_a) / (n_p * n_a)))
lmer(
       # g weighted by SE is predicted by intercept and inverse SE
       # with random intercept by sample
       I(yi / sei_corrected) ~ 1 + I(1 / sei_corrected) + (1 | lab_id),
       data = multimmeta
) |>
       # estimate of interest is the intercept
       summary(correlation = FALSE)

## Test of Moderators --------------

# Categorical variable (i.e., publication)
rma.mv(
       yi = yi,
       V = vi,
       random = ~ 1 | lab_id / ID,
       # Specify categorical moderator (i.e., journal article)
       subset = (publication == "Journal article"),
       data = multimmeta
)
rma.mv(
       yi = yi,
       V = vi,
       random = ~ 1 | lab_id / ID,
       # Specify categorical moderator (i.e., thesis/dissertation)
       subset = (publication == "Thesis/dissertation"),
       data = multimmeta
)
rma.mv(
       yi = yi,
       V = vi,
       random = ~ 1 | lab_id / ID,
       # Specify categorical moderator (i.e., conference)
       subset = (publication == "Conference"),
       data = multimmeta
)

# Continuous variable (i.e., female proportion)
rma.mv(
       yi = yi,
       V = vi,
       random = ~ 1 | lab_id / ID,
       # Specify categorical moderator (i.e., sex)
       mods = ~female_proportion,
       method = "REML",
       data = multimmeta
) |>
       summary()

### Forest Plot of Moderators --------------

# pdf function starts the graphics device driver to create PDF files
# Name the file of the forest plot
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
       ylim = c(-1, 147),

       # Add sample size information for presence of smartphone (n_p) and absence of smartphone (n_a) group into forest plot
       ilab = cbind(n_p, n_a),
       ilab.xpos = c(x = -4.2, x = -3.6),

       # Label studies on the forest plot
       slab = paste(author, year_published, sep = ", "),

       # Add x-axis limits
       xlim = c(-7, 4),

       # Add confidence interval limits
       # Adjust intervals based on the number of steps
       alim = c(-2.5, 2.5),
       steps = 11,

       # Change size of polygons
       efac = 0.5,

       # Remove headers (if any), for manual input
       header = FALSE
)

# Add text labels for moderator (type of publication)
# Adjust the position of the labels with x (horizontal) and y (vertical) function
# Labels for different creativity task types (Moderator Analysis)
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
# subset argument ensures only the relevant rows are used
res.j = rma(
       yi,
       vi,
       random = ~ 1 | lab_id / ID,
       subset = (publication == "Journal article"),
       data = multimmeta
)
res.t = rma(
       yi,
       vi,
       random = ~ 1 | lab_id / ID,
       subset = (publication == "Thesis/dissertation"),
       data = multimmeta
)
res.c = rma(
       yi, 
       vi, 
       random = ~ 1 | lab_id / ID,
       subset = (publication == "Conference"), 
       data = multimmeta
)

# Add summary effect sizes for each of the moderators
# addpoly argument adds a summary effect size (diamond shape) for each moderator
# row argument places the summary at the corresponding position in the plot
addpoly(res.c, row = 1) # summary effect for "Conference" group
addpoly(res.t, row = 6) # summary effect for "Thesis/Dissertation" group
addpoly(res.j, row = 78) # summary effect for "Journal article" group

# Add"Author(s) Year" header
text(x = -6.5, y = 146, "Author(s) Year", font = 2) # text function includes text within the plot

# Add “Sample Size” header
text(x = -3.9, y = 147, "Sample Size", font = 2)

# Add specific sample size column headers, “Presence” (Presence of Smartphones Group) and “Absence” (Absence of Smartphones Group)
text(c(x = -4.2, x = -3.6), y = 146, c("Presence", "Absence"), font = 2)

# Add "g [95% CI]" header
text(x = 3.6, y = 146, "g [95% CI]", font = 2) # text function includes text within the plot

# Close the forest plot and finalise it as a saved file
dev.off()

#### END OF CODE ####
