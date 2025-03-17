#### START OF CODE ####

### Set Up --------------

# Set working directory to that of script's current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# R version 4.4.2 
# Load packages
library(metafor)  # version 4.8-0
library(lmer)     # version 1.1-36

# Display settings (to disable scientific notation)
options(scipen = 9999, digits = 4)

# Read in data drawn from Hartanto et al. 2024
# Original paper: https://doi.org/10.1037/tmb0000123
mlmmeta = read.csv("spc_clean_V2.csv")

# Create new column with unique IDs
mlmmeta$ID <- 1:nrow(mlmmeta)

### Prepare Data --------------

# escalc function is used to compute effect sizes
multimmeta = escalc(
  # Type of effect size measure
  # See ?escalc for more information
  measure = "SMD",
  
  # Columns for sample size of each group
  n1i = n_p, n2i = n_a,
  
  # Columns for means of each group
  m1i = cog_M_p, m2i = cog_M_a,
  
  # Columns for standard deviation of each group
  sd1i = cog_SD_p, sd2i = cog_SD_a,
  
  # Specify data.frame that the information will be extracted from
  data = mlmmeta) 

# Fix row order of data frame 
# For easier plotting later
# Convert Publication to a factor with specified levels
mlmmeta$publication = factor(mlmmeta$publication, levels = c("Journal article", "Thesis/dissertation", "Conference"))
# Order data frame based on publication
mlmmeta = mlmmeta[order(mlmmeta$publication), ] 

### Calculate overall effect size, using the effect size (yi) and sampling variances (vi)  --------------
# Using traditional meta code (refer to meta_traditional.R)
usingsimplemetacode = rma(
  yi = yi, vi = vi,
  data = multimmeta)
summary(usingsimplemetacode)

# Using multilevel code
mlmmetaresults = rma.mv(
  # Effect size estimates
  yi = yi, 
  # Sampling variances
  V = vi,
  # Include random effects for grouping variable (i.e., lab)
  random = ~ 1 | lab_id/unique_id,
  # Specify where to get the data from
  data = multimmeta)

# summary function used to provide detailed results of the meta-analysis
summary(mlmmetaresults)

###### Forest Plot ###### ### Forest Plot -------------- 

# pdf function starts the graphics device driver to create PDF files
# Name the file of the forest plot 
# Adjust the width and height of the pdf file
pdf(file = "mlmforestplot.pdf", width = 7, height = 14)

# Start creating the forest plot itself
forest(mlmmetaresults, # Specify dataset
       
       # Arrangement of studies
       # "obs" to arrange by effect sizes
       order = "obs",

       # Add y-axis limits
       ylim = c(-1, 16),
       
       # Add sample size information for dyslexia (n_p) and control (n_a) group into forest plot
       # Adjust positioning of sample size information with x (horizontal) function
       # x values represents the x-coordinates where the sample size values of the experimental and control groups will be placed
       # cbind function combines the columns indicating the sample size of the groups (n_p and n_a)
       # ilab.xpos specifies the horizontal arrangement of the columns
       ilab = cbind(n_p, n_a), 
       ilab.xpos = c(x = -4, z = -3),
       
       # Label studies on the forest plot 
       # Extracts info from the "author" and "year_published" column of data
       # slab argument used to label each effect size with its respective study
       # paste function creates the label 
       # sep function used to separate the columns apart with ','
       slab = paste(author, year_published, sep=", "),
       
       # Add x-axis limits
       xlim = c(-8, 4),
       
       # Add confidence interval limits
       # Adjust intervals based on the number of steps
       alim = c(-2.5, 2.5), steps = 11,
       
       # Show (TRUE) or hide (FALSE) default headers
       # Hide when we want to manually specify our own headers
       headers = FALSE)

# For the following lines of code, 
# use text function to manually include text within the plot

# Add "Author(s) Year" header
# Include desired text of header within the double prime symbol ""
# Adjust the position of the header with the x (horizontal) and y (vertical) function
# Adjust font size of header with the font function 
text(x = -7.2, y = 14.5, "Author(s) Year", font = 2) 

# Add "Sample Size" header
# Include desired text of header within the double prime symbol ""
# Adjust the position of the header with the x (horizontal) and y (vertical) function
# Adjust font size of header with the font function 
text(-3.5, y = 15, "Sample Size", cex=0.5, font=2)

# Add specific sample size column headers, “Presence” and “Absence” 
# Include desired text of header within the double prime symbol ""
# Adjust the position of the header with the x (horizontal) and y (vertical) function
# x values represents the x-coordinates of where the "Presence" and "Absence" headers will be placed
# Adjust font size of header with the font function 
text(c(-4, -3), y = 14.5, c("Presence", "Absence"), font=2)

# Add "g [95% CI]" header
# Include desired text of header within the double prime symbol ""
# Adjust the position of the header with the x (horizontal) and y (vertical) function
# Adjust font size of header with the font function 
text(x = 3.5, y = 14.5, "g [95% CI]", font = 2) # text function includes text within the plot

# Close the forest plot and finalise it as a saved file
dev.off()

### Checking for Publication Bias -------------- 

# Funnel Plot

# pdf function starts the graphics device driver to create PDF files
# Name the file of the funnel
# Adjust the width and height of the pdf file
pdf(file = "mlmfunnelplot.pdf", width = 8, height = 5)
# funnel argument to create the funnel plot, specify the data to create the plot
funnel(mlmmetaresults,
       legend = TRUE)
#Close the funnel plot and finalise it as a saved file
dev.off()

# Rank Correlation Test\

# ranktest argument to compute kendall tau value
ranktest(mlmmetaresults)

# Egger' Test

# Calculate standard error (SE)
multimmeta$sei_corrected = with(multimmeta, sqrt((n_p+n_a)/(n_p*n_a)))
lmer(
  # g weighted by SE is predicted by intercept and inverse SE
  # with random intercept by sample
  I(yi/sei_corrected) ~ 1 + I(1/sei_corrected) + (1 | lab_id),
  data = multimmeta
) |> 
  # estimate of interest is the intercept
  summary(correlation = FALSE)

## Test of Moderators --------------

# Categorical variable (i.e., publication)
rma.mv(yi = yi, V = vi, 
       random = ~ 1 | lab_id/unique_id,
       # Specify categorical moderator (i.e., journal article)
       subset=(publication == "Journal article"),
       data = multimmeta)
rma.mv(yi = yi, V = vi, 
       random = ~ 1 | lab_id/unique_id,
       # Specify categorical moderator (i.e., thesis/dissertation)
       subset=(publication == "Thesis/dissertation"),
       data = multimmeta)
rma.mv(yi = yi, V = vi, 
       random = ~ 1 | lab_id/unique_id,
       # Specify categorical moderator (i.e., conference)
       subset=(publication == "Conference"),
       data = multimmeta)

# Continuous variable (i.e., female proportion)
rma.mv(yi = yi, V = vi,
       random = ~ 1 | lab_id/unique_id,
       # Specify categorical moderator (i.e., sex)
       mods =~ female_proportion,
       method = "REML", data = multimmeta) |>
  summary()
