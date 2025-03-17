#### START OF CODE ####

### Set Up --------------

# Set working directory to that of script's current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# R version 4.4.2 
# Load packages
library(metafor)  # version 4.8-0

# Display settings (to disable scientific notation)
options(scipen = 9999, digits = 4)

# Read in data drawn from Majeed et al. 2021
# Original paper: https://doi.org/10.1002/dys.1677
tradmeta_raw = read.csv("DYSCRE META Database Search _ Summary Table.csv")

### Prepare Data --------------

# escalc function is used to compute effect sizes
tradmeta = escalc(  
  # Type of effect size measure
  # See ?escalc for more information
  measure = "SMD",
  
  # Columns for sample size of each group
  n1i = n_dys, n2i = n_control,
  
  # Columns for means of each group
  m1i = Mean_CRE_dys, m2i = Mean_CRE_control,
  
  # Columns for standard deviation of each group
  sd1i = SD_CRE_dys, sd2i = SD_CRE_control,
  
  # Specify data.frame that the information will be extracted from
  data = tradmeta_raw)

# Fix row order of data frame 
# For easier plotting later
# Convert Creativity.Measure_type to a factor with specified levels
tradmeta$Creativity.Measure_type = factor(tradmeta$Creativity.Measure_type, levels = c("Verbal", "Mixed", "Non-verbal"))
# Order data frame based on type of creativity measure 
tradmeta = tradmeta[order(tradmeta$Creativity.Measure_type), ] 

### Calculate overall effect size, using the effect size (yi) and sampling variances (vi)  --------------

# rma function used to estimate the overall effect size
tradmetaresults = rma(   
  # Effect size estimates
  yi = yi, 
  # Sampling variance 
  vi = vi,
  # REML specifies that the Restricted Maximum Likelihood (REML) method is used to estimate heterogeneity
  method = "REML",      
  # Specify where to get the data from
  data = tradmeta)

# summary function used to provide detailed results of the meta-analysis
summary(tradmetaresults) 

### Forest Plot -------------- 

# pdf function starts the graphics device driver to create PDF files
# Name the file of the forest plot 
# Adjust the width and height of the pdf file
pdf(file = "tradforestplot.pdf", width = 11, height = 9) 

# Start creating the forest plot itself
forest(tradmetaresults, # Specify dataset
       
       # Arrangement of studies
       # "obs" to arrange by effect sizes 
       # To organise by column, replace "obs" with the specific column in the data frame
       order = "obs", 
       
       # Add y-axis limits 
       ylim = c(-1, 16),
       
       # Add sample size information for dyslexia (n_dys) and control (n_control) group into forest plot
       # Adjust positioning of sample size information with x (horizontal) function
       # x values represents the x-coordinates where the sample size values of the dyslexia and control groups will be placed
       # cbind function combines the columns indicating the sample size of the groups (dys and control)
       # ilab.xpos specifies the horizontal arrangement of the columns
       ilab = cbind(n_dys, n_control), 
       ilab.xpos = c(x = -3.8, x = -4.2),  
    
       # Label studies on the forest plot 
       # Extracts info from the "Paper" and "Study" column of data
       # slab argument used to label each effect size with its respective study
       # paste function creates the label 
       # sep function used to separate the columns apart with ','
       slab = paste(Paper, paste("Study", Study), sep=", "),
      
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

# Add “Sample Size” header
# Include desired text of header within the double prime symbol ""
# Adjust the position of the header with the x (horizontal) and y (vertical) function
# Adjust font size of header with the font function 
text(x = -4.0, y = 15, "Sample Size", font = 2) # text function includes text within the plot

# Add specific sample size column headers, “Dslx” (Dyslexia Group) and “Ctrl” (Control Group)
# Include desired text of header within the double prime symbol ""
# Adjust the position of the header with the x (horizontal) and y (vertical) function
# x values represents the x-coordinates of where the "Dslx" and "Ctrl" headers will be placed
# Adjust font size of header with the font function 
text(c(x = -3.8, x = -4.2), y = 14.5, c("Dslx", "Ctrl"), font=2) 

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
pdf(file = "tradfunnelplot.pdf", width = 8, height = 5)
# funnel argument to create the funnel plot, specify the data to create the plot 
funnel(tradmetaresults, 
      legend = TRUE) 
# Close the funnel plot and finalise it as a saved file 
dev.off()

# Rank Correlation Test 

# ranktest argument to compute kendall tau value 
ranktest(tradmetaresults) 

# Egger's Test 

# Calculate standard error (SE) 
tradmeta$sei_corrected = with(tradmeta, sqrt((n_dys+n_control)/(n_dys*n_control)))
rma(               
  # Effect size estimates
  yi = yi, 
  # Sampling variance 
  vi = vi,
  # Indicate moderator which is SE/sei_corrected
  mods = ~ sei_corrected,
  # Indicate weight which is inverse SE^2
  weights = 1/sei_corrected^2,
  # Specify dataset 
  data = tradmeta
) |>
  # Estimate of interest is the slope
  summary()

### Test of Moderators --------------

# Continuous variable (i.e., female proportion)
rma(yi = yi, vi = vi,
    # Specify continuous moderator (i.e., sex)
    mods =~ Proportion.of.female,
    method = "REML", data = tradmeta)

# Categorical variable (i.e., type of creativity measure)
rma(yi = yi, vi = vi,
    # Specify categorical moderator (i.e., verbal)
    subset = (Creativity.Measure_type == "Verbal"),
    method = "REML", 
    data = tradmeta)
rma(yi = yi, vi = vi,
    # Specify categorical moderator (i.e., mixed)
    subset = (Creativity.Measure_type == "Mixed"),  
    method = "REML", 
    data = tradmeta)
rma(yi = yi, vi = vi,
    # Specify categorical moderator (i.e., non-verbal)
    subset = (Creativity.Measure_type == "Non-verbal"), 
    method = "REML", 
    data = tradmeta)

### Forest Plot of Moderators -------------- 

# Name the file of the forest plot 
# Adjust the width and height of the pdf file
pdf(file = "tradforestplotwithmoderators.pdf", width = 14, height = 10) 
forest(tradmetaresults,
       
       # Manually arrange effect sizes by creativity measure type
       # - Mixed: Rows 20 to 18
       # - Verbal: Rows 14 to 10
       # - Non-verbal: Rows 6 to 2
       # The arrangement must consider spacing and must end at row 2 
       rows = c(20:18, 14:10, 6:2), 
       
       # Add y-axis limits 
       ylim = c(-1,24),
       
       # Add sample size information for dyslexia (n_dys) and control (n_control) group into forest plot
       ilab = cbind(n_dys, n_control), 
       ilab.xpos = c(x=-3.8,x=-4.2), 
       
       # Label studies on the forest plot 
       slab = paste(Paper, paste("Study", Study), sep=", "), 
       
       # Add x-axis limits
       xlim = c(-7,4), 
       
       # Add confidence interval limits
       # Adjust intervals based on the number of steps
       alim = c(-2.5, 2.5), steps = 11,
       
       # Remove headers (if any), for manual input
       headers = FALSE)

# Add text labels for moderator (type of creativity task)
# Adjust the position of the labels with x (horizontal) and y (vertical) function
# Labels for different creativity task types (Moderator Analysis)
# - "Verbal" at y = 7
# - "Mixed" at y = 15
# - "Non-verbal" at y = 21
# `pos = 4` ensures left alignment of the text
text(x=-7, pos=4, y=c(7, 15, 21), c("Verbal", "Mixed","Non-verbal"), font=2)

# Moderation analysis
# subset argument ensures only the relevant rows are used 
res.v <- rma(yi, vi, subset = (Creativity.Measure_type == "Verbal"), data = tradmeta)
res.n <- rma(yi, vi, subset = (Creativity.Measure_type == "Non-verbal"), data = tradmeta)
res.m <- rma(yi, vi, subset = (Creativity.Measure_type == "Mixed"), data = tradmeta)

# Add summary effect sizes for each of the moderators 
# addpoly argument adds a summary effect size (diamond shape) for each moderator
# row argument places the summary at the corresponding position in the plot
addpoly(res.v, row=1) # summary effect for "verbal" group
addpoly(res.n, row=9) # summary effect for "non-verbal" group
addpoly(res.m, row=17) # summary effect for "mixed" group

# Add"Author(s) Year" header
text(x=-6.5, y=23, "Author(s) Year", font=2) # text function includes text within the plot

# Add “Sample Size” header
text(x=-4.0, y=23.8, "Sample Size", font=2)

# Add specific sample size column headers, “Dslx” (Dyslexia Group) and “Ctrl” (Control Group)
text(x=c(-4.2, -3.8), y=23, c("Dslx","Ctrl"), font=2)

# Add "g [95% CI]" header
text(x=3.5, y=23, "g [95% CI]", font=2) # text function includes text within the plot

# Close the forest plot and finalise it as a saved file
dev.off()

#### END OF CODE ####
