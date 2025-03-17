### Set Up --------------

# Set working directory to that of script's current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### read in data
mlmmeta = read.csv("spc_clean_V2.csv")

### Coding for Multilevel Meta-Analysis ------------------
multimmeta = escalc(
  # type of effect size measure
  measure = "SMD",
  
  # columns for sample of each group
  n1i = n_p, n2i = n_a,
  
  # columns for effect size of each group
  m1i = cog_M_p, m2i = cog_M_a,
  
  # columns for standard deviation of each group
  sd1i = cog_SD_p, sd2i = cog_SD_a,
  
  mlmmeta$ID <- seq_along(mlmmeta[,1]),
  
  # specify data file that the information will be extracted from
  data = mlmmeta) 

## calculates the overall effect size using the effect size and sampling variances
# using simple meta code
usingsimplemetacode = rma(
  yi = yi, vi = vi,
  data = multimmeta)
summary(usingsimplemetacode)

# using mlm code
mlmmetaresults = rma.mv(
  yi = yi, V = vi,
  # include random effects for grouping variable (i.e., lab)
  random = ~ 1 | lab_id/unique_id,
  data = multimmeta)
summary(mlmmetaresults)

###### Forest Plot ###### 
pdf(file = "mlmforestplot.pdf", width = 7, height = 14)
forest(mlmmetaresults, 
       # label left column
       header="Author(s), (Year)",
       
       # arrange effect sizes in ascending order
       order = "obs",
       
       # add sample size info into forest plot and specify location of info
       ilab = cbind(n_p, n_a), ilab.xpos = c(-4, -3),
       
       # extracts info from the "author" column of data 
       slab = paste(author, year_published, sep=", "),
       cex=0.5)
# specify location of text
text(-3.5, mlmmetaresults$k+3, "Sample Size", cex=0.5, font=2)
text(c(-4, -3), mlmmetaresults$k+2.02, c("Presence", "Absence"), cex=0.5, font=2)
dev.off()

###### Checking for Publication Bias ###### 
# Funnel Plot
pdf(file = "mlmfunnelplot.pdf", width = 8, height = 5)
funnel(mlmmetaresults,
       legend = TRUE)
dev.off()

# Rank Correlation Test
ranktest(mlmmetaresults)

# Egger' Test
multimmeta$sei_corrected = with(multimmeta, sqrt((n_p+n_a)/(n_p*n_a)))
lmer(
  I(yi/sei_corrected) ~ 1 + I(1/sei_corrected) + (1 | lab_id),
  data = multimmeta
) |> 
  summary(correlation = FALSE)

### Test of Moderators (i.e., type of research design and female proportion) --------------
##### Categorical 
rma.mv(yi = yi, V = vi, 
       random = ~ 1 | lab_id/unique_id,
       # specify categorical moderator (i.e., journal article)
       subset=(publication == "Journal article"),
       data = multimmeta)
rma.mv(yi = yi, V = vi, 
       random = ~ 1 | lab_id/unique_id,
       subset=(publication == "Thesis/dissertation"),
       data = multimmeta)
rma.mv(yi = yi, V = vi, 
       random = ~ 1 | lab_id/unique_id,
       subset=(publication == "Conference"),
       data = multimmeta)

##### Continuous
rma.mv(yi = yi, V = vi,
       random = ~ 1 | lab_id/unique_id,
       # specify categorical moderator (i.e., sex)
       mods =~ female_proportion,
       method = "REML", data = multimmeta) |>
  summary()
