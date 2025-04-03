##########################################################################
##                                                                      ##
##            NEED FOR COGNITION AND WELLBEING META-ANALYSIS            ##
##                                                                      ##
##########################################################################


########## admin ########## 

# set wd to where current R script is located at

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load library

library(metafor)

# read raw data file

d = read.csv("NFCWB_Data_Final_0923.csv")


# clean data file (reverse correlation for negative well-being and calculate sei)

d.new = d %>% dplyr::mutate(
  corr_nfcwb = ifelse(wellbeing_category == "Negative well-being", -corr_nfcwb, corr_nfcwb)) 
d.new = d.new %>% escalc(
  measure = "ZCOR",
  ri = corr_nfcwb,
  ni = sample_size,
  data = .) %>%
  dplyr::mutate(
    sei = sqrt(vi) / sqrt(sample_size)
  )


########## analysis ##########  

# publication bias -----------

# egger's test for publication bias [ include SE for each effect size as moderator ]

EGGERS = rma.mv( yi = yi, V = vi, mods = sei,
                 random = ~ 1 | sample_id/meta_id, 
                 method = "REML",
                 data = d.new, slab = label) # no evidence for pub bias
EGGERS %>% summary()


# published vs unpublished records

rma.mv(yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(publication_type != "Journal article" & publication_type != "Panel data"), slab = label)
psych::fisherz2r(c(0.2180,0.1605,0.2754))

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(publication_type == "Journal article" | publication_type == "Panel data"), slab = label)
psych::fisherz2r(c(0.1886,0.1361,0.2411))


# year of publication as moderator

rma.mv( yi = yi, V = vi, mods = as.numeric(substr(d.new$year, 1,4)),
        random = ~ 1 | sample_id/meta_id,
        method = "REML",
        data = d.new, slab = label) # no evidence for pub bias


# mod = nfc scale used

rma.mv( yi = yi, V = vi,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(nfc_measure == "NCS-34"), slab = label)
psych::fisherz2r(c(0.2083,0.1167,0.2998))

rma.mv( yi = yi, V = vi,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(nfc_measure == "NCS-18"), slab = label)
psych::fisherz2r(c(0.1997,0.1489,0.2505))

rma.mv( yi = yi, V = vi,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(nfc_measure == "NCS-6"), slab = label)
psych::fisherz2r(c(0.2767,0.1868,0.3665))

rma.mv( yi = yi, V = vi,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(nfc_measure == "Other"), slab = label)
psych::fisherz2r(c(0.1704,0.0966,0.2442))


# calculate overall meta analytic effect size ---------------

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new, slab = label)
psych::fisherz2r(c(0.1977,0.1585,0.2369))


# calculate overall meta analytic effect size for subsets of well-being categories ---------------

# Life satisfaction

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Life satisfaction"), slab = label)
psych::fisherz2r(c(0.1434,0.0798,0.2069))
table(d.new %>% filter(wellbeing_indicator == "Life satisfaction") %>% select(wellbeing_measure))


# Positive affect

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Positive affect"), slab = label)
psych::fisherz2r(c(0.3619,0.2999,0.4238))
table(d.new %>% filter(wellbeing_indicator == "Positive affect") %>% select(wellbeing_measure))


# Negative affect

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Negative affect"), slab = label)
psych::fisherz2r(c(0.0894,0.0274,0.1513))
table(d.new %>% filter(wellbeing_indicator == "Negative affect") %>% select(wellbeing_measure))


# Purpose in life

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Purpose in life"), slab = label)
psych::fisherz2r(c(0.2628,0.1561,0.3694))
table(d.new %>% filter(wellbeing_indicator == "Purpose in life") %>% select(wellbeing_measure))


# Self-acceptance

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Self-acceptance"), slab = label)
psych::fisherz2r(c(0.2003,0.1632,0.2373))
table(d.new %>% filter(wellbeing_indicator == "Self-acceptance") %>% select(wellbeing_measure))


# Personal growth

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Personal growth"), slab = label)
psych::fisherz2r(c(0.4902,0.4046,0.5757))
table(d.new %>% filter(wellbeing_indicator == "Personal growth") %>% select(wellbeing_measure))


# Environmental mastery

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Environmental mastery"), slab = label)
psych::fisherz2r(c(0.3462,0.2376,0.4549))
table(d.new %>% filter(wellbeing_indicator == "Environmental mastery") %>% select(wellbeing_measure))


# Positive relations with others

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Positive relations with others"), slab = label)
psych::fisherz2r(c(0.1611,0.0659,0.2564))
table(d.new %>% filter(wellbeing_indicator == "Positive relations with others") %>% select(wellbeing_measure))


# Autonomy

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Autonomy"), slab = label)
psych::fisherz2r(c(0.4088,0.2655,0.5522))
table(d.new %>% filter(wellbeing_indicator == "Autonomy") %>% select(wellbeing_measure))


# Depression

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Depression"), slab = label)
psych::fisherz2r(c(0.1899,0.1256,0.2543))
table(d.new %>% filter(wellbeing_indicator == "Depression") %>% select(wellbeing_measure))


# Anxiety

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Anxiety"), slab = label)
psych::fisherz2r(c(0.1840,0.1078,0.2603))
table(d.new %>% filter(wellbeing_indicator == "Anxiety") %>% select(wellbeing_measure))


# Stress

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Stress"), slab = label)
psych::fisherz2r(c(0.0690,0.0024,0.1355))
table(d.new %>% filter(wellbeing_indicator == "Stress") %>% select(wellbeing_measure))



# test for moderators: well-being valence ---------------

rma.mv( yi = yi, V = vi,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(wellbeing_category == "Positive well-being"), slab = label)
psych::fisherz2r(c(0.2218,0.1737,0.2698))

rma.mv( yi = yi, V = vi,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(wellbeing_category == "Negative well-being"), slab = label)
psych::fisherz2r(c(0.1348,0.0941,0.1755))


# test for moderators: student sample ---------------

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(sample_student == "Yes (Children)" | sample_student == "Yes (Undergrads/ Graduates)"), slab = label)
psych::fisherz2r(c(0.2340,0.1841,0.2840))

rma.mv( yi = yi, V = vi,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(sample_student == "No"), slab = label)
psych::fisherz2r(c(0.1319,0.0700,0.1938))



# test for moderators: student sample [subgroup analyses] --------------

# Life satisfaction

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(sample_student == "Yes (Children)" | sample_student == "Yes (Undergrads/ Graduates)") %>% filter(wellbeing_indicator == "Life satisfaction"), slab = label)
psych::fisherz2r(c(0.1648,0.0612,0.2684))  

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(sample_student == "No") %>% filter(wellbeing_indicator == "Life satisfaction"), slab = label)
psych::fisherz2r(c(0.1039,0.0428,0.1649))


# Anxiety
rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(sample_student == "Yes (Children)" | sample_student == "Yes (Undergrads/ Graduates)") %>% filter(wellbeing_indicator == "Anxiety"), slab = label)
psych::fisherz2r(c(0.2279,0.1643,0.2916))

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(sample_student == "No") %>% filter(wellbeing_indicator == "Anxiety"), slab = label)
psych::fisherz2r(c(0.0450,-0.1439,0.2339))


# Purpose in life
rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(sample_student == "Yes (Children)" | sample_student == "Yes (Undergrads/ Graduates)") %>% filter(wellbeing_indicator == "Purpose in life"), slab = label)
psych::fisherz2r(c(0.3600,0.2349,0.4851))

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% dplyr::filter(sample_student == "No") %>% filter(wellbeing_indicator == "Purpose in life"), slab = label)
psych::fisherz2r(c(0.1292,-0.0087,0.2670))

# test for moderators: age ---------------

rma.mv( yi = yi, V = vi, mods = age_m,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new, slab = label)

table(psych::describe(d.new$age_m))

rma.mv( yi = yi, V = vi, mods = (age_m-28.7899602150538)/13.220699917234,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new, slab = label)

# test for moderators: age [subgroup analyses] ---------------  

# Life satisfaction
rma.mv( yi = yi, V = vi, mods = age_m,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Life satisfaction"), slab = label)

# Anxiety
rma.mv( yi = yi, V = vi, mods = age_m,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Anxiety"), slab = label)

# Purpose in life
rma.mv( yi = yi, V = vi, mods = age_m,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Purpose in life"), slab = label)

# test for moderators: gender ---------------    

rma.mv( yi = yi, V = vi, mods = female_proportion,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new, slab = label)

# test for moderators: gender [subgroup analyses] ---------------  

# Life satisfaction
rma.mv( yi = yi, V = vi, mods = female_proportion,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Life satisfaction"), slab = label)

# Anxiety
rma.mv( yi = yi, V = vi, mods = female_proportion,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Anxiety"), slab = label)

# Purpose in life
rma.mv( yi = yi, V = vi, mods = female_proportion,
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.new %>% filter(wellbeing_indicator == "Purpose in life"), slab = label)

########## sensitivity analyses ##########

d.psq = d.new %>% dplyr::filter(
  inclusion.in.the.sample != "No" & inclusion.in.the.sample != "Unclear",
  study.subjects.and.the.setting != "No" & study.subjects.and.the.setting != "Unclear",
  need.for.cognition.measure != "No" & need.for.cognition.measure != "Unclear",
  outcomes.measure != "No" & outcomes.measure != "Unclear" )

# overall meta-analytic effect size ---------

rma.mv( yi = yi, V = vi, 
        random = ~ 1 | sample_id/meta_id, 
        method = "REML",
        data = d.psq, slab = label)
psych::fisherz2r(c(0.1952,0.1425,0.2479))

########## generate figures ##########  

# eggers' test funnel plot ---------------

EGGERS %>% funnel(xlab = "Need for Cognition - Well-Being Effect Size")

# sub-group forest plot for NFC scale as moderator funnel plot ---------------

extract_subgroup_info = function(lapply_output, rlabs) {
  est = sapply(lapply_output, coef); names(est) = NULL
  se  = sapply(lapply_output, function(x){x$se})
  nk  = sapply(lapply_output, function(x){x$s.nlevels})
  r   = psych::fisherz2r(est) %>% round(2) %>% format(nsmall = 2)
  return(data.frame(rowlabel = rlabs, ri = r, esti = est, sei = se, ni = nk[1,], ki = nk[2,]))
}

sg1 = c("NCS-34", "NCS-18", "NCS-6","Other"); sg1.res = lapply(
  sg1,
  function(x){
    rma.mv(
      yi, vi,
      random = ~ 1 | sample_id/meta_id,
      method = "REML",
      data = d.new, subset = nfc_measure == x)
  }
) %>% extract_subgroup_info(sg1); rm(sg1)

ilab.xpos =  c(-2.5, -2.0, -1.5)

try(dev.off()); pdf(file="NFCWB_nfcsubgroups_forest_210923.pdf", width=7.5, height=4)
leftmost = -5.5
forest(
  x = sg1.res$esti,
  sei = sg1.res$sei,
  slab = sg1.res$rowlabel,
  top = 2, ylim = c(0, 7), rows = c(5:2),
  alim = c(-0.5, 0.5), steps = 3, xlim = c(-5.4, 2.25),
  xlab = "Need for Cognition - Well-Being Effect Size",
  cex = 0.8,
  header = FALSE,
  ilab = sg1.res %>% dplyr::select(ni, ki, ri),
  ilab.xpos =  ilab.xpos,
)
text("Need for Cognition Measure", y = 7, x = leftmost, pos = 4, font = 2)
abline(h = 1)
headertext = c("m", "k", "r"); for(i in 1:length(ilab.xpos)){
  text(ilab.xpos[i], 7, headertext[i], adj = c(0.5, 0.5), font = 2, cex = 0.8)
}; rm(headertext); rm(i)
dev.off(); rm(leftmost)


# sub-group forest plot for wellbeing subcategories funnel plot ---------------

sg0 = c("Life satisfaction", "Positive affect", "Negative affect","Purpose in life", "Self-acceptance", "Personal growth", "Environmental mastery", "Positive relations with others", "Autonomy", "Depression", "Anxiety", "Stress"); sg0.res = lapply(
  sg0,
  function(x){
    rma.mv(
      yi, vi,
      random = ~ 1 | sample_id/meta_id,
      method = "REML",
      data = d.new, subset = wellbeing_indicator == x)
  }
) %>% extract_subgroup_info(sg0); rm(sg0)

ilab.xpos =  c(-2.5, -2.0, -1.5)

try(dev.off()); pdf(file="NFCWB_wbsubgroups_forest_210923.pdf", width=7.5, height=6)
leftmost = -5.4
forest(
  x = sg0.res$esti,
  sei = sg0.res$sei,
  slab = sg0.res$rowlabel,
  top = 2, ylim = c(2, 16), rows = c(14:3),
  alim = c(-0.75, 0.75), steps = 3, xlim = c(-5.4, 2.25),
  xlab = "Need for Cognition - Well-Being Effect Size",
  cex = 0.8,
  header = FALSE,
  ilab = sg0.res %>% dplyr::select(ni, ki, ri),
  ilab.xpos =  ilab.xpos,
)
text("Well-Being Type", y = 16, x = leftmost, pos = 4, font = 2)
abline(h = 2)
headertext = c("m", "k", "r"); for(i in 1:length(ilab.xpos)){
  text(ilab.xpos[i], 16, headertext[i], adj = c(0.5, 0.5), font = 2, cex = 0.8)
}; rm(headertext); rm(i)
dev.off(); rm(leftmost)

########## end of script ########## 