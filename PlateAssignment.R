# This program balances samples on two epigenetics arrays, so that batch effects are minimised 
# not just between the plate but beteween the row as well.

# This program runs as a sequel to case control matching
# For input it requires the output file of the matched case-controls
# which comes from CaseControlMatching.R


HOME_DIRECTORY <- "/Users/timothystone/Desktop/"
matchedControl <- read.csv(paste0(HOME_DIRECTORY, "Final_Set_of_Matched_Controls.csv"))
if (nrow(matchedControl) != 192) { stop ("You do not have 192 matched controls\n")}

# Relative weightings between the covariates

GENDER_WEIGHT <- 4
CANCERGENDER_WEIGHT <- 4
DIAGNOSIS_WEIGHT <- 5
AGE_WEIGHT <- 3
CIG_WEIGHT <- 3
PPI_WEIGHT <- 1
HEARTBURN_WEIGHT <- 1

cat("Starting")

# Initialise variables

count <- 0
bestPValue <- 0

# This program loops forever
# But the text output should indicate that when no new improvements occur after a long search interval, you can probably stop
# We ran our actual final matching over a long weekend

while (1) {
  
  count <- count + 1
  
  # matchCopy will be split into two plates, and tested
  matchCopy <- matchedControl
  
  # Write a messsage after every 100 counts
  if ((count/100) == round(count/100,0)) {
    cat("Attempt number ", count, "\n")
  }
  
  # Create seperate DataFrames for each diagnosis
  cancers <- filter(matchedControl, FinalDiagnosis == "Cancer")
  im <- filter(matchedControl, FinalDiagnosis == "IM")
  hg <- filter(matchedControl, FinalDiagnosis == "HG")
  ndbe <- filter(matchedControl, FinalDiagnosis == "NDBE")
  normal <- filter(matchedControl, FinalDiagnosis == "Normal")
  hv <- filter(matchedControl, FinalDiagnosis == "HV")
  
  # This is the "jittering" section
  # There are two jitters involved, the cancer jitter and the control jitter
  # The cancer samples have a more limited jitter than the various classes of control
  # Consequently we have a two-stage jitter
  
  # First we compute the 50% size of the samples, samples will jitter around that value
  samplesize <- c(nrow(cancers)/2, nrow(im) / 2, nrow(hg) / 2, nrow(ndbe) / 2, nrow(normal) / 2, nrow(hv) / 2)
  
  # This section jitters the cancers, swapping either an adenocarcinoma for an additional intramucosal
  # or vice versa
  cancerChange <- sample(c(0.5, -0.5), 2, replace=F)
  samplesize[1:2] <- samplesize[1:2] - cancerChange
  
  # Now to jitter the controls
  # 2 Randomly selected of the 4 classes of control will be jittered
  # The sampling favours jittering of the people with no oesophageal abnormality: prob = c(1, 2, 4, 4) 
  # Healthy volunteers and "normals" are 4 times more likely to be jittered than the high-grade dysplasias
  # This is because there are more of them, they have a larger permutation space
  # This sampling is repeated 10 times, subjects will therefore "random walk" their jittering
  
  for (i in 1:10) {
    
    jitterChange <- sample(1:4, 2, prob = c(1, 2, 4, 4), replace = F)
    samplesize[2 + jitterChange[1]] <- samplesize[2 + jitterChange[1]] + 0.5
    samplesize[2 + jitterChange[2]] <- samplesize[2 + jitterChange[2]] - 0.5
  }
  
  # This section implements the jittering of the controls, and also removes any non-integers in Samplesize
  adjust <- which(samplesize != round(samplesize,0))
  adjustChange <- sample(c(0.5, -0.5), 2, replace = F)
  samplesize[adjust] <- samplesize[adjust] + adjustChange
  
  # This is a precautionary stop, in case the jittering somehow failed
  # The sum of all the individual group samples should equal 96
  if (sum(samplesize) != 96 ) { stop( "Sample Size error\n")}
  
  # This assembles the final candidate plates
  
  # We have the sample sizes of the groups, we randomly selected from the groups, 
  # These are 96 samples, half the total, they will be assigned to Plate number 1
  Csampl <- sample_n(cancers, samplesize[1]) %>% select(Subject.number) %>% unlist() %>% as.character()
  IMsampl <- sample_n(im, samplesize[2]) %>% select(Subject.number) %>% unlist() %>% as.character()
  HGsampl <- sample_n(hg, samplesize[3]) %>% select(Subject.number) %>% unlist() %>% as.character()
  NDBEsampl <- sample_n(ndbe, samplesize[4]) %>% select(Subject.number) %>% unlist() %>% as.character()
  normsampl <- sample_n(normal, samplesize[5]) %>% select(Subject.number) %>% unlist() %>% as.character()
  hvsampl <- sample_n(hv, samplesize[6]) %>% select(Subject.number) %>% unlist() %>% as.character()
  
  # Make the default plate assignment = 2
  matchCopy[,'Plate'] <- 2
  
  # Change the plate assignmeent of the samples in samplesize to plate = 1
  matchCopy <- mutate(matchCopy,  Plate = replace(Plate, Subject.number %in% Csampl, 1))
  matchCopy <- mutate(matchCopy,  Plate = replace(Plate, Subject.number %in% IMsampl, 1))
  matchCopy <- mutate(matchCopy,  Plate = replace(Plate, Subject.number %in% HGsampl, 1))
  matchCopy <- mutate(matchCopy,  Plate = replace(Plate, Subject.number %in% NDBEsampl, 1))
  matchCopy <- mutate(matchCopy,  Plate = replace(Plate, Subject.number %in% normsampl, 1))
  matchCopy <- mutate(matchCopy,  Plate = replace(Plate, Subject.number %in% hvsampl, 1))
  
  ## THIS IS THE COMPOSITE P-VALUE SECTION
  # The aim is to create a composite, weighted p-value that is maximised
  # The larger the composite p-value, the greater the homogeneity between the plates
  # This is the complete reverse of what a standard statistical test normally does, where a low value is normally desirable
  # In these case we require as high a p-value as possible, so that the base populations are most homogeneous
  
  # Cancer diagnosis requires a Chi-Square (categorical)
  # This requires an additional filter line to remove the non-cancers
  CancerGenderP <- filter(matchCopy, FinalDiagnosis == "Cancer") %>% 
    summarise(pval = chisq.test(Gender, Plate$p.value))
  if (CancerGenderP < 0.05) { 
    #    cat("Gender fail\n")
    next()
  }
  
  # Gender, use a chi-square test
  genderP <- matchCopy %>% 
    summarise(pval = chisq.test(Gender, Plate)$p.value)
  if (genderP < 0.05) { 
    #    cat("Gender fail\n")
    next()
  }
  
  # Diagnosis (which includes the sub-categories) is a chi-square
  DiagnosisP <- matchCopy %>% 
    summarise(pval = chisq.test(FinalDiagnosis, Plate)$p.value)
  
  if (DiagnosisP < 0.05) { 
    #    cat("Diagnosis fail\n")
    next()
  }
  
  
  # Age is continuous, t-test
  ageP <- t.test(age ~ Plate, data = matchCopy)$p.value
  if (ageT < 0.05) {
    #cat ("Age T fail\n")
    next()
  }
  
  # bmi is continuous, t-test
  bmiP <- t.test(bmi ~ Plate, data = matchCopy)$p.value
  if (bmiT < 0.05) {
    #cat("BMI t fail\n")
    next()
  }
  
  # Cigarettes are continuous, t-test
  cigP <- t.test(CigsSmoked ~ Plate, data = matchCopy)$p.value
  if (cigT < 0.05) {
    #cat("CIG t fail\n")
    next()
  }
  
  # PPI and heartburn have been calculated in a manner analagous to the "pack-years" of smoking
  # PPI consumption is a continuous variable, is a t-test
  PPIP <- t.test(PPI ~ Plate, data = matchCopy)$p.value
  if (PPIT < 0.05) {
    #cat("PPI t fail\n")
    next()
  }
  
  # heartburn is continuous, t-test
  hbP <- t.test(hb ~ Plate, data = matchCopy)$p.value
  if (hbT < 0.05) {
    #cat("Heartburn t fail\n")
    next()
  }
  
  # Generate the composite p-score
  candidateP <- genderP * GENDER_WEIGHT + 
    CancerGenderP * CANCERGENDER_WEIGHT + 
    DiagnosisP * DIAGNOSIS_WEIGHT + 
    ageP * AGE_WEIGHT + 
    bmiP * BMI_WEIGHT + 
    cigP * CIG_WEIGHT + 
    PPIP * PPI_WEIGHT + 
    hbP * HEARTBURN_WEIGHT
  
  # Did we improve on the last score?
  if (candidagteP > bestPValue) {
    finalPlate <- matchCopy
    cat("We have a better match!!!\n")
    # Write the file with the improved match and update totalP
    write.csv(finalPlate, paste0(HOME_DIRECTORY, "BestMatch.csv"), row.names=F)
    bestPValue <- candidateP
  }
}