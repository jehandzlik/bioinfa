## 10-May-2022
## Normalization

## Redistribute percentages of expressed gene 
redistribute <- function(percentages){
  
  return(100 * percentages/ sum(percentages))
  
}

##########################################################
## A) Same genes, different library size
## DIFFERENCE IN LIBRARY SIZE
##########################################################

# control genes
control = c(10,300,500,20,100,8)

# treatment genes, remain the same, different library size
treatment_diffLibSize = c(10,300,500,20,100,8) * 20

## Solution: RPM
## The different library size can be addressed with RPM
10^6 * control/sum(control)

10^6 * treatment_diffLibSize/sum(treatment_diffLibSize)

##########################################################
## B) Same library size, differentially expressed genes
## DIFFERENCE IN LIBRARY COMPOSITION
##########################################################

## nr of transcripts
N = 50000

## Percentage of total N assigned to genes
gene1 = 5
gene2 = 20
gene3 = 30
gene4 = 9
gene5 = 1
gene6 = 35

## Remove the last gene
control_prob = c(gene1, gene2, gene3, gene4, gene5,0)

control = c(gene1, gene2, gene3, gene4, gene5, gene6)*N/100

treatment = redistribute(control_prob)*N/100

## RPM doesn't work anymore
10^6 * control/sum(control)
10^6 * treatment/sum(treatment)

df <- data.frame(sample1 = control, sample2 = treatment, sample3 = treatment*20)
df <- data.frame(sample1 = control, sample2 = treatment)

DEseq2_norm <- function(df){
  
  # step 1: log each gene
  df_log <- log(df)
  
  # step 2: average each gene
  log_aver <- data.frame(Average = apply(df_log, 1, mean))
  
  # step 3: filter out genes with infinity - in theory, this helps
  # focus the scaling factors on the house keeping genes -
  # genes transcribed at similar levels regardless of cell type
  df_log <- df_log[-which(log_aver$Average == -Inf),]
  log_aver <- log_aver[-which(log_aver$Average == -Inf),]
  
  # step 4: substract the average log value from the log(counts)
  # substracting logs also means division so what is checked here is 
  # the ratio of the reads in each sample to the average across samples
  # log(reads for gene X / average for gene X)
  df_log_minus_aver <- df_log - log_aver
  
  # step 5: calculate the median of the ratios for each sample 
  log_medians <- apply(df_log_minus_aver, 2, median)
  
  # step 6: convert the medians to "normal numbers" to get the scaling 
  # factors
  scaling_factors <- exp(1)^log_medians
  
  return(scaling_factors)
}

scaling_factors <- DEseq2_norm(df)

## This is exactly what we wanted!
control/scaling_factors[1]
treatment/scaling_factors[2]

sample1 = control
sample2 = treatment
sample3 = treatment*20

sample1/scaling_factors[1]
sample2/scaling_factors[2]
sample3/scaling_factors[3]

control_diffLibSize = c(10,300,500,20,100,8) * 20

# treatment genes
# two genes are de, one expressed more and one expressed less
treatment = c(10,300*20,500,20,100*0.1,8)

## B) Different library size, different genes
treatment_diffLibSize = c(10,300*20,500,20,100*0.1,8) * 20

##########################################################

## nr of transcripts
N = 50000

## Percentage of total N assigned to genes
gene1 = 5
gene2 = 20
gene3 = 30
gene4 = 9
gene5 = 1
gene6 = 35

## Remove the last gene
control_prob = c(gene1, gene2, gene3, gene4, gene5)

control = c(gene1, gene2, gene3, gene4, gene5, gene6)*N/100

treatment = redistribute(control_prob)*N/100


gene_lengths <- c(10^3, 2*10^2, 4*10^3, 1000, 1000, 1000)
####################################################
## RPM
####################################################

A <- 10^6 * control/sum(control)

B <- 10^6 * control_diffLibSize/sum(control_diffLibSize)

C <- 10^6 * treatment/sum(treatment)

D <- 10^6 * treatment_diffLibSize/sum(treatment_diffLibSize)

## A == B and that happens only when the cells are of similar type with
## no DE genes. Here the difference in library size can be accounted by
## the normalization by the library length.

####################################################
## RPKM
####################################################

## Here the genes are divided by their transcript length in kilobases. 
A <- (10^6 * control/sum(control)) / gene_lengths
## etc....


