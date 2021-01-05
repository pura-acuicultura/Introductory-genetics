##### DEVELOPING POPULATIONS IN ALPHASIMR #####

library(AlphaSimR)
library(tidyverse)
library(Matrix)
library(ape)

# vignette("intro", package = "AlphaSimR")

# 1. CREATE FOUNDER HAPLOTYPES #####
#__ 1.1 Founder haplotypes with quickHaplo #####
pop_haplos <- quickHaplo(nInd = 10, # number of individuals
                         nChr = 1, # number of chromosomes
                         segSites = 10, # number of segregation sites
                         genLen = 1, # genetic length of chromosomes
                         ploidy = 2L, # ploidy level of organism
                         inbred = TRUE) # are founders inbred?

#__ 1.2 General haplotype calls #####
str(pop_haplos) # slots for MapPop-class
pop_haplos@genMap[1] # returns the location of the segregation sites
pop_haplos@centromere[1] # returns the location of the centromere
pop_haplos@nInd # number of individuals
pop_haplos@nChr # number of chromosomes
pop_haplos@ploidy # ploidy of individuals
pop_haplos@nLoci[1] # returns vector of loci per chromosome
pop_haplos@geno[1] # what is returned here, a matrix for each individual

# matrix of haplotypes, loci x individual (0/1 on each diploid chromosome)
haplos <- pullSegSiteHaplo(pop = pop_haplos, # Pop-class object
                           haplo = "all", # "all" or use 1 for males and 2 for females
                           chr = 1) # chromosome number, NULL = all
class(haplos)
str(haplos)
haplos

# 2. SET SIMULATION PARAMETERS #####

SP = SimParam$new(pop_haplos) # create a variable for new simulation parameters
str(SP)

#__ 2.1 Additive traits #####
SP$addTraitA(nQtlPerChr = 20, # number of QTL for trait
             mean = 500, # mean genetic value of the trait
             var = 50, # variance of trait
             corA = NULL, # matrix of correlations between additive effects
             gamma = FALSE, # to use a gamma distribution in place of a normal distribution
             shape = 1, # shape parameter for gamma distribution only
             force = FALSE) # keep false until this is understood!

# set environmental variance to create phenotypes #
SP$setVarE(h2 = 0.5, # vector of narrow sense heritability
           H2 = NULL, # vector of broad sense heritability
           varE = NULL) # vector of error variances 

#__ 2.2 SNP chips #####
SP$addSnpChip(nSnpPerChr = 10, # number of SNPs per chromosome
              minSnpFreq = NULL, # minimum allowable frequency for SNP loci
              refPop = NULL) # reference population for calculating SNP frequency

#__ 2.3 Gender ####
# this has seemingly been changed to setSexes
SP$setGender("yes_sys") # assigns gender systematically (male, female, male...)
SP$setGender("yes_rand") # randomly assigns gender
SP$setGender("no") # no gender assignment is the default

#__ 2.4 Recombination ratio #####
SP$setRecombRatio(2) # alters female recombination ration relative to males

#__ 2.5 Parameter tracking #####
SP$setTrackPed(TRUE) # tracks pedigrees
SP$setTrackRec(TRUE) # tracks recombination

# 3. CREATE A FOUNDER POPULATION FROM HAPLOTYPES AND SIM PARAMETERS #####

pop_0 <- newPop(rawPop = pop_haplos, # object of MapPop-class
                mother = NULL, # optional id for mothers
                father = NULL, # optional id for fathers
                origM = NULL, # optional alternative id for mothers
                origF = NULL, # optional alternative id for fathers
                isDH = FALSE, # indicate if double haploids and/or inbred
                simParam = SP) # SimParam objects

# SET PHENOTPYE OF THE POPULATION # [this may not be useful initially]
pop_0 <- setPheno(pop = pop_0, # Pop-class object 
                  varE = 10, # error variances for phenotype
                  reps = 1, # for replicated data
                  fixEff = NULL, # fixed effects to assign to population in GS models
                  p = NULL, # p value for environmental covariate of GE traits
                  onlyPheno = FALSE, # False returns pop-class, true returns matrix of phenotypes
                  simParam = SP) # SimParam object

# 4. OBSERVING POPULATIONS #####
#__ 4.1 General Pop-class calls #####
pop_0@id # vector of individual IDs
pop_0@mother # vector of dam IDs     
pop_0@father # vector of sire IDs
pop_0@gender # vector of individual genders
pop_0@nTraits # number of traits
pop_0@gv # vector of individual genetic values 
pop_0@pheno # vector of individual phenotypes
pop_0@ebv # vector of estimated breeding values (need to be added through SP)
pop_0@gxe[1] # list of GE effects (need to be added through SP)
pop_0@fixEff # vector of fixed effects (1 unless added through SP)
pop_0@reps # vector of replicates

#__ 4.2 Assessment of heritability #####
df <- data.frame(x = pop_0@pheno, y = pop_0@gv)
df %>% ggplot(aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y ~ x)

h2 >- cov(df$x, df$y) / var(df$x) # equal to heritability and slope of graph
h2

#__ 4.3 Observing segregation sites #####
# matrix of haplotypes, loci x individual (0/1 on each diploid chromosome)
pullSegSiteHaplo(pop = pop_0, # Pop-class object
                 haplo = "all", # "all" or use 1 for males and 2 for females
                 chr = 1, # chromosome number, NULL = all
                 simParam = SP) # SimParam object

# matrix of genotypes, loci x individual (NB: 0/0 = 0, 0/1 and 1/0 = 1, 1/1 = 2)
pullSegSiteGeno(pop = pop_0, # Pop-class object
                chr = 1, # chromosome number, NULL = all
                simParam = SP) # SimParam object

# assessment of relatedness
streeb = nj(dist.gene(pullSegSiteGeno(pop = pop_0, chr = 1, simParam = SP)))
plot(streeb)

#__ 4.4 Observing quantitative trait loci (QTL) #####
# dataframe of QTL id, chromosome and physical location in morgans
getQtlMap(trait = 1, # trait number
          gender = "A", # "A" for average, "F" for female, "M" for male
          simParam = SP) # SimParam object

# matrix of QTL haplotypes, loci x individual (0/1 on each diploid chromosome)
pullQtlHaplo(pop = pop_0, # Pop-class object
             trait = 1, # trait number
             haplo = "all", # "all" or use 1 for males and 2 for females
             chr = 1, # chromosome number, NULL = all
             simParam = SP) # SimParam object

# matrix of QTL genotypes, loci x individual (NB: 0/0 = 0, 0/1 and 1/0 = 1, 1/1 = 2)
pullQtlGeno(pop = pop_0, # Pop-class object
            trait = 1, # trait number
            chr = 1, # chromosome number, NULL = all
            simParam = SP) # SimParam object

#__ 4.5 Observing single nucleotide polymorphisms (SNP) #####
# dataframe of SNP id, chromosome and physical location in morgans
getSnpMap(snpChip = 1,
          gender = "A", # A = average, F = female, M = male
          simParam = SP)

# matrix of SNP haplotypes, loci x individual (0/1 on each diploid chromosome)
pullSnpHaplo(pop = pop_0, # Pop-class object
             snpChip = 1, # which chip to use
             haplo = "all", # also 1 = female or 2 = male
             chr = 1, # chromosome number, NULL = all
             simParam = SP) # SimParam object

# matrix of SNP genotypes, loci x individual (NB: 0/0 = 0, 0/1 and 1/0 = 1, 1/1 = 2)
pullSnpGeno(pop = pop_0, # Pop-class object
            snpChip = 1, # which chip to use
            chr = 1, # chromosome number, NULL = all
            simParam = SP) # SimParam object

#__ 4.6 Call function #####

gv(pop_0) # matrix of trait x individual genetic values
meanG(pop_0) # mean genetic value of the population
mean(gv(pop_0)) # verification of meanG() calculation

pheno(pop_0) # matrix of trait x individual phenotypic values
meanP(pop_0) # mean phenotypic value of the population
mean(pheno(pop_0)) # verification of the meanP() calculation

### individual measures
bv(pop = pop_0,
   simParam = SP) # vector of "TRUE" breeding values for traits
gv(pop_0) - mean(gv(pop_0)) # verification of bv() calculation

aa(pop = pop_0,
   simParam = SP) # vector of epistatic deviations
dd(pop = pop_0,
   simParam = SP) # dominance deviations

# population level genic variance (maybe expressed a percentages?)
genicVarA(pop_0,
          simParam = SP) # additive genic variance
genicVarD(pop_0,
          simParam = SP) # dominance genic variance
genicVarAA(pop_0,
           simParam = SP) # epistatic genic variance
genicVarG(pop_0,
          simParam = SP) # sum of A, D, and AA

# population level trait variance
varA(pop_0, simParam = SP) # additive variance
varD(pop_0, simParam = SP) # dominance variance
varAA(pop_0, simParam = SP) # epistatic variance
varG(pop_0) #  total genetic variance

varP(pop_0) #  total phenotypic variance
var(pheno(pop_0)) * (nInd(pop_0) - 1)/nInd(pop_0) # verification of variance calculation
var(pheno(pop_0))

# heritability calculations

var(gv(pop_0)) / var(pheno(pop_0))
varA(pop_0, simParam = SP) / varP(pop_0)

# call all genetic parameters
genParam(pop_0, simParam = SP)

# 5. MODIFYING POPULATIONS #####

#__ 5.1 Assigning traditional breeding values #####
# traditional breeding value (cannot be directly added - needs to be RRsol-class object)

q_gen_bvs <- (pheno(pop_0) - mean(pheno(pop_0))) * h2

#__ 5.2 Assigning BLUP from pedigree #####

#__ 5.3 Assigning BLUP from genomic prediction #####

ans <- RRBLUP(pop = pop_0, # Pop-class object
              traits = 1, # number of traits to model
              use = "pheno", # train model using pheno, gv, ebv, bv, or rand
              snpChip = 1, # which chip to use if more than 1
              useQtl = FALSE, # should QTL genotype be used instead of SNP chip
              maxIter = 1000L, # only used when more than 1 trait is being modeled
              useReps = FALSE, # should reps be used to model heterogeneous error
              simParam = SP) # SimParam object

ans
ans@gv[[1]]@addEff - mean(ans@gv[[1]]@addEff)


pop_0 <- setEBV(pop = pop_0, # Pop-class object
              solution = ans, # RRsol-class object
              value = "bv", # can be either gv, bv, male, or female
              targetPop = NULL, # an optional target population for bv, male, female
              append = FALSE, # replaces existing EBV by default, true adds an additional column
              simParam = SP) # SimParam object

bv(pop_0) # unaffected
ebv(pop_0) # as gv seems to add intercept (actuall different - just tiny), bv adds an unlisted criterion
gv(pop_0) # unaffected

cor(gv(pop_0), ebv(pop_0))
plot(ebv(pop_0), gv(pop_0))

# understanding RRsol-class
ans # two main lists 'gv' and 'bv'

ans@gv[[1]]
ans@gv[[1]]@addEff
ans@gv[[1]]@intercept
ans@gv[[1]]@nLoci
ans@gv[[1]]@lociPerChr
ans@gv[[1]]@lociLoc

ans@bv[[1]]
ans@bv[[1]]@addEff
ans@bv[[1]]@intercept
ans@bv[[1]]@nLoci
ans@bv[[1]]@lociPerChr
ans@bv[[1]]@lociLoc

# 6. SELECTING FROM POPULATIONS #####

# selecting individuals (p. 68)
pop_2 <- selectInd(pop = pop_0,  # Pop-class object
          nInd = 5, # number of individuals to select
          trait = 1, # trait for selection
          use = "pheno", # select on gv, ebv, bv, pheno, or rand
          gender = "B", # B for both; F for female; M for males
          selectTop = TRUE, # True = highest values; False = lowest values
          returnPop = TRUE, # True = results returned as Pop-cass; False = results as index of selected individuals
          candidates = NULL, # optional vector of selection candidates
          simParam = NULL) # SimParam object

# 7. MATING NEW POPULATIONS #####
#__ 7.1 Hybrid crossing (p. 19) #####
# only used when new pop are double haploids or inbred

pop_1 <- hybridCross(females = pop_0f, # female Pop-class
                     males = pop_0m, # male Pop-class
                     crossPlan = "testcross", # full factorial cross; alternatively a two-column matrix can be supplied
                     returnHybridPop = FALSE, # returns a Pop-class; TRUE used for fully inbred crosses
                     simParam = SP) # SimParam object

#__ 7.2 Random crossing (p. 42) #####

# selects parental combinations for all possible combinations (p. 42)
pop_1 <- randCross(pop = pop_0, # Pop-class object
                   nCrosses = 25, # number of crosses to make
                   nProgeny = 1, # number of progeny per cross
                   balance = TRUE, # balance the progeny per parent if using gender
                   parents = NULL, # vector of allowable parents
                   ignoreGender = FALSE, # can gender be ignored?
                   simParam = SP) # simParam object

# randCross2

#__ 7.3 User designed crossing #####

# make a cross plan for user designed crossing
crossPlan <- matrix(c(1, 10),
                    nrow = 1,
                    ncol = 2) 

# from a single population (p. 21)
pop_1 <- makeCross(pop = pop_0, # Pop-class object
                   crossPlan = crossPlan, # two column matrix indicating female/male cross made
                   nProgeny = 10, # number of progeny per cross
                   simParam = SP) # simParam object

# from separate female/male populations (p. 22)
pop_1 <- makeCross(pop = pop_0f, # female Pop-class object
                   pop = pop_0m, # male Pop-class object
                   crossPlan = crossPlan, # two column matrix indicating female/male cross made
                   nProgeny = 10, # number of progeny per cross
                   simParam = SP) # simParam object

#__ 7.4 To be completed #####
str(pop_2)
pedigree <- data.frame(ind = pop_2@id,
                       sex = pop_2@gender,
                       sire = pop_2@father,
                       dam = pop_2@mother)
pedigree

# selectCross
# pedigreeCross
