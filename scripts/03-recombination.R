##### INTRODUCTION TO RECOMBINATION #####

library(AlphaSimR)
library(tidyverse)
library(Matrix)
library(ape)
library(data.table)
library(chromPlot)
library(gridExtra)

##### 1. CREATE PARENTS #####

founderPop <- quickHaplo(nInd = 2,       # number of individuals
                         nChr = 1,       # number of chromosomes
                         segSites = 1000,  # number of segregation sites
                         genLen = 0.1,     # genetic length of chromosomes
                         ploidy = 2L,    # ploidy level of organism
                         inbred = FALSE) # are founders inbred?

SP <- SimParam$new(founderPop) # create a variable for new simulation parameters

SP$setTrackRec(TRUE)

pop <- newPop(founderPop, 
              simParam = SP)

##### 2. TABULATE GENE FREQUENCIES IN PARENTS #####

gfreq <- pullSegSiteHaplo(pop = founderPop, 
                          haplo = "all",                                     chr = 1,                   # chromosome number, NULL = all
                          simParam = SP)            
gfreq <- data.frame(gfreq)
gfreq <- setDT(gfreq, keep.rownames = TRUE)[]

colnames(gfreq) <- c("rn", 1:10)
gfreq <- gather(data = gfreq,
                key = "loci",
                value = "g_type", -rn)
gfreq <- gfreq %>% separate(rn, into = c("ind", "allele"), sep = "_")
gfreq_1 <- gfreq %>% 
  group_by(loci) %>%
  count(g_type)
gfreq_1 <- as.data.frame(gfreq_1)
gfreq_1 <- spread(data = gfreq_1,
                  key = "g_type",
                  value = "n",
                  replace_na(0))
colnames(gfreq_1) <- c("loci", "p", "q")
gfreq_1 <- gfreq_1 %>% 
  mutate(p_o = p/(p+q)) %>%
  mutate(q_o = q/(p+q))
gfreq_1

##### 3. CROSS PARENTS #####
crossPlan <- matrix(c(1, 2), nrow = 1, ncol = 2) # develop a cross plan
pop2 <- makeCross(pop = pop,
                  crossPlan = crossPlan,
                  nProgeny = 100,
                  simParam = SP)

##### 4. TABULATE GENE FREQUENCIES IN OFFSPRING #####
gfreq <- pullSegSiteHaplo(pop = pop2, 
                          haplo = "all",                                     chr = 1,                   # chromosome number, NULL = all
                          simParam = SP)            
gfreq <- data.frame(gfreq)
gfreq <- setDT(gfreq, keep.rownames = TRUE)[]

colnames(gfreq) <- c("rn", 1:10)
gfreq <- gather(data = gfreq,
                key = "loci",
                value = "g_type", -rn)
gfreq <- gfreq %>% separate(rn, into = c("ind", "allele"), sep = "_")
gfreq_2 <- gfreq %>% 
  group_by(loci) %>%
  count(g_type)
gfreq_2 <- as.data.frame(gfreq_2)
gfreq_2 <- spread(data = gfreq_2,
                  key = "g_type",
                  value = "n",
                  replace_na(0))
colnames(gfreq_2) <- c("loci", "p", "q")
gfreq_2 <- gfreq_2 %>% 
  mutate(p_o = p/(p+q)) %>%
  mutate(q_o = q/(p+q))
gfreq_2

##### 5. ASSESS EFFECTS OF RECOMBINATION IN OFFSPRING #####

pop <- c(rep("parents", 10), rep("offspring", 10))
pop

a <- gfreq_1[,-c(2,3)]
b <- gfreq_2[,-c(2,3)]
c <- merge(a, b, by = "loci")
d <- c %>% mutate(p_dif = p_o.x - p_o.y) %>%
           mutate(q_dif = q_o.x - q_o.y)
d <- d[, -c(2, 3, 4, 5)]
d$loci <- as.numeric(d$loci)
d <- gather(data = d,
            key = "allele",
            value = "freq",
            -loci)
d %>% ggplot(aes(x = loci, y =freq, color = allele)) +
  geom_point()
d
gfreq_1
gfreq_2

##### 6. VERSION 1 ATTEMPT TO MAP RECOMBINATION #####
SP$pedigree
SP$isTrackPed
SP$recHist
str(SP$recHist) # not exactly sure what is being returned here...

# access parent haplotypes 
f0 <- pullSegSiteHaplo(pop = founderPop, 
                       haplo = "all",                                     chr = 1,                   # chromosome number, NULL = all
                       simParam = SP)
f0

# access the offspring haplotypes
f1 <- pullSegSiteHaplo(pop = pop2, 
                       haplo = "all",                                     chr = 1,                   # chromosome number, NULL = all
                       simParam = SP)
f1

# make vectors of paternal haplotypes 1 and 2
male_gam_1 <- as.vector(f0[1,])
male_gam_1 <- matrix(male_gam_1, nrow = 10, ncol = 10, byrow = TRUE)
male_gam_1

male_gam_2 <- as.vector(f0[2,])
male_gam_2 <- matrix(male_gam_2, nrow = 10, ncol = 10, byrow = TRUE)
male_gam_2

# make vectors of offspring paternal haplotypes
pad_hap <- f1[c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19),]
pad_hap <- matrix(pad_hap, nrow = 10, ncol = 10)
pad_hap

# determine the paternal haplotype of the offspring
a <- rowSums((male_gam_1 - pad_hap)^2)
b <- rowSums((male_gam_2 - pad_hap)^2)

c <- b < a # TRUE = paternal haplotype 2; FALSE = paternal haplotype 1
c <- as.integer(as.logical(c))
c <- c + 1
c

# create a hybrid matrix of paternal haplotypes

hyb_haps <- matrix(rep(as.vector(f0[c,])), nrow = 10, ncol = 10)

# determine the number of recombination events at each loci
colSums((hyb_haps - pad_hap)^2)

hyb_haps - pad_hap

##### 7. VERSION 2 ATTEMPT TO MAP RECOMBINATION #####

### step 1. make a parent haplotypes data frame ### 

# set up data.frame of parent qualifiers
HAP <- 1:(founderPop@nChr * founderPop@nInd * founderPop@ploidy)
CHR <- rep(1:founderPop@nChr, founderPop@nInd * founderPop@ploidy)
ID <- rep(1:founderPop@nInd, each = founderPop@ploidy)
ALE <-rep(1:founderPop@ploidy, founderPop@nInd)
PAR_Q <- data.frame(HAP, CHR, ID, ALE)

# set up matrix of parent haplotypes
par_H <- pullSegSiteHaplo(pop = founderPop, 
                          haplo = "all",                                     chr = 1,                   # chromosome number, NULL = all
                          simParam = SP)
par_H <- matrix(par_H,
                nrow = founderPop@nInd * founderPop@ploidy, 
                ncol = founderPop@nLoci, 
                byrow = FALSE)
colnames(par_H) <- paste("L",
                          1:founderPop@nLoci,
                          sep = ".")
PAR_H <- as.data.frame(cbind(PAR_Q, par_H))
par_H
PAR_H

### step 2. make an offsping haplotypes data frame ###

# set up data.frame of offspring qualifiers
HAP <- 1:(pop2@nChr * pop2@nInd * pop2@ploidy)
CHR <- rep(1:pop2@nChr, pop2@nInd * pop2@ploidy)
ID <- rep(1:pop2@nInd, each = pop2@ploidy)
ALE <-rep(1:pop2@ploidy, pop2@nInd)
OFF_Q <- data.frame(HAP, CHR, ID, ALE)
OFF_Q

# set up matrix of offspring haplotypes
off_H <- pullSegSiteHaplo(pop = pop2, 
                          haplo = "all",                                     chr = 1,                   # chromosome number, NULL = all
                          simParam = SP)
off_H <- matrix(off_H,
                nrow = pop2@nInd * pop2@ploidy, 
                ncol = pop2@nLoci, 
                byrow = FALSE)
colnames(off_H) <- paste("L",
                         1:pop2@nLoci,
                         sep = ".")
OFF_H <- as.data.frame(cbind(OFF_Q, off_H))
off_H
OFF_H

### step 3. assign parental haplotypes to offspring haplotypes ###

par_SEL <- as.data.frame(matrix(0, nrow = pop2@nInd * pop2@ploidy,
                                ncol = founderPop@nInd * founderPop@ploidy))
colnames(par_SEL) <- 1:(founderPop@nInd * founderPop@ploidy)

par_SEL

for(i in 1:(founderPop@nInd * founderPop@ploidy)){
  par_SEL[,i] <- rowSums((off_H - matrix(rep(par_H[i,], 
                                    pop2@nInd * pop2@ploidy),
                                nrow = pop2@nInd * pop2@ploidy,
                                ncol = pop2@nLoci,
                                byrow = TRUE))^2) 
}

PAR_SEL <- as.numeric(colnames(par_SEL)[apply(par_SEL,
                     1,
                     which.min)])



REC <- as.data.frame(matrix(0, nrow = pop2@nInd * pop2@ploidy,
                            ncol = pop2@nLoci))

for(ii in 1:(pop2@nInd * pop2@ploidy)){
  REC[ii,] <- off_H[ii,] - par_H[PAR_SEL[ii],]
}

colnames(REC) <- 1:pop2@nLoci
REC <- cbind(OFF_Q, REC)

rec <- gather(data = REC,
       key = "loci",
       value = "rec", -HAP, -CHR, -ID, -ALE) 

rec_map <- rec %>% filter(ALE == 1) %>%
                    group_by(loci) %>%
                    summarise(Rec = sum(rec^2))

rec_map <- as.data.frame(rec_map)
rec_map$loci <- as.numeric(rec_map$loci)

plot_male <- rec_map %>% ggplot(aes(x = loci, y = Rec)) +
              geom_line()

rec_map <- rec %>% filter(ALE == 2) %>%
  group_by(loci) %>%
  summarise(Rec = sum(rec^2))

rec_map <- as.data.frame(rec_map)
rec_map$loci <- as.numeric(rec_map$loci)

plot_female <- rec_map %>% ggplot(aes(x = loci, y = Rec)) +
  geom_line()

grid.arrange(plot_male, plot_female, nrow = 2)
