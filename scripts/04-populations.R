##### NEUTRAL MARKERS AND POPULATION GENETICS #####

library(AlphaSimR)
library(tidyverse)
library(Matrix)
library(ape)
library(ggtree)
library(phangorn)

##### 1. CREATE TEST POPULATION #####

pop_haps <- quickHaplo(nInd = 10,      # number of individuals
                         nChr = 1,       # number of chromosomes
                         segSites = 10,  # number of segregation sites
                         genLen = 1,     # genetic length of chromosomes
                         ploidy = 2L,    # ploidy level of organism
                         inbred = FALSE) # are founders inbred?

SP = SimParam$new(pop_haps) # create a variable for new simulation parameters

SP$addSnpChip(nSnpPerChr = 5,    # number of SNPs per chromosome
              minSnpFreq = NULL, # minimum allowable frequency for SNP loci
              refPop = NULL)     # reference population for calculating SNP frequency

pop1 <- newPop(pop_haps, 
               simParam = SP)

##### 2. GENE and SNP MARKER FREQUENCIES #####

# tabulate actual gene frequencies
a1 <- pullSegSiteHaplo(pop = pop1,     
                 haplo = "all",
                 chr = 1,       
                 simParam = SP) # Pull loci haplotypes
a2 <- matrix(a1,
             nrow = nrow(a1), 
             ncol = ncol(a1), 
             byrow = FALSE) # format loci haplotypes as an unlabeled matrix
a3 <- colSums(a2 == 0) / nrow(a1) # calculate frequency of '1' allele at all loci
a4 <- data.frame(origin = "loci",
                 pos = round((unlist(pop_haps@genMap)), 2),
                 p = a3,
                 q = 1 - a3) # set up a data frame of loci position and '1' allele frequency
head(a4)

# tabulate SNP marker frequencies
b1 <- pullSnpHaplo(pop = pop1,
             snpChip = 1,
             haplo = 'all',
             chr = NULL,
             simParam = SP) # pull SNP haplotypes
b2 <- matrix(b1,
             nrow = nrow(a1), 
             ncol = ncol(b1), 
             byrow = FALSE) # format SNP haplotypes as a matrix
b3 <- colSums(b2 == 0) / nrow(b1) # calculate frequency of '0' allele at all SNP
b4 <- getSnpMap(snpChip = 1,
                simParam = SP) # pull SNP map
b5 <- data.frame(origin = "snp",
                 pos = b4[3],
                 p = b3,
                 q = 1 - b3) # set up a data frame of SNP position and '1' allele frequency
head(b5)

# combine gene and snp frequencies, and plot
a5 <- rbind(a4, b4) # combine loci and SNP data frames
a5 %>% ggplot(aes(x = pos, y = gf, color = origin)) +
  geom_point() +
  geom_line() # graph to identify SNP locations

mean(a2) # compare the loci and SNP mean frequencies (NB: is this a measure of heterozygocity?)
mean(b2)

##### 3. GENOTYPE FREQUENCIES, HETEROZYGOCITY, HWE, AND F STATS #####
c1 <- pullSegSiteGeno(pop = pop1,      # Pop-class object
                      chr = 1,        # chromosome number, NULL = all
                      simParam = SP)

c1 <- matrix(c1,
             nrow = pop1@nInd, 
             ncol = pop1@nLoci, 
             byrow = FALSE)

obs_hetZ <- colSums(c1 == 1) / nrow(c1)
p_x <- (2 * colSums(c1 == 0) + colSums(c1 == 1)) / (2 * nrow(c1))
q_x <- 1 - p_x
exp_hetZ <- (2 * p_x * q_x)

hetZ <- data.frame(loci = 1:pop1@nLoci, 
                   f_st = 1 - obs_hetZ / exp_hetZ)

hetZ %>% ggplot(aes(x = loci, y = f_st)) +
  geom_point(color = 'blue') +
  geom_hline(yintercept = mean(hetZ$f_st),
             linetype = 'dashed',
             color = 'black',
             size = 1.5)
# p_x
# q_x
mean(obs_hetZ)
mean(exp_hetz)

# histogram of HWE 
e1 <- data.frame(origin = "observed",
                 loci = 1:pop1@nLoci,
                 hetZ = obs_hetZ)
e2 <- data.frame(origin = "expected",
                 loci = 1:pop1@nLoci,
                 hetZ = exp_hetZ)
e3 <- rbind(e1, e2)

e3 %>% ggplot(aes(x = hetZ, ..count.., fill = origin)) +
  geom_density(adjust = 2, alpha = 0.5) +
  geom_vline(xintercept = mean(e2$hetZ),
             linetype = 'dashed',
             color = 'red',
             size = 1) +
  geom_vline(xintercept = mean(e1$hetZ),
             linetype = 'dashed',
             color = 'dark green',
             size = 1)

e3 %>% ggplot(aes(x = hetZ, fill = origin)) +
  geom_histogram(binwidth = 0.1, position="identity", alpha = 0.5) +
  geom_vline(xintercept = mean(e2$hetZ),
             linetype = 'dashed',
             color = 'red',
             size = 1) +
  geom_vline(xintercept = mean(e1$hetZ),
             linetype = 'dashed',
             color = 'dark green',
             size = 1)


# data = subset(e3, origin == 'observed'), fill = "green", alpha = 0.2, 
# data = subset(e3, origin == 'expected'), fill = "blue", alpha = 0.2)

##### 4. INFERENCE FROM SNP MARKERS #####
c1 <- pullSegSiteGeno(pop = pop1,      # Pop-class object
                      chr = 1,        # chromosome number, NULL = all
                      simParam = SP)
streec1 = nj(dist.gene(c1))
plot(streec1)

d1 <- pullSnpGeno(pop = pop1,      # Pop-class object
                  chr = 1,        # chromosome number, NULL = all
                  simParam = SP)
streed1 = nj(dist.gene(d1))
plot(streed1)

# tree comparison tools
treedist(streec1, streed1, check.labels = TRUE)
RF.dist(streec1, streed1)
sprdist(streec1, streed1)

# observing both plots together (APE)
d3 <- comparePhylo(streec1, streed1,
                   plot = TRUE,
                   force.rooted = TRUE,
                   use.edge.length = FALSE)

d4 <- matrix(c(1:10, 1:10), nrow = 10, ncol = 2)
cophyloplot(streec1, streed1, d4)

# observing both plots together (ggtree)
x <- streec1
y <- streed1
p1 <- ggtree(x)
p2 <- ggtree(y)

d1 <- p1$data
d2 <- p2$data

## reverse x-axis and set offset to make the tree in the right hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1

pp <- p1 + 
  geom_tiplab() + 
  geom_tree(data=d2) + 
  geom_tiplab(data = d2, hjust=1)

dd <- bind_rows(d1, d2) %>% 
  filter(!is.na(label))

pp + geom_line(aes(x, y, group=label), data=dd, color='grey')

