##### SIMULATING A TILAPIA LIKE GENOME #####

library(AlphaSimR)
library(tidyverse)
library(chromPlot)

##### 1. DEVELOP PRELIMINARY MAP DATA.FRAME #####

# 1A. Load tilapia genome data taken from NCBI genome
a <-read.csv(file = 'simulated tilapia genome.txt', header = TRUE)
a <- a[,c(1, 6, 12, 14)] # eliminate a few unnecessary columns
a <- a[-c(23,24),] # eliminate a few unnecessary rows

# 2A. Pull some data needed to guide haploid development
Loci <- as.numeric(gsub(",", "", a$Gene))
Morgans <- a$Length..cM. / 100

# 3A. Create the founder population haploids
pop_haps <- quickHaplo(nInd = 10,              # number of individuals
                       nChr = length(a$Chrom), # number of chromosomes
                       segSites = Loci,        # number of segregation sites
                       genLen = Morgans,       # genetic length of chromosomes
                       ploidy = 2L,            # ploidy level of organism
                       inbred = FALSE)         # are founders inbred?

SP = SimParam$new(pop_haps)
SP$addSnpChip(nSnpPerChr = 3,   # number of SNPs per chromosome
              minSnpFreq = NULL, # minimum allowable frequency for SNP loci
              refPop = NULL)
pop_0 <- newPop(pop_haps, 
                simParam = SP) 

Chrom <- a$Chrom
Start <- rep(0, nrow(a))
End <- a$Size..Mb.* 1000000
Name <- rep('contig', nrow(a))

tg_map <- data.frame(Chrom, Start, End, Name)
chromPlot(gaps = tg_map)

##### 2. ADD TO THE MAP FILE #####

# 2A. Add centromeres
a <- cbind(a, pop_haps@centromere)
colnames(a) <- c("Chrom", "Mb", "Genes", "cM", "Centro")
a <- a %>% mutate(cen_Mb = Centro * (Mb * 1000000) / (cM / 100))

tg_centro <- data.frame(Chrom,
                        'Start' = a$cen_Mb,
                        'End' = a$cen_Mb,
                        'Name' = "centromere")
tg_map <- rbind(tg_map, tg_centro)
chromPlot(gaps = tg_map)

# 2B. Add SNP to map
b <- a[, c(1, 2, 4)]

c <- getSnpMap(snpChip = 1,
               simParam = SP)

colnames(c) <- c("snp", "Chrom", "Pos")

d <- merge(c, b, by = 'Chrom')
d <- d %>% mutate(snp_bp = Pos * (Mb * 1000000) / (cM / 100))
head(d)

Chrom <- d$Chrom
Start <- d$snp_bp
End <- d$snp_bp + 1
ID <- paste("SNP", d$snp, sep = "_")
Colors <- rep("darkgreen", length(Chrom))

AIMS <- data.frame(Chrom, Start, End, ID, Colors)
head(AIMS)

chromPlot(gaps = tg_map,
          stat = AIMS,
          statCol = "Value",
          chr = c(1:8))

