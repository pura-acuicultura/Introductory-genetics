##### DEVELOPING CHROMOSOME MAPS IN SUPPORT OF ALPHASIMR #####

library(AlphaSimR)
library(tidyverse)
library(chromPlot)
library(data.table)

##### MAKING FOUNDER MAP-POPCLASS #####
founderPop <- quickHaplo(nInd = 10,      # number of individuals
                         nChr = 1,       # number of chromosomes
                         segSites = 10,  # number of segregation sites
                         genLen = 10, # genetic length of chromosomes
                         ploidy = 2L,    # ploidy level of organism
                         inbred = FALSE) # are founders inbred?
str(founderPop)

##### CHANGING LOCATION OF CENTROMERE #####
SP = SimParam$new(founderPop)
SP$switchGenMap(genMap = founderPop@genMap[1],
                centromere = 50)
str(SP)
SP$centromere

##### UNDERSTANDING THE MAP-POPCLASS #####

founderPop@genMap[1] 
founderPop@centromere[1] 
founderPop@nLoci
founderPop@nChr

##### PREPARING MAP LOCATIONS #####
loci <- c(1:founderPop@nLoci) # vector of loci IDs
seg_st <- unlist(founderPop@genMap[1]) # vector of segregation site locations
loci_ID <- sprintf("LOCI_%s", 1:founderPop@nLoci)
phys_map <- data.frame(loci_ID, seg_st, loci)
phys_map

##### GENE FREQUENCY OF POPULATION #####
gfreq <- pullSegSiteHaplo(pop = founderPop, # Pop-class object
                 haplo = "all",             # "all" or use 1 for males and 2 for females
                 chr = 1,                   # chromosome number, NULL = all
                 simParam = SP)            # SimParam object
gfreq <- data.frame(gfreq)
gfreq <- setDT(gfreq, keep.rownames = TRUE)[]
colnames(gfreq) <- c("ind", 1:10)
gfreq <- gather(data = gfreq,
                key = "loci",
                value = "g_type", -ind)
gfreq <- gfreq %>% separate(ind, into = c("ind", "allele"), sep = "_")

gfreq_1 <- gfreq %>% 
           group_by(loci) %>%
           count(g_type)
gfreq_1 <- as.data.frame(gfreq_1)
gfreq_1 <- spread(data = gfreq_1,
              key = "g_type",
              value = "n")
colnames(gfreq_1) <- c("loci", "p", "q")
gfreq_1 <- gfreq_1 %>% 
            mutate(p_o = p/(p+q)) %>%
            mutate(q_o = q/(p+q)) %>% 
            mutate(p2_x = p_o^2) %>%
            mutate(pq_x = 2 * p_o * q_o) %>%
            mutate(q2_x = q_o^2) 
gfreq_1



##### GENOTYPIC DISTRIBUTION OF POPULATIONS #####
gtype <- pullSegSiteGeno(pop = founderPop,    # Pop-class object
                     chr = 1,              # chromosome number, NULL = all
                     simParam = SP)        # SimParam object
gtype <- data.frame(gtype)
gtype <- setDT(gtype, keep.rownames = TRUE)[]
colnames(gtype) <- c("ind", 1:10)
gtype<- gather(data = gtype,
                key = "loci",
                value = "g_type",
                -ind)
gtype
gtype_1 <- gtype %>% group_by(loci) %>%
              count(g_type)
gtype_1 <- as.data.frame(gtype_1)
gtype_1 <- spread(data = gtype_1,
              key = "g_type",
              value = "n",
              replace_na(0))

colnames(gtype_1) <- c("loci", "p2", "pq", "q2")
gtype_1 <- gtype_1 %>%
           mutate(p2_o = p2/10) %>%
           mutate(pq_o = pq/10) %>%
           mutate(q2_o = q2/10) %>%
           mutate(p_x = (2 * p2 + pq)/20) %>%
           mutate(q_x = (2 * q2 + pq)/20)
gtype_1

##### HARDY-WEINBEG EQUILIBRIUM OF POPULATION #####
a <- merge(phys_map, gfreq_1, by = "loci")
b <- merge(a, gtype_1, by = "loci")
c <- b[,c(1, 8, 9, 10, 14, 15, 16)]
c %>% mutate(p2_d = p2_o - p2_x) %>%
      mutate(pq_d = pq_o - pq_x) %>%
      mutate(q2_d = q2_o - q2_x)

b
d <- b[, c(1, 6, 7, 17, 18)]
d





##### MAPPING CHROMOSOMES #####
Chrom <- rep(1, 2)
Start <- c(0, SP$centromere)
End <- c(232, SP$centromere)
Name <- c("contig", "centromere")

df <- data.frame(Chrom, Start, End, Name)
df
chromPlot(gaps = df)


##### CROMPLOT EXAMPLE DATA #####
data(hg_gap)
head(hg_gap)

chromPlot(gaps = hg_gap)
