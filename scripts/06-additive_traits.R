##### UNDERSTANDING ADDITIVE TRAITS IN ALPHASIM_R #####

# 1. INTRODUCTION #####

library(AlphaSimR)
library(tidyverse)
library(gridExtra)
library(chromPlot)
library(ggrepel)
library(MASS)
library(fitdistrplus)

# 2. SIMULATING A POPULATION TO EXPLORE ADDITIVE TRAITS #####
#__ 2.1 Define founder haplotypes #####
pop_haps <- quickHaplo(nInd = 1000,      
                       nChr = 1,      
                       segSites = 100, 
                       genLen = 1,   
                       ploidy = 2L,   
                       inbred = FALSE)

#__ 2.2 Set simulation parameters #####
SP <- SimParam$new(pop_haps)

SP$addTraitA(nQtlPerChr = 10, # number of QTL for trait
             mean = 100, # mean genetic value of the trait
             var = 10, # variance of trait
             corA = NULL, # matrix of correlations between additive effects
             gamma = FALSE,  # to use a gamma distribution in place of a normal distribution
             shape = 1, # shape parameter for gamma distribution only
             force = FALSE) # keep false until this is understood!

SP$setVarE(h2 = 0.5, # narrow sense heritability
           H2 = NULL, # broad sense heritability
           varE = NULL) # error variance 

SP$addSnpChip(nSnpPerChr = 10, # number of SNPs per chromosome
              minSnpFreq = NULL, # minimum allowable frequency for SNP loci
              refPop = NULL) 

#__ 2.3 Create population - generation 0 #####
pop_0 <- newPop(pop_haps, 
                simParam = SP)

# 3. BREEDING, GENETIC, AND PHENOTYPIC VALUES #####
#__ 3.1 Individual level measures #####
pv_gram <- pheno(pop_0)
gv_gram <- gv(pop_0)
bv_gram <- bv(pop_0) + mean(gv(pop_0))
pv_norm <- pheno(pop_0) - mean(pheno(pop_0))
gv_norm <- gv(pop_0) - mean(gv(pop_0))
bv_norm <- bv(pop_0)

ind_vals <- data.frame(pv_gram, gv_gram, bv_gram, pv_norm, gv_norm, bv_norm)
head(round(ind_vals, 2))

plot_gram <- ind_vals %>% ggplot(aes(x = pv_gram, y = bv_gram)) +
                                     geom_point() +
                                     geom_smooth(method='lm', formula= y ~ x)

plot_norm <- ind_vals %>% ggplot(aes(x = pv_norm, y = bv_norm)) +
                                     geom_point() +
                                     geom_smooth(method='lm', formula= y ~ x)

grid.arrange(plot_gram, plot_norm, ncol = 2)

#__ 3.2 Population level measures #####

gram_mod <- lm(bv_gram ~ pv_gram, data = ind_vals)
norm_mod <- lm(bv_norm ~ pv_norm, data = ind_vals)

y_int <- c(coef(gram_mod)[1], coef(norm_mod)[1])
slope <- c(coef(gram_mod)[2], coef(norm_mod)[2])
herit <- c(cov(ind_vals$pv_gram, ind_vals$bv_gram) / var(ind_vals$pv_gram),
           cov(ind_vals$pv_norm, ind_vals$bv_norm) / var(ind_vals$pv_norm))
accur <- c(cor(ind_vals$pv_gram, ind_vals$bv_gram),
           cor(ind_vals$pv_norm, ind_vals$bv_norm))

pop_stats <- data.frame(y_int, slope, herit, accur)
round(pop_stats, 3)

# 4. HERITABILITY AND ADDITIVE TRAITS #####
#__ 4.1 Monte Carlo function observed vs simulated heritability #####
h2_MC <- function(h2){
         pop_haps <- quickHaplo(nInd = 10,      
                                nChr = 1,      
                                segSites = 100, 
                                genLen = 1,   
                                ploidy = 2L,   
                                inbred = FALSE)
        
         SP = SimParam$new(pop_haps)
        
         SP$addTraitA(nQtlPerChr = 10,
                      mean = 100, 
                      var = 10,
                      corA = NULL,
                      gamma = FALSE,
                      shape = 1, 
                      force = FALSE)
        
         SP$setVarE(h2 = h2,
                    H2 = NULL, 
                    varE = NULL)
        
         pop_0 <- newPop(pop_haps, 
                         simParam = SP)
        
         return(list(cov(pop_0@pheno, pop_0@gv) / var(pop_0@pheno),
                     cor(pop_0@pheno, pop_0@gv)))}

#__ 4.2 Monte Carlo analysis #####
mc_results <- lapply(rep(0.5, 1000), h2_MC)
h2 <- unlist(sapply(mc_results, function(x) x[1]))
r2 <- unlist(sapply(mc_results, function(x) x[2]))
mc_df <- data.frame(h2, r2)
head(mc_df)

plot_h2 <- mc_df %>% ggplot() +
                     geom_histogram(aes(x = h2), 
                                    binwidth = 0.05,
                                    color = 'black',
                                    fill = 'blue',
                                    alpha = 0.5) +
                     xlim(c(0, 1))
plot_r2 <- mc_df %>% ggplot() +
                     geom_histogram(aes(x = r2),
                                    binwidth = 0.05,
                                    color = 'black',
                                    fill = "light blue",
                                    alpha = 0.5) +
                     xlim(c(0, 1))

grid.arrange(plot_h2, plot_r2, nrow = 2)

mc_df %>% ggplot(aes(x = h2, y = r2)) +
  geom_point()



?fitdist
fit <- fitdist(mc_df$h2, "norm")
fit
plot(fit, las = 1)

#__ 4.3 heritability effects on accuracy #####
h.9 <- data.frame(h2_mean = mean(mc_df$h2),
                h2_sd = sd(mc_df$h2),
                h2_cv = sd(mc_df$h2) / mean(mc_df$h2),
                h2_var = var(mc_df$h2),
                r2_mean = mean(mc_df$r2),
                r2_sd = sd(mc_df$r2),
                r2_cv = sd(mc_df$r2) / mean(mc_df$r2),
                r2_var = var(mc_df$r2))

df1 <- round(rbind(h.1, h.2, h.3, h.4, h.5, h.6, h.7, h.8, h.9),2)
df1 %>% ggplot(aes(x = h2_mean, y = r2_mean)) +
  geom_point()

# 5. QTL MAPPING #####
#__ 5.1 make a chromosome template #####
MB_to_cM <- 1

c1 <- unlist(pop_haps@genMap[1]) * 100 * (MB_to_cM * 1000000)
c2 <- unlist(pop_haps@centromere) * 100 * (MB_to_cM * 1000000)

Chrom <- rep(1, 2)
Start <- c(0, c2) 
End <- c(max(c1), c2)
Name <- c("contig", "centromere")

at_map <- data.frame(Chrom, Start, End, Name)
chromPlot(gaps = at_map)

#__ 5.2 Add SNP as genomic elements #####

c3 <- getSnpMap(snpChip = 1,
                gender = "A",
                simParam = SP)
c4 <- c3$pos * 100 * (MB_to_cM * 1000000)
c4

Chrom <- rep(1, length(c4))
Start <- c4
End <- c4 + 1000000
Name <- paste("SNP", 1:length(c4), sep = "_")
Colors <- rep("red", length(c4))

at_map_1 <- data.frame(Chrom, Start, End, Name, Colors)

chromPlot(gaps = at_map,
          bands = at_map_1)
at_map_1

#__ 5.3 add QTL as labels #####

c5 <- getQtlMap(trait = 1,
                gender = "A",
                simParam = SP)
c6 <- c5$pos * 100 * (MB_to_cM * 1000000)
c6

Chrom <- rep(1, length(c6))
Start <- c6
End <- c6 + 1
ID <- paste("QTL", 1:length(c6), sep = "_")
Colors <- rep("darkgreen", length(c6))

AIMS <- data.frame(Chrom, Start, End, ID, Colors)
head(AIMS)

chromPlot(gaps = at_map,
          bands = at_map_1,
          stat = AIMS,
          statCol = "Value")

# 6. INDIVIDUAL QTL EFFECTS #####
#__ 6.1 Pull genotypes, format into matrix, and solve linear model #####
qgeno <- pullQtlGeno(pop = pop_0,
                     trait = 1,
                     chr = 1,
                     simParam = SP) #
QTL <- matrix(qgeno,
              nrow = nrow(qgeno),
              ncol = ncol(qgeno),
              byrow = FALSE)

#__ 6.2 Genetic value model #####

d1 <- lm(pop_0@gv ~ QTL) # fit model

d2 <- rep(as.vector(d1$coefficients[1]), nrow(QTL)) # mean of the model
d3 <- as.vector(d1$coefficients[2:(ncol(QTL) + 1)]) # coefficients of model
d3[is.na(d3)] <- 0
d4 <- rowSums(t(t(QTL) * d3)) + d2 # model solutions
d5 <- data.frame(actual = pop_0@gv, predicted = d4) # compare with actual
d5

#__ 6.3 Phenotype model #####

e1 <- lm(pop_0@pheno ~ QTL) # fit model

e2 <- rep(as.vector(e1$coefficients[1]), nrow(QTL)) # mean of the model
e3 <- as.vector(e1$coefficients[2:(ncol(QTL) + 1)]) # coefficients of model
e3[is.na(e3)] <- 0
e4 <- rowSums(t(t(QTL) * e3)) + e2 # model solutions
e5 <- data.frame(actual = pop_0@pheno, predicted = e4) # compare with actual
e5

#__ 6.4 Comparing coefficients - why are they different? #####
f1 <- data.frame(g_coef = d3, p_coef = e3)
f1

f1 %>% ggplot(aes(x = p_coef, y = g_coef)) +
                  geom_point() +
                  geom_smooth(method='lm', formula= y ~ x)

cov(f1$g_coef, f1$p_coef) / var(f1$p_coef)
cor(f1$g_coef, f1$p_coef)

# 7. GENOME WIDE ASSOCIATION, A MAP OF QTL EFFECTS #####
#__ 7.1 Call genotypes at each loci #####
gwa_1 <- pullSegSiteGeno(pop = pop_0, # Pop-class object
                        chr = 1, # chromosome number, NULL = all
                        simParam = SP) # SimParam object
#__ 7.2 Compare mean phenotypic values as a function genotype #####
gwa_2 <- data.frame(loci = as.factor(gwa_1[,1]), 
                 pheno = pheno(pop_0))
gwa_plot_1 <- gwa_2 %>%  ggplot(aes(y = pheno, x = loci)) +
                                geom_boxplot()

gwa_2 <- data.frame(loci = as.factor(gwa_1[,66]), 
                    pheno = pheno(pop_0))
gwa_plot_2 <- gwa_2 %>%  ggplot(aes(y = pheno, x = loci)) +
                                geom_boxplot()

grid.arrange(gwa_plot_1, gwa_plot_2, ncol = 2)

#__ 7.3 Genome wide analysis using -log p values #####

gWAS <- function(loci){
              g1 <- data.frame(loci_1 = gwas[,loci], 
                               pheno = pheno(pop_0))
              g2 <- aov(pheno ~ loci_1, 
                        data = g1)
              -log10(summary(g2)[[1]][[1,"Pr(>F)"]])
              }
gwas_1 <- sapply(1:100, FUN = gWAS)
gwas_2 <- data.frame(x = 1:100, y = gwas_1)
gwas_2 %>% ggplot(aes(x = x, y = y)) +
  geom_point() +
  geom_line() +
  geom_label_repel(data = gwas_3,
                   label = gwas_3$x)

labels <- as.integer(1 + c5$pos*100)
labels
gwas_3 <- gwas_2 %>% filter(x %in% labels)

#__ 7.4 Combing QTL effects with chromosome map #####

# MISCELLANEOUS SCRIPTS TO SORT OUT #####


qtl_gen <- function(x, y) rep(c(rep(0, x), rep(1, x), rep(2, x)), y)

qtl <- 7

x <- c(1, cumprod(rep(3, (qtl-1))))
y <- rev(x)

q1 <- mapply(x, FUN = qtl_gen, y = y)

hist(rowSums(q1))



genicVarA(pop_0, simParam = SP) # additive genetic variance
var(pop_0@gv)
var(pop_0@pheno)

varA(pop_0, simParam = SP)
varG(pop_0)
varA(pop_0, simParam = SP) / varG(pop_0)


### where do breeding values come from?
bv <- bv(pop_0,
         simParam = SP) # calculates breeding values for traits
bv
bv_1 <- scale(pop_0@gv)
plot(bv ~ bv_1)
abline(lm(bv ~ bv_1))

cov(bv, bv_1) / var(bv_1)
cor(bv, bv_1)

plot(bv ~ pop@gv)
abline(lm(bv ~ pop@gv))

plot(bv ~ pop@pheno)
abline(lm(bv ~ pop@pheno))

# looking at standard units
pt_scale <- scale(pop_0@pheno)

plot(pt_scale ~ bv_1)
abline(lm(pt_scale ~ bv_1))

cov(pt_scale, bv_1) / var(bv_1)
cor(pt_scale, bv_1)


### OBSERVING QTL ###

# returns QTL locations on chromosome (0-1 centimorgans)
getQtlMap(trait = 1,         # trait number
          gender = "A",      # "A" for average, "F" for female, "M" for male
          simParam = SP)    # SimParam object

# returns alleles of each QTL (0/1 on each diploid chromosome)
pullQtlHaplo(pop = pop,      # Pop-class object
             trait = 1,      # trait number
             haplo = "all",  # "all" or use 1 for males and 2 for females
             chr = 1,        # chromosome number, NULL = all
             simParam = SP)  # SimParam object

# returns the genotype of each QTL (NB: 0/0 = 0, 0/1 and 1/0 = 1, 1/1 = 2)
qgeno <- pullQtlGeno(pop = pop,
                     trait = 1,
                     chr = 1,
                     simParam = SP)

# summing up the QTL of individuals #
qgeno <- as.data.frame(qgeno)
aeff <- qgeno %>%
  mutate(sum = rowSums(.))

head(aeff$sum)
plot(aeff$sum, pop@gv)
abline(lm(pop@gv~aeff$sum))

plot(aeff$sum, pop@pheno)
abline(lm(pop@pheno~aeff$sum))


# 4.5 Comparing QTL effects - working with observed values?

p_max <- c(0, 0, 2, 2, 2, 0, 0, 0, 2, 0) # maximum gv possible
sum(p_max * e3) + e2[1]
p_min <- c(2, 2, 0, 0, 0, 2, 2, 2, 0, 2) # minimum gv possible
sum(p_min * e3) + e2[1]

