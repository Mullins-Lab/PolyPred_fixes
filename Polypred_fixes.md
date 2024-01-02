## Protocol for testing fine-mapping-informed polygenic risk scores via PRS-CS and PolyPred
Authors:
Maria Koromina maria.koromina@mssm.edu
Niamh Mullins niamh.mullins@mssm.edu

### Materials:
PGC3 BD GWAS and fine-mapping results (formatted for PRS testing) all in build hg19: 

- Formatted PRS-CS results of PGC3 BD GWAS: pgc3_bip_forprscs.txt

- Genome-wide fine-mapping results of PGC3 BD GWAS via SuSiE: pgc3_bip_susie_hrc.txt.gz

- Genome-wide fine-mapping results of PGC3 BD GWAS via Polyfun SuSiE: pgc3_bip_polyfun_susie_hrc.txt.gz

Link to download these files (checksums included) can be found [here](https://www.dropbox.com/sh/8lkzx1mmapik0ms/AABoRcqclvrgeJ9p8gzSF9Oca?dl=0). 

Genotypes of target dataset: We used best guess genotypes (imputed allele dosages converted to hard calls) for this analysis, with QC as appropriate for GWAS. (For data processed with Ricopili use .bgn files). 

Covariates for target dataset: As per a GWAS, please include the relevant principal components for the target dataset (and factors capturing site of recruitment/ genotyping batch if necessary). Please do not include age and sex.

PRS-CS software: [installation instructions](https://github.com/getian107/PRScs).

LD reference panel: LD reference panel of the appropriate ancestry provided with PRS-CS can be found [here](https://github.com/getian107/PRScs). 

PolyPred software: [installation instructions](https://github.com/omerwe/polyfun).

r2redux software: [usage instructions](https://github.com/mommy003/r2redux). 

### Methods:
### 1. Polygenic risk scoring via PRS-CS method
### 1.1 Calculate weights from PGC BD GWAS summary statistics using PRS-CS

Example script:
```
python PRScs.py --ref_dir=/path/to/ldblk_ukbb_eur \
--bim_prefix=/path/to/targetdatasetname \
--sst_file= pgc3_bip_forprscs.txt \
--n_gwas=411551 \ #Full sample size for PGC3 BD GWAS
--phi=1e-2 \
--out_dir=targetdatasetname_prscsweights #should include ~1 million SNPs
```


### 1.2 Calculate PRS for individuals in the target dataset using PLINK

Example script:
```
plink --bfile /path/to/targetdatasetname \
--score targetdatasetname_prscsweights \
--out targetdatasetname_forprscs_scores
```

### 1.3 Test association between PRS and BD vs control status in the target dataset
Below is R code as an example. 

```
install.packages("fmsb")
library(fmsb)

###TEST PRS IN BD TARGET COHORT
a=read.table("targetdatasetname_forprscs_scores.prs", h=T)  

#PHENO FILE FOR TARGET COHORT
b=read.table("targetdatasetname_full_pheno.txt", h=T)

#COV FILE FOR TARGET COHORT
c=read.table("targetdatasetname_covfile", h=T)

d=merge(a,b,by=c("IID"))
x=merge(d,c,by=c("IID"))

table(x$PHENO)
names(x)[5]="BD"
x=subset(x, !(is.na(BD)))
table(x$BD) #check n cases and controls here
x$BD=as.factor(x$BD) #must be a factor with two levels, this code assumes control status is the baseline 
x$PRS<-scale(x$PRS) #standardize PRS to mean 0 and sd 1
#TEST PRS 
model1<-glm(BD ~ PRS +C1 +C2 +C3 +C4 +C5, x, family=binomial) #(logit link)
summary(model1)
model2<-glm(BD ~ C1 +C2 +C3 +C4 +C5, x, family=binomial)
summary(model2)
anova(model1, model2, test="Chi")$P  #extract P value for PRS 
```

### 1.4 Please report all the following results:

``PRS`` 	``Target`` ``cohort name``	``N SNPs in PRS``	``Target cohort Ncases``	``Target cohort Ncontrols``	``Beta of PRS on case status``	``P value of PRS``

#for PRS-CS						


### 2. Polygenic risk scoring via PolyPred method
Note section 2 needs to be done twice – once using SuSie fine-mapping and once using Polyfun-SuSie.

### 2.1 File preparation
Note: In this step, we match the orientation of A1 and A2 alleles in the GW fine-mapping and PRS-CS weights file (separately!) with the ones in the target cohort .bim file for PolyPred to run correctly. 

Below is some hacky R code for orienting the alleles to match the .bim file. (Feel free to use your own): 

```
#align A1 and A2 between GW fine-mapping file and bim file

x=read.delim(gzfile("pgc3_bip_polyfun_susie_hrc.txt.gz"), h=T) #GW fine-mapping file
table(duplicated(x$SNP)) 

y=read.table("targetdatasetname.bim", h=F)
names(y)[1]="CHR"
names(y)[2]="SNP"
names(y)[4]="BP"
names(y)[5]="A1"
names(y)[6]="A2"

z=merge(x,y, by=c("SNP"))
table(duplicated(z$SNP)) 

##check alleles
z$A1A1=z$A1.x==z$A1.y
z$A2A2=z$A2.x==z$A2.y
table(z$A1A1, z$A2A2)

z=z[,c(2,1,3:13)]

names(z)[1]="CHR"
names(z)[3]="BP"
names(z)[4]="A1"
names(z)[5]="A2"
names(z)[11]="BETA"

ok <- z[ which(z$A1A1=="TRUE" & z$A2A2=="TRUE"), ] #SNPs correctly aligned
ok=ok[,c(2,1,3:13)]
names(ok)[1]="CHR"
names(ok)[3]="BP"
names(ok)[4]="A1"
names(ok)[5]="A2"
names(ok)[11]="BETA"

flip<- z[ which(z$A1A1=="FALSE" & z$A2A2=="FALSE"), ] #SNPs for allele switching
flip=flip[,c(2,1,3:13)]
names(flip)[1]="CHR"
names(flip)[3]="BP"
names(flip)[5]="A1" #switch alleles
names(flip)[4]="A2" #switch alleles
flip$BETA.x=flip$BETA.x*-1 #recalculate beta for other allele
names(flip)[11]="BETA"
flip=flip[,c(1:3,5,4,6:13)]

z=rbind(ok, flip) #join all SNPs
dim(z) 

write.table(z,"polyfun_susie_hrc_targetdatasetname_prscssnps.allelesaligned.txt", row.names=F, col.names=T, quote=F, sep="\t") # Polyfun-susie results with alleles aligned, ready for use with PolyPred. (~ 7.5 million SNPs)


#align A1 and A2 between prscs weights file and bim file
x=read.table("targetdatasetname.bim", h=F)
names(x)[1]="CHR"
names(x)[2]="SNP"
names(x)[4]="BP"
names(x)[5]="A1"
names(x)[6]="A2"
 
##y file is a merge of all chrom .txt files from PRS-CS with CHR, BP, A2 columns from the GWAS sumstats added
y=read.table("targetdatasetname_prscsweights.allcols", h=T) #PRS-CS weights file
table(duplicated(y$SNP)) 

z=merge(x,y, by=c("SNP"))
table(duplicated(z$SNP)) 

##check alleles
z$A1A1=z$A1.x==z$A1.y
z$A2A2=z$A2.x==z$A2.y
table(z$A1A1, z$A2A2)

z=z[,c(2,1,3:13)]

names(z)[1]="CHR"
names(z)[3]="BP"
names(z)[4]="A1"
names(z)[5]="A2"
names(z)[11]="BETA"

ok <- z[ which(z$A1A1=="TRUE" & z$A2A2=="TRUE"), ] #SNPs correctly aligned
ok=ok[,c(2,1,3:13)]
names(ok)[1]="CHR"
names(ok)[3]="BP"
names(ok)[4]="A1"
names(ok)[5]="A2"
names(ok)[11]="BETA"

flip<- z[ which(z$A1A1=="FALSE" & z$A2A2=="FALSE"), ] #SNPs for allele switching
flip=flip[,c(2,1,3:13)]
names(flip)[1]="CHR"
names(flip)[3]="BP"
names(flip)[5]="A1" #switch alleles
names(flip)[4]="A2" #switch alleles
flip$BETA.x=flip$BETA.x*-1 #recalculate beta for other allele
names(flip)[11]="BETA"
flip=flip[,c(1:3,5,4,6:13)]

z=rbind(ok, flip) #join all SNPs
dim(z) 

write.table(z, " targetdatasetname_prscsweights.allcols.allelesaligned.txt", row.names=F, col.names=T, quote=F, sep="\t"). #prscs weights file with alleles aligned, ready for use with PolyPred. (~ 1 million SNPs)
```

### 2.2 Mix weights from PRS-CS and fine-mapping using PolyPred
In this step, we will mix the weights (effect sizes) from fine-mapping (either SuSiE or Polyfun-SuSiE) and PRS-CS using PolyPred to train the optimal PRS weights. 

The .pheno file below should be a plink style pheno file including only the unrelated cases and controls for the PRS analysis (related individuals may remain in the plink binary file. 

Important note: the PolyPred.py script will not run if there are duplicate variants in the target cohort, so these will need to be removed from the plink files (.bed, .bim).

Make sure to activate the ‘polyfun’ conda environment: ``` conda activate polyfun ```
Example script: 
```
python PolyPred.py \
--combine-betas \
--betas polyfun_susie_hrc.targetdatasetname_prscssnps.allelesaligned.txt, targetdatasetname_prscsweights.allcols.allelesaligned.txt \
--pheno targetdatasetname_full_pheno.txt \
--output-prefix combine_targetdatasetname_polyfun_susie_prscs \
--plink-exe plink targetdatasetname_nodups.bed
```

### 2.3 Compute PRS with these mixed weights using procedures in 1.2.

### 2.4 Test association between PRS and BD vs control status in the target using procedures in 1.3.

### 2.5 Report results as in Table in 1.4, and additionally the weight that was attributed to each input file by PolyPred, which can be found in the targetdatasetname_polyfun_susie_prscs.prsfile.


```PRS``` 	```Target cohort name```	```N SNPs in PRS```	```Target cohort Ncases```	```Target cohort Ncontrols```	```Beta of PRS on case status```	```P value of PRS```	```Mix weights Fine-mapping (PolyPred only)```	```Mix weights PRS-CS (PolyPred only)```
#for PRS-CS						

#for PolyPred (SuSiE)								

#for PolyPred (Polyfun-SuSiE)								


### 3. Use r2redux software to calculate R2, 95%CI for each R2 and pvalues for the R2diff between PRS-CS and PolyPred based method

### 3.1 Slightly format your PRS results files using `awk` or `grep` (PolyPred, PRS-CS) in order to include the IIDs and the Score columns.

### 3.2 Use the following code in order to calculate R20, and S.E and pvalues for the R2diff between PRSCS and PolyPred based method

```
library(tidyverse)
library(dplyr)
library(r2redux)

#load prscs results for one cohort
prs_cs <- read.table('targetcohort_r2redux_prscs.txt', h=T)
head(prs_cs)
dim(prs_cs)

#load prscs polyfun-susie results for one cohort
polyfun_susie <- read.table('targetcohort_r2redux_prscs_polyfun_susie_pgc3.txt', h=T)
head(polyfun_susie)
dim(polyfun_susie)

#load prscs susie results for one cohort
susie <- read.table('targetcohort_r2redux_prscs_susie_pgc3.txt', h=T)
head(susie)
dim(susie)

#rename the PRS columns
colnames(susie)[2] = "PRS_SUSIE"
colnames(polyfun_susie)[2] = "PRS_POLYFUN_SUSIE"
colnames(prs_cs)[2] = "PRS_PRSCS"

#load the pheno file
pheno <- read.table('targetdatasetname_full_pheno.txt', h=T)

#load the covariate file
cov <- read.table(' targetdatasetname_covfile', h=T)

##adjust Phenos for PCs and residualise
merge_cov_pheno <- merge(cov, pheno, by="IID")
merge_cov_pheno$PHENO=as.factor(merge_cov_pheno$PHENO)
pheno_pc_model <- glm(PHENO ~ C1 + C2 + C3 + C4 + C5, merge_cov_pheno, family=binomial)
merge_cov_pheno$corrected_pheno <- residuals(pheno_pc_model)

##merge the dataframes in rounds of two using the IID column
merge_1 <- merge(merge_cov_pheno, prs_cs, by="IID")
merge_2 <- merge(merge_1, susie, by="IID")
merge_3 <- merge(merge_2, polyfun_susie, by="IID")

##standardize PRS to mean 0 and sd 1
merge_3$PRS_PRSCS <- scale(merge_3$PRS_PRSCS)
merge_3$PRS_SUSIE <- scale(merge_3$PRS_SUSIE)
merge_3$PRS_POLYFUN_SUSIE <- scale(merge_3$PRS_POLYFUN_SUSIE)

#prepare the dataframe for r2redux
filter_targetcohort <- merge_3 %>% select(corrected_pheno, PRS_PRSCS, PRS_SUSIE, PRS_POLYFUN_SUSIE)

#calculate R2 in PRS pairs per cohort
nv=length(filter_for2107$corrected_pheno)
v1=c(1)
v2=c(2)
output=r2_diff(filter_targetcohort, v1, v2, nv)
print(output)
write.table(output, targetcohort_r2redux_prscs_susie.txt')

nv=length(filter_targetcohort $corrected_pheno)
v1=c(1)
v2=c(3)
output2=r2_diff(filter_targetcohort, v1, v2, nv)
print(output2)
write.table(output2,'targetcohort_r2redux_prscs_polyfun_susie.txt'). ##the output includes R20 values for the PRS methods, SE for each method, and p values for the R2 diff 
```

### 3.2 Convert the R20 into R2 liability scale

```
R2O= output$rsq2 #R2 from output 3.1 
K = 0.02 #prevalence of BD in the population – we assume 2% 
P=Ncases/(Ncases+Ncontrols) #prevalence of disorder in the target dataset 
thd = -qnorm(K,0,1) #the threshold on the normal distribution which truncates the proportion of disease prevalence 
zv = dnorm(thd) #z (normal density) 
mv = zv/K #mean liability for case 
theta = mv*(P-K)/(1-K)*(mv*(P-K)/(1-K)-thd) #theta in equation  
cv = K*(1-K)/zv^2*K*(1-K)/(P*(1-P)) #C in equation  
R2 = R2O*cv/(1+R2O*theta*cv)*100 
R2 #R2 on the liability scale

##use the var2 (SE) from the 3.1 outputs to compute 95%CI for the R2 liab scale values
Lower limit = -1.96*(var2*100)+ R2 liab
Upper limit = 1.96*(var2*100) + R2 liab
```

#use the r2_diff function from r2 redux to convert SE and 95%CI to the liab scale too.

#results table should like:
``PRS`` 	``Target cohort name``	``R2 liability for PRS (not as %)``	``95%CI of R2 liab``	``P value of R2 diff`` 	``N SNPs in PRS``	``Target cohort Ncases``	``Target cohort Ncontrols``	``Beta of PRS on case status``	``P value of PRS``	``Mix weights Fine-mapping (PolyPred only)``	``Mix weights PRS-CS (PolyPred only)``

#for PRS-CS			                                                                							                                                                                                                

#for PolyPred (SuSiE)											

#for PolyPred (Polyfun-SuSiE)											



