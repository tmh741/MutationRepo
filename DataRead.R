library(tidyverse)
library(maftools)
library(BiocManager)

#Read data
test <- read.maf("TCGA.ESCA.mutect.7f8e1e7c-621c-4dfd-8fad-af07c739dbfc.DR-10.0.somatic.maf.gz")

#Visually summarize data
plotmafSummary(maf=test, rmOutlier=T,addStat="median",dashboard=T,titvRaw=F)

#Look at alleles
test@data %>% group_by(Allele) %>% summarize(n=n())

# Filter out SNV mutations in data.
mutations <- test@data %>% filter(VARIANT_CLASS=="SNV") %>%
  select(Chromosome, Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,HGVSc,Allele,TRANSCRIPT_STRAND) %>% 
  mutate(Mutation = gsub('[^A-Z]|>', '', HGVSc)) %>% ungroup()

# Filter out SNV mutations for silent mutations.
mutations2 <- test@maf.silent %>% filter(VARIANT_CLASS=="SNV") %>% 
  select(Chromosome, Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,HGVSc,Allele,TRANSCRIPT_STRAND) %>% 
mutate(Mutation = gsub('[^A-Z]|>', '', HGVSc)) %>% ungroup()

# Combine both, fill in blanks, and subtitute to put them in form.
mutationsboth <- rbind(mutations,mutations2)
mutationsboth[which(mutationsboth$Mutation==""),]$Mutation <- paste(mutationsboth[which(mutationsboth$Mutation==""),]$Tumor_Seq_Allele1,mutationsboth[which(mutationsboth$Mutation==""),]$Tumor_Seq_Allele2,sep="")
mutationsboth[mutationsboth$Mutation=="AG",]$Mutation <- "TC"
mutationsboth[mutationsboth$Mutation=="AC",]$Mutation <- "TG"
mutationsboth[mutationsboth$Mutation=="AT",]$Mutation <- "TA"
mutationsboth[mutationsboth$Mutation=="GC",]$Mutation <- "CG"
mutationsboth[mutationsboth$Mutation=="GA",]$Mutation <- "CT"
mutationsboth[mutationsboth$Mutation=="GT",]$Mutation <- "CA"

# Create a list in case there are some weird ones remaining to be filtered out.
list = c("TG","TA","TC","CT","CG","CA","AG","AC","AT","GC","GA","GT")

# Put them in table & calculate proportions.
mutationstable <- mutationsboth %>% 
  group_by(Mutation) %>%
  summarize(n=n()) %>% ungroup()
totalval = sum(mutationstable$n)
mutationstable$nadjust = mutationstable$n/totalval

mutationstable %>% ggplot() + aes(x=Mutation,y=nadjust) + geom_bar(stat="identity") + 
  geom_text(aes(label=n),position=position_dodge(width=1),vjust=-0.25) +
  ylim(y=c(0,1))

mutations %>% group_by(Mutation) %>% summarize(n=n()) %>% ungroup()



mutationstable <- mutationsboth %>% 
  group_by(Mutation,Chromosome) %>%
  summarize(n=n()) %>% ungroup() %>% group_by(Chromosome) %>%  
  mutate(tot = sum(n),nadjust = round(n/tot,2)) %>% ungroup()

mutationstable %>% ggplot() + aes(x=Mutation,y=nadjust) + geom_bar(stat="identity",aes(fill=tot)) + 
  geom_text(aes(label=nadjust),position=position_dodge(width=1),vjust=-0.25) +
  ylim(y=c(0,1)) + facet_wrap(~Chromosome)

