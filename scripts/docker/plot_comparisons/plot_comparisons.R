library(reshape2)
library(plyr)
library(ggplot2)
library(dplyr)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

args <- commandArgs(trailingOnly = TRUE)
comparison_list <- args[1]
happy_fof <- args[2]
sv_fof <- args[3]
output_prefix <- args[4]

comparisons <- read.table(comparison_list, sep="\t",
                          header=FALSE,
                          col.names=c("Sample1", "Sample2", "ComparisonType"))

happy_files <- read.table(happy_fof, header=FALSE, col.names=c("file"))
sv_files <- read.table(sv_fof, header=FALSE, col.names=c("file"))
happy <- data.frame()
for (h in happy_files$file) {
    df <- read.table(file=h, header=TRUE, sep=",")
    happy <- rbind(happy, df)
}

sv <- data.frame()
for (s in sv_files$file) {
    df <- read.table(file=s, header=TRUE, sep="\t")
    sv <- rbind(sv, df)
}
happy$denominator <- happy$TRUTH.TP + happy$TRUTH.FN + happy$QUERY.FP
happy$numerator <- happy$TRUTH.FN + happy$QUERY.FP
happy$Rate <- happy$numerator/happy$denominator

sv_cast <- dcast(sv, Sample1+Sample2+Subset~Class, value.var="Count")
sv_cast$Type <- "SV"
sv_cast$Filter <- "PASS"
sv_cast$Subtype <- "*"
sv_cast$Genotype <- "*"
sv_cast$denominator <- sv_cast$match + sv_cast$discordant + sv_cast$match_discordant_type + sv_cast$discordant_discordant_type + sv_cast$"0-only" + sv_cast$"1-only"
sv_cast$numerator <- sv_cast$discordant + sv_cast$discordant_discordant_type + sv_cast$"0-only" + sv_cast$"1-only"
sv_cast$Rate <- sv_cast$numerator/sv_cast$denominator

all <- rbind(happy[,c("Sample1", "Sample2", "Subset", "Type", "Filter", "Subtype", "Genotype", "Rate")],
             sv_cast[,c("Sample1", "Sample2", "Subset", "Type", "Filter", "Subtype", "Genotype", "Rate")])
all <- merge(all, comparisons)

all$Type <- revalue(all$Type, c("SNP"="SNV", "INDEL"="indel"))
all$Type <- factor(all$Type, levels=c("SNV", "indel", "SV"))
all$Subset <- factor(all$Subset, levels=c("easy", "medium", "hard", "TS_boundary", "*"))

ggplot(all %>% filter(Subtype=="*" & Filter=="PASS"), aes(x=ComparisonType, y=Rate, fill=ComparisonType)) +
    geom_boxplot() +
    facet_grid(Type~Subset, scales="free") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), panel.background = element_blank()) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    scale_fill_manual(values=cbPalette) +
    ylab("Variant discordance rate") +
    xlab("Comparison type")
ggsave(filename=paste(output_prefix, "png", sep="."))
