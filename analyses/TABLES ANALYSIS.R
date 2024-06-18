library(dplyr)

################################################################################
# MARTA VS TEMERICS DATA
################################################################################
TEmerics_marta <- fread("CHIMERICS_v1.tab")
marta <- fread("chimerics.tab")

# HEAD
TEmerics_head <- TEmerics_marta[Temerics_marta$Condition == "head_AKA-017"]
#495
marta_head <- marta[marta$tissue == "head" & marta$strain == "AKA-017"]
#492
head <- inner_join(TEmerics_head, marta_head, by = c("transcript_trinity", "transcript_stringtie", "type"))
num_matches_head <- nrow(head)
num_matches_head / 495 * 100
# 98.18

# GUT
TEmerics_gut <- TEmerics_marta[Temerics_marta$Condition == "gut_AKA-017"]
#373
marta_gut <- marta[marta$tissue == "gut" & marta$strain == "AKA-017"]
#368
gut <- inner_join(TEmerics_gut, marta_gut, by = c("transcript_trinity", "transcript_stringtie", "type"))
num_matches_gut <- nrow(gut)
num_matches_gut / 373 * 100
# 97.85


################################################################################
# CHIMERATE VS TEMERICS DATA
################################################################################

# PARAQUAT
################################################################################
# BRG_P
################################################################################
BRG_P <- read.table("projects/BRG/chimreads_evidence_FINAL.tsv", sep = "\t", quote = "")
colnames(BRG_P) <-  c("transcript_stringtie", "TE_consensus", "chim_reads", "transcript_ID", "transcript_FPKM")
BRG_P$Condition <- "BRG_P"
BRG_P <- BRG_P %>% relocate(Condition)

################################################################################
# GRG_P
################################################################################
GRG_P <- read.table("projects/GRG/chimreads_evidence_FINAL.tsv", sep = "\t", quote = "")
colnames(GRG_P) <-  c("transcript_stringtie", "TE_consensus", "chim_reads", "transcript_ID", "transcript_FPKM")
GRG_P$Condition <- "GRG_P"
GRG_P <- GRG_P %>% relocate(Condition)

################################################################################
# DRG_P
################################################################################
DRG_P <- read.table("projects/DRG/chimreads_evidence_FINAL.tsv", sep = "\t", quote = "")
colnames(DRG_P) <-  c("transcript_stringtie", "TE_consensus", "chim_reads", "transcript_ID", "transcript_FPKM")
DRG_P$Condition <- "DRG_P"
DRG_P <- DRG_P %>% relocate(Condition)

################################################################################
# FRG_P
################################################################################
FRG_P <- read.table("projects/FRG/chimreads_evidence_FINAL.tsv", sep = "\t", quote = "")
colnames(FRG_P) <-  c("transcript_stringtie", "TE_consensus", "chim_reads", "transcript_ID", "transcript_FPKM")
FRG_P$Condition <- "FRG_P"
FRG_P <- FRG_P %>% relocate(Condition)


# CONTROL
################################################################################
# BRG_C
################################################################################
BRG_C <- read.table("projects/BRG_C/chimreads_evidence_FINAL.tsv", sep = "\t", quote = "")
colnames(BRG_C) <-  c("transcript_stringtie", "TE_consensus", "chim_reads", "transcript_ID", "transcript_FPKM")
BRG_C$Condition <- "BRG_C"
BRG_C <- BRG_C %>% relocate(Condition)

################################################################################
# FRG_C
################################################################################
FRG_C <- read.table("projects/FRG_C/chimreads_evidence_FINAL.tsv", sep = "\t", quote = "")
colnames(FRG_C) <-  c("transcript_stringtie", "TE_consensus", "chim_reads", "transcript_ID", "transcript_FPKM")
FRG_C$Condition <- "FRG_C"
FRG_C <- FRG_C %>% relocate(Condition)

################################################################################
# DRG_C
################################################################################
DRG_C <- read.table("projects/DRG_C/chimreads_evidence_FINAL.tsv", sep = "\t", quote = "")
colnames(DRG_C) <-  c("transcript_stringtie", "TE_consensus", "chim_reads", "transcript_ID", "transcript_FPKM")
DRG_C$Condition <- "DRG_C"
DRG_C <- DRG_C %>% relocate(Condition)

################################################################################
# GRG_C
################################################################################
GRG_C <- read.table("projects/GRG_C/chimreads_evidence_FINAL.tsv", sep = "\t", quote = "")
colnames(GRG_C) <-  c("transcript_stringtie", "TE_consensus", "chim_reads", "transcript_ID", "transcript_FPKM")
GRG_C$Condition <- "GRG_C"
GRG_C <- GRG_C %>% relocate(Condition)

# MERGE
chimeraTE <- rbind(BRG_C, BRG_P, DRG_C, DRG_P, FRG_C, FRG_P, GRG_P, GRG_P)

# TEmerics classification
TEmerics <- read.table("Chimerics_v2.tab", sep = "\t")
TEmerics_BRG_C <- TEmerics[TEmerics$Condition == "BRGC",]
TEmerics_BRG_P <- TEmerics[TEmerics$Condition == "BRGP",]
TEmerics_CRG_C <- TEmerics[TEmerics$Condition == "CRGC",]
TEmerics_CRG_P <- TEmerics[TEmerics$Condition == "CRGP",]
TEmerics_DRG_C <- TEmerics[TEmerics$Condition == "DRGC",]
TEmerics_DRG_P <- TEmerics[TEmerics$Condition == "DRGP",]
TEmerics_ERG_C <- TEmerics[TEmerics$Condition == "ERGC",]
TEmerics_ERG_P <- TEmerics[TEmerics$Condition == "ERGP",]
TEmerics_FRG_C <- TEmerics[TEmerics$Condition == "FRGC",]
TEmerics_FRG_P <- TEmerics[TEmerics$Condition == "FRGP",]
TEmerics_GRG_C <- TEmerics[TEmerics$Condition == "GRGC",]
TEmerics_GRG_P <- TEmerics[TEmerics$Condition == "GRGP",]

# CALCULATIONS

# BRG
tot_BRG_C  <- inner_join(TEmerics_BRG_C, BRG_C, by = c("TE_consensus", "transcript_stringtie"))
num_matches_BRC_C <- nrow(tot_BRG_C)
num_matches_BRC_C / length(BRG_C) * 100
# 95.1
tot_BRG_P  <- inner_join(TEmerics_GRG_P, BRG_P, by = c("TE_consensus", "transcript_stringtie"))
num_matches_BRC_P <- nrow(tot_BRG_P)
num_matches_BRC_P / length(BRG_P) * 100
# 95.95


# DRG
tot_DRG_C  <- inner_join(TEmerics_DRG_C, DRG_C, by = c("TE_consensus", "transcript_stringtie"))
num_matches_DRC_C <- nrow(tot_DRG_C)
num_matches_DRC_C / length(DRG_C) * 100
# 98.48
tot_DRG_P  <- inner_join(TEmerics_DRG_P, DRG_P, by = c("TE_consensus", "transcript_stringtie"))
num_matches_DRC_P <- nrow(tot_DRG_P)
num_matches_DRC_P / length(DRG_P) * 100
# 97.61


# FRG
tot_FRG_C  <- inner_join(TEmerics_FRG_C, FRG_C, by = c("TE_consensus", "transcript_stringtie"))
num_matches_FRC_C <- nrow(tot_FRG_C)
num_matches_FRC_C / length(FRG_C) * 100
# 98.59
tot_FRG_P  <- inner_join(TEmerics_FRG_P, FRG_P, by = c("TE_consensus", "transcript_stringtie"))
num_matches_FRC_P <- nrow(tot_FRG_P)
num_matches_FRC_P / length(FRG_P) * 100
# 99.1


# GRG
tot_GRG_C  <- inner_join(TEmerics_GRG_C, GRG_C, by = c("TE_consensus", "transcript_stringtie"))
num_matches_GRC_C <- nrow(tot_GRG_C)
num_matches_GRC_C / length(GRG_C) * 100
# 95.23
tot_GRG_P  <- inner_join(TEmerics_GRG_P, GRG_P, by = c("TE_consensus", "transcript_stringtie"))
num_matches_GRC_P <- nrow(tot_GRG_P)
num_matches_GRC_P / length(GRG_P) * 100
# 93.51


###############################################################################
# TEMERICS DATA
################################################################################

# BRG
tot_BRG  <- inner_join(TEmerics_BRG_C, TEmerics_BRG_P, by = c("transcript_trinity", "transcript_stringtie", "TE_family"))
num_matches_BRC <- nrow(tot_BRG)
num_matches_BRC / length(TEmerics_BRG_C) * 100
# 98.55

# CRG
tot_CRG  <- inner_join(TEmerics_CRG_C, TEmerics_CRG_P, by = c("transcript_trinity", "transcript_stringtie", "TE_family"))
num_matches_CRC <- nrow(tot_CRG)
num_matches_CRC / length(TEmerics_CRG_C) * 100
# 99.05

# DRG
tot_DRG  <- inner_join(TEmerics_DRG_C, TEmerics_DRG_P, by = c("transcript_trinity", "transcript_stringtie", "TE_family"))
num_matches_DRC <- nrow(tot_DRG)
num_matches_DRC / length(TEmerics_DRG_C) * 100
# 95.59

# ERG
tot_ERG  <- inner_join(TEmerics_ERG_C, TEmerics_ERG_P, by = c("transcript_trinity", "transcript_stringtie", "TE_family"))
num_matches_ERC <- nrow(tot_ERG)
num_matches_ERC / length(TEmerics_ERG_C) * 100
# 99.26

# FRG
tot_FRG  <- inner_join(TEmerics_FRG_C, TEmerics_FRG_P, by = c("transcript_trinity", "transcript_stringtie", "TE_family"))
num_matches_FRC <- nrow(tot_FRG)
num_matches_FRC / length(TEmerics_FRG_C) * 100
# 98.94

# GRG
tot_GRG  <- inner_join(TEmerics_GRG_C, TEmerics_GRG_P, by = c("transcript_trinity", "transcript_stringtie", "TE_family"))
num_matches_GRC <- nrow(tot_GRG)
num_matches_GRC / length(TEmerics_GRG_C) * 100
# 97.79

