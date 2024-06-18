library(ggplot2)
library(data.table)
library(dplyr)
library(DT)
library(cowplot)
library(stringr)
library(tidyr)
library(ggrepel)
library(simplecolors)
library(envalysis)

marta <- fread("chimerics.tab")
marta_gut <- marta[marta$tissue == "gut" & marta$strain == "AKA-017"]

######################################################################################

chimerics <- read.table("Chimerics_v2.tab", sep = "\t")
all_transcripts <- read.table("all_transcripts_expr_1TPM_v2.tab", sep = "\t")

all_transcripts_m <- read.table("all_transcripts_expr_1TPM.tab", sep = "\t")

colnames(all_transcripts_m) <- c("tissue", "strain", "transcript_trinity", 
                                 "transcript_stringtie", "gene", "transcript", 
                                 "expression", "type", "roo_type", "group")
all_transcripts_m <- all_transcripts_m[all_transcripts_m$tissue == "gut" & all_transcripts_m$strain == "AKA-017",]
all_transcripts_m <- merge(all_transcripts_m, 
                         unique(marta_gut[,c("transcript_stringtie",
                                             "transcript_trinity")]), 
                         by=c("transcript_stringtie",
                              "transcript_trinity"), all = T)

all_transcripts_m$comb <- ""
all_transcripts_m <- remove_missing(all_transcripts_m, vars = "type")
all_transcripts_m[all_transcripts_m$type=="chimeric",]$comb <- paste0(all_transcripts_m[all_transcripts_m$type=="chimeric",]$transcript_stringtie,"-",
                                                                  all_transcripts_m[all_transcripts_m$type=="chimeric",]$id)
all_transcripts_m[all_transcripts_m$type=="non-chimeric",]$comb <- paste0(all_transcripts_m[all_transcripts_m$type=="non-chimeric",]$transcript_stringtie,"-",
                                                                      all_transcripts_m[all_transcripts_m$type=="non-chimeric",]$type)



all_transcripts <- merge(all_transcripts, 
                         unique(chimerics[,c("Condition",
                                             "transcript_stringtie",
                                             "transcript_trinity","id")]), 
                         by=c("Condition","transcript_stringtie",
                              "transcript_trinity"), all = T)
all_transcripts$comb <- ""
all_transcripts <- remove_missing(all_transcripts, vars = "type")
all_transcripts[all_transcripts$type=="chimeric",]$comb <- paste0(all_transcripts[all_transcripts$type=="chimeric",]$transcript_stringtie,"-",
                                                                  all_transcripts[all_transcripts$type=="chimeric",]$id)
all_transcripts[all_transcripts$type=="non-chimeric",]$comb <- paste0(all_transcripts[all_transcripts$type=="non-chimeric",]$transcript_stringtie,"-",
                                                                      all_transcripts[all_transcripts$type=="non-chimeric",]$type)

results_M <- data.frame(sample_id = character(),
                        type = character(),
                        class = character(),
                        value = numeric(),
                        percentage = numeric(),
                        yPos = numeric(),
                        stringsAsFactors = FALSE)

results_C <- data.frame(sample_id = character(),
                      type = character(),
                      class = character(),
                      value = numeric(),
                      percentage = numeric(),
                      yPos = numeric(),
                      stringsAsFactors = FALSE)

results_P <- data.frame(sample_id = character(),
                        type = character(),
                        class = character(),
                        value = numeric(),
                        percentage = numeric(),
                        yPos = numeric(),
                        stringsAsFactors = FALSE)

sample_id_C <- c("BRG_C", "CRG_C", "DRG_C", "ERG_C", "FRG_C", "GRG_C")
sample_id_P <- c("BRG_P", "CRG_P", "DRG_P", "ERG_P", "FRG_P", "GRG_P")
tissue <- "gut"

for (s in sample_id_C) {
  # compute the percentage of chimeric transcripts for this strain
  percentage <- all_transcripts %>%
    filter(type == "chimeric" & Condition == s) %>%
    summarize(percentage = round(length(unique(comb))/length(unique(all_transcripts[all_transcripts$Condition==s,]$comb))*100, 1)) %>%
    pull(percentage)
  
  # compute the number of unique chimeric transcripts for this strain
  unique_count_chimerics <- all_transcripts %>%
    filter(type == "chimeric" & Condition == s) %>%
    summarize(unique_count = n_distinct(comb)) %>%
    pull(unique_count)
  
  unique_count_all <- all_transcripts %>%
    filter(Condition == s) %>%
    summarize(unique_count = n_distinct(comb)) %>%
    pull(unique_count)
  
  # add the results to the data frame
  results_C <- results_C %>%
    add_row(sample_id = s, type = "differenceTotalTranscripts", 
            class = "Total transcriptome", 
            value = unique_count_all-unique_count_chimerics, 
            percentage = percentage, 
            yPos = unique_count_all-unique_count_chimerics) %>%
    add_row(sample_id = s, type = "totalChimerics", 
            class = "Total transcriptome", value = unique_count_chimerics, 
            percentage = NA, yPos = unique_count_all-unique_count_chimerics)
}

for (s in sample_id_P) {
  # compute the percentage of chimeric transcripts for this strain
  percentage <- all_transcripts %>%
    filter(type == "chimeric" & Condition == s) %>%
    summarize(percentage = round(length(unique(comb))/length(unique(all_transcripts[all_transcripts$Condition==s,]$comb))*100, 1)) %>%
    pull(percentage)
  
  # compute the number of unique chimeric transcripts for this strain
  unique_count_chimerics <- all_transcripts %>%
    filter(type == "chimeric" & Condition == s) %>%
    summarize(unique_count = n_distinct(comb)) %>%
    pull(unique_count)
  
  unique_count_all <- all_transcripts %>%
    filter(Condition == s) %>%
    summarize(unique_count = n_distinct(comb)) %>%
    pull(unique_count)
  
  # add the results to the data frame
  results_P <- results_P %>%
    add_row(sample_id = s, type = "differenceTotalTranscripts", 
            class = "Total transcriptome", 
            value = unique_count_all-unique_count_chimerics, 
            percentage = percentage, 
            yPos = unique_count_all-unique_count_chimerics) %>%
    add_row(sample_id = s, type = "totalChimerics", 
            class = "Total transcriptome", value = unique_count_chimerics, 
            percentage = NA, yPos = unique_count_all-unique_count_chimerics)
}

for (s in tissue) {
  # compute the percentage of chimeric transcripts for this strain
  percentage <- all_transcripts_m %>%
    filter(type == "chimeric" & tissue == s) %>%
    summarize(percentage = round(length(unique(comb))/length(unique(all_transcripts_m[all_transcripts_m$tissue==s,]$comb))*100, 1)) %>%
    pull(percentage)
  
  # compute the number of unique chimeric transcripts for this strain
  unique_count_chimerics <- all_transcripts_m %>%
    filter(type == "chimeric" & tissue == s) %>%
    summarize(unique_count = n_distinct(comb)) %>%
    pull(unique_count)
  
  unique_count_all <- all_transcripts_m %>%
    filter(tissue == s) %>%
    summarize(unique_count = n_distinct(comb)) %>%
    pull(unique_count)
  
  # add the results to the data frame
  results_M <- results_M %>%
    add_row(sample_id = s, type = "differenceTotalTranscripts", 
            class = "Total transcriptome", 
            value = unique_count_all-unique_count_chimerics, 
            percentage = percentage, 
            yPos = unique_count_all-unique_count_chimerics) %>%
    add_row(sample_id = s, type = "totalChimerics", 
            class = "Total transcriptome", value = unique_count_chimerics, 
            percentage = NA, yPos = unique_count_all-unique_count_chimerics)
}



results_C$class <- factor(results_C$class, 
                          levels=c("Total transcriptome", 
                                   "Body-part specific transcriptome"))
results_C$type <- factor(results_C$type, 
                         levels=c("totalChimerics", 
                                  "differenceTotalTranscripts"))
results_P$class <- factor(results_P$class, 
                          levels=c("Total transcriptome", 
                                   "Body-part specific transcriptome"))
results_P$type <- factor(results_P$type, 
                         levels=c("totalChimerics", 
                                  "differenceTotalTranscripts"))
results_M$class <- factor(results_M$class, 
                          levels=c("Total transcriptome", 
                                   "Body-part specific transcriptome"))
results_M$type <- factor(results_M$type, 
                         levels=c("totalChimerics", 
                                  "differenceTotalTranscripts"))

results_M$percentage <- c(3.6, NA)
results_M$value <- c(7502, 272)
results_M$yPos <- c(7502,7502)
results_C <- rbind(results_C, results_M)

######################################################################################
### FIGURE
a <- ggplot(results_C[results_C$class == "Total transcriptome",], 
            aes(x = sample_id, y = value, fill = type)) + 
  geom_col() + 
  theme_publish() +
  facet_grid(. ~ class, scales = "free") + 
  geom_text(aes(label = percentage, y = yPos), hjust = 1.1, size = 6) + 
  coord_flip() + 
  scale_fill_manual(values = c("#D51410", "#f2dba9"), guide = "none") +  
  labs(y = "", x = "") + 
  scale_x_discrete(labels = c("AKA-017", "RAL-176", "RAL-855", "RAL-091", 
                              "RAL-059", "RAL-737", "RAL-426")) +
  theme(
    axis.text = element_text(size = 24), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 24), 
    plot.title = element_text(size = 24), 
    strip.text.x = element_blank(),
    strip.background = element_rect(fill = "white", color = "white")
  )
a

b <- ggplot(results_P[results_P$class == "Total transcriptome",], 
            aes(x = sample_id, y = value, fill = type)) + 
  geom_col() + 
  theme_publish() +
  facet_grid(. ~ class, scales = "free") + 
  geom_text(aes(label = percentage, y = yPos), hjust = 1.1, size = 6) + 
  coord_flip() + 
  scale_fill_manual(values = c("#D51410", "#f2dba9"), guide = "none") +  
  labs(y = "Number of transcripts", x = "") + 
  scale_x_discrete(labels = c("RAL-176", "RAL-855", "RAL-091", 
                              "RAL-059", "RAL-737", "RAL-426")) +
  theme(
    axis.text = element_text(size = 24), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 24), 
    plot.title = element_text(size = 24), 
    strip.text.x = element_blank(),
    strip.background = element_rect(fill = "white", color = "white")
  )
b

b <- b + ggtitle("Stress")
a <- a + ggtitle("Control")
combined_plot <- cowplot::plot_grid(a,b,ncol = 1)
title <- ggdraw() + 
  draw_label("Contribution of chimeric gene–TE transcripts to the total transcriptome ", 
             fontface = 'bold', size = 20)
cowplot::plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1, 0.1))

######################################################################################

expression_data <- read.table("all_transcripts_expr_1TPM_v2.tab", sep = "\t")
chimerics_id <- read.table("id_transcripts_v2.tab", sep = "\t")

expression_data_C <- expression_data[
  expression_data$Condition == "BRG_C" |
    expression_data$Condition == "CRG_C" | 
    expression_data$Condition == "DRG_C" |
    expression_data$Condition == "ERG_C" |
    expression_data$Condition == "FRG_C" |
    expression_data$Condition == "GRG_C",]

expression_data_P <- expression_data[
  expression_data$Condition == "BRG_P" |
    expression_data$Condition == "CRG_P" | 
    expression_data$Condition == "DRG_P" |
    expression_data$Condition == "ERG_P" |
    expression_data$Condition == "FRG_P" |
    expression_data$Condition == "GRG_P",]

expression_data_df_C<-merge(expression_data, chimerics_id, 
                            by=c("Condition", "transcript_stringtie"), all.x=T)
expression_data_df_P<-merge(expression_data, chimerics_id, 
                            by=c("Condition", "transcript_stringtie"), all.x=T)

chimerics_pos <- unique(chimerics[,c("Condition", "transcript_trinity", 
                                     "transcript_stringtie", "gene", 
                                     "transcript", "id", "exon"),]) %>% 
  group_by(Condition, transcript_trinity, transcript_stringtie, gene, 
           transcript, id) %>%  
  summarize(exon_pos = str_c(exon, collapse = ","))
chimerics_pos$type<-"chimeric"

expression_data_df_C <- merge(expression_data_df_C, chimerics_pos, 
                              by= c("Condition", "transcript_trinity", 
                                    "transcript_stringtie", "gene", 
                                    "transcript", "id", "type"),all.x=TRUE)
expression_data_df_P <- merge(expression_data_df_P, chimerics_pos, 
                              by= c("Condition", "transcript_trinity", 
                                    "transcript_stringtie", "gene", 
                                    "transcript", "id", "type"),all.x=TRUE)

non_chimerics_expression_P <- expression_data_df_P %>% 
  filter(type == "non-chimeric") %>% 
  select(Condition, transcript_stringtie, expr) %>% unique() 
non_chimerics_expression_P$type <- "Non-chimeric transcripts"
non_chimerics_expression_P<-rename(non_chimerics_expression_P, 
                                   c("id"="transcript_stringtie"))
non_chimerics_expression_C <- expression_data_df_C %>%
  filter(type == "non-chimeric") %>% 
  select(Condition, transcript_stringtie, expr) %>% unique() 
non_chimerics_expression_C$type <- "Non-chimeric transcripts"
non_chimerics_expression_C<-rename(non_chimerics_expression_C, 
                                   c("id"="transcript_stringtie"))

chimerics_expression_P <- expression_data_df_P %>% 
  filter(type == "chimeric") %>% select(Condition, id, expr) %>% 
  unique() 

chimerics_expression_P$type <- "Chimeric transcripts"
chimerics_expression_C <- expression_data_df_C %>% 
  filter(type == "chimeric") %>% select(Condition, id, expr) %>% 
  unique() 
chimerics_expression_C$type <- "Chimeric transcripts"

fig5a_df_C <- rbind(non_chimerics_expression_C, chimerics_expression_C)
fig5a_df_P <- rbind(non_chimerics_expression_P, chimerics_expression_P)

fig5a_df_C$type <- factor(fig5a_df_C$type, 
                          levels=rev(c("Non-chimeric transcripts", 
                                       "Chimeric transcripts")))
fig5a_df_P$type <- factor(fig5a_df_P$type, 
                          levels=rev(c("Non-chimeric transcripts", 
                                       "Chimeric transcripts")))

fig5a_C<-ggplot(fig5a_df_C, aes(x=type, y = log2(expr))) + 
  geom_boxplot(fill=c("firebrick3", "gray80")) + 
  labs(y=expression(log[2](TPM)),x="")  + coord_flip() +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))   + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), text = element_text(size=16), 
        plot.title = element_text(size = 16))
fig5a_C

fig5a_P<-ggplot(fig5a_df_P, aes(x=type, y = log2(expr))) + 
  geom_boxplot(fill=c("firebrick3", "gray80")) + 
  labs(y=expression(log[2](TPM)),x="")  + coord_flip() +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), text = element_text(size=16), 
        plot.title = element_text(size = 16))
fig5a_P

fig5a_P <- fig5a_P + ggtitle("Stress")
fig5a_C <- fig5a_C + ggtitle("Control")
combined_plot <- cowplot::plot_grid(fig5a_C,fig5a_P,ncol = 2)
final_plot <- cowplot::plot_grid(combined_plot, ncol = 1, 
                                 rel_heights = c(1, 0.1))
title <- ggdraw() + draw_label("Boxplots for the expression levels", 
                               fontface = 'bold', size = 20)
cowplot::plot_grid(title, final_plot, ncol = 1, rel_heights = c(0.1, 1, 0.1))
################################################################################
# FIGURE
ids <- chimerics %>% group_by(id) %>%
  summarise(n=length(unique(TE_family)))  %>%
  filter(n==1) %>% pull(id)

chimerics_1TE <- chimerics[chimerics$id %in% ids,]

ids <- chimerics_1TE %>% 
  group_by(id) %>%
  summarise(n=length(unique(exon)))  %>%
  filter(n==1) %>% pull(id)

chimerics_1TE_1pos <- chimerics_1TE[chimerics_1TE$id %in% ids,]

chimerics_1TE_1pos_C <- chimerics_1TE_1pos[
  chimerics_1TE_1pos$Condition == "BRGC" |
    chimerics_1TE_1pos$Condition == "CRGC" |
    chimerics_1TE_1pos$Condition == "DRGC" |
    chimerics_1TE_1pos$Condition == "ERGC" |
    chimerics_1TE_1pos$Condition == "FRGC" |
    chimerics_1TE_1pos$Condition == "GRGC",
]

chimerics_1TE_1pos_P <- chimerics_1TE_1pos[
  chimerics_1TE_1pos$Condition == "BRGP" |
    chimerics_1TE_1pos$Condition == "CRGP" |
    chimerics_1TE_1pos$Condition == "DRGP" |
    chimerics_1TE_1pos$Condition == "ERGP" |
    chimerics_1TE_1pos$Condition == "FRGP" |
    chimerics_1TE_1pos$Condition == "GRGP",
]

frequency <- c("Strain-specific", "Shared across\n2-5 strains")
position <- c("3'UTR", "5'UTR", "Internal exon")
params <- expand.grid(frequency, position)

fig_s3a_df_C <- data.frame(position=as.character(), 
                           chimeric_transcripts=as.numeric())
fig_s3a_df_P <- data.frame(position=as.character(), 
                           chimeric_transcripts=as.numeric())

for (i in seq_len(nrow(params))) {
  freq <- params$Var1[i]
  pos <- params$Var2[i]
  
  if (freq == "Strain-specific") {
    data <- chimerics_1TE_1pos_C %>% 
      group_by(id) %>%
      mutate(n=length(unique(Condition)))  %>%
      filter(n==1)  %>% ungroup() 
    
  } else if (freq == "Shared across\n2-5 strains") {
    data <- chimerics_1TE_1pos_C %>% 
      group_by(id) %>%
      mutate(n=length(unique(Condition)))  %>%
      filter(n>1,n<=5) %>% ungroup() 
    
  }
  
  if (pos == "3'UTR") {
    n <- data %>% filter(exon=="Last exon") %>% summarize(n_distinct(id)) %>% 
      pull()
  } else if (pos == "5'UTR") {
    n <- data %>% filter(exon=="First exon") %>% summarize(n_distinct(id)) %>% 
      pull()
  } else if (pos == "Internal exon") {
    n <- data %>% filter(exon=="Middle exon") %>% summarize(n_distinct(id)) %>% 
      pull()
  }
  
  result <- data.frame(position=pos, chimeric_transcripts=n)
  
  fig_s3a_df_C <- rbind(fig_s3a_df_C, result)  
  
}

fig_s3a_df_C

for (i in seq_len(nrow(params))) {
  freq <- params$Var1[i]
  pos <- params$Var2[i]
  
  if (freq == "Strain-specific") {
    data <- chimerics_1TE_1pos_P %>% 
      group_by(id) %>%
      mutate(n=length(unique(Condition)))  %>%
      filter(n==1)  %>% ungroup() 
    
  } else if (freq == "Shared across\n2-5 strains") {
    data <- chimerics_1TE_1pos_P %>% 
      group_by(id) %>%
      mutate(n=length(unique(Condition)))  %>%
      filter(n>1,n<=5) %>% ungroup() 
    
  }
  
  if (pos == "3'UTR") {
    n <- data %>% filter(exon=="Last exon") %>% summarize(n_distinct(id)) %>% 
      pull()
  } else if (pos == "5'UTR") {
    n <- data %>% filter(exon=="First exon") %>% summarize(n_distinct(id)) %>% 
      pull()
  } else if (pos == "Internal exon") {
    n <- data %>% filter(exon=="Middle exon") %>% summarize(n_distinct(id)) %>% 
      pull()
  }
  
  result <- data.frame(position=pos, chimeric_transcripts=n)
  
  fig_s3a_df_P <- rbind(fig_s3a_df_P, result)  
  
}
fig_s3a_df_P

fig_s3a_df_C$position <- factor(fig_s3a_df_C$position, levels=c("5'UTR", 
                                                                "Internal exon",
                                                                "3'UTR"))
fig_s3a_df_P$position <- factor(fig_s3a_df_P$position, levels=c("5'UTR", 
                                                                "Internal exon",
                                                                "3'UTR"))

f3a_C<-ggplot(fig_s3a_df_C, aes(y=chimeric_transcripts, x = position)) + 
  geom_col(aes(fill=position)) + theme_publish() + 
  scale_fill_manual(values= c("#6495ed", "#6f8faf" , "#5d3fd3")) + 
  labs(x = "", y = "") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), text = element_text(size=16), 
        plot.title = element_text(size = 16))
f3a_P<-ggplot(fig_s3a_df_P, aes(y=chimeric_transcripts, x = position)) + 
  geom_col(aes(fill=position)) + theme_publish() + 
  scale_fill_manual(values= c("#6495ed", "#6f8faf" , "#5d3fd3")) + 
  labs(x = "", y = "") +
  scale_y_continuous(expand = c(0, 0))  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), text = element_text(size=16), 
        plot.title = element_text(size = 16))

f3a_P <- f3a_P + ggtitle("Stress") + theme(legend.position = "none")
legend <- cowplot::get_legend(f3a_C)
f3a_C <- f3a_C + ggtitle("Control") + theme(legend.position = "none")
combined_plot <- cowplot::plot_grid(f3a_C,f3a_P,ncol = 2)
final_plot <- cowplot::plot_grid(combined_plot, legend, ncol = 1, 
                                 rel_heights = c(1, 0.1))
title <- ggdraw() + draw_label("Number of gene–TE chimeric transcripts", 
                               fontface = 'bold', size = 20)
cowplot::plot_grid(title, final_plot, ncol = 1, rel_heights = c(0.1, 1, 0.1))


