library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(ggplot2)

# Professional Pathway Enrichment Workflow
# Mapping differentially expressed gene symbols to Entrez database

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# [RESEARCH SENSITIVE: CUSTOM PATHWAY FILTERING]
# Specialized selection of Biological Process (BP) terms focused 
# exclusively on vascular remodeling and neuro-inflammatory dynamics.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Standard Enrichment Execution (BP, MF, CC)
go_enrich_BP <- enrichGO(gene = entrez_ids, OrgDb = org.Rn.eg.db, ont = "BP", pvalueCutoff = 0.05)

# Visualization for Research Communication
dotplot(go_enrich_BP, showCategory = 20) + 
  ggtitle("VDA Core Pathfinding: Biological Process Enrichment")
  
ggsave("GO_BP_Pathways.png", width = 10, height = 8)
