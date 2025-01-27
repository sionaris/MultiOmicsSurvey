cosmic_pre = read.csv("Resources/COSMIC_CGC_Breast_pre.csv")
cosmic = cosmic_pre[grep("breast", cosmic_pre$Tumour.Types.Somatic.), ]
write.csv(cosmic, "Resources/COSMIC_CGC_Breast_somatic.csv")
