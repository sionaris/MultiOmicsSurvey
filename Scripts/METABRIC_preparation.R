# A script for the preprocessing of METABRIC data for validation purposes

# Libraries #####
library(data.table)
library(readxl)
library(dplyr)

# Import #####

# Import molecular data from the corresponding Resources subdirectory
METABRIC_CNV = fread("Resources/METABRIC/METABRIC_data_CNA.txt")
METABRIC_RNA = fread("Resources/METABRIC/METABRIC_data_Expression.txt")
METABRIC_MUT = fread("Resources/METABRIC/METABRIC_mutation_data.txt")

# Clinical data
clin1 = fread("Resources/METABRIC/TableS6.txt")
clin2 = read.csv("Resources/METABRIC/TableS6_1.csv")
clin3 = fread("Resources/METABRIC/TableS7.txt")
full_index = read_excel("Resources/METABRIC/TableS8.xls", sheet = 3)

# Processing #####
# Filter clinical data
clindata = clin3 %>% 
  dplyr::select(Sample.ID = METABRIC.ID, Normal.Sample.ID = MATCHED.NORMAL.METABRIC.ID,
                Cohort, Age.At.Diagnosis, Date.Of.Diagnosis, Last.Followup.Status,
                NPI, ER.Status, Lymph.Nodes.Positive, Chemo = CT, Endo = HT, Radio = RT,
                Grade, Size, Stage, Histological.Type, Bone.metastasis = BONE,
                Breast.distant.metastasis = BREAST, Liver.metastasis = LIVER,
                Lymph.node.metastasis = LNS, Ovary.metastasis = OVARY, Pleura.metastasis = PLEURA,
                Ascites.metastasis = ASCITES, Skin.metastasis = SKIN, Lung.metastasis = LUNG,
                Pericardial.metastasis = PERICARDIAL, Soft.tissues.metastasis = SOFT_TISSUES,
                Brain.metastasis = BRAIN, Mediastinal.metastasis = MEDIASTINAL, 
                Meningeal.metastasis = MENINGEAL, Adrenal.metastasis = ADRENALS,
                Abdomen.metastasis = ABDOMEN, GI.metastasis = GI_TRACT, Bladder.metastasis = BLADDER,
                Eye.metastasis = EYE, Gallbladder.metastasis = GALLBLADDER,
                Uterus.and.fallopian.tube.metastasis = UTERUS_AND_FALLOPIAN_TUBE,
                Kidney.metastasis = KIDNEY, Peritoneum.metastasis = PERITONEUM,
                Pancreas.metastasis = PANCREAS, Other.metastasis = OTHER,
                Unspecified.metastasis = UNSPECIFIED, Death.of.breast.cancer = DeathBreast,
                Death, Days.to.death.or.last.followup = `T`, 
                Days.to.first.locoregional.relapse.or.last.followup = TLR,
                Locoregional.relapse = LR,
                Days.to.first.distant.relapse.or.last.followup = TDR,
                Distant.relapse = DR, Relapse.type = TYPE.RELAPSE,
                Relapse.site = SITE, Relapse.time = TIME.RELAPSE) %>%
  distinct()

# Initiate the object that will be exported
input_METABRIC = list()
