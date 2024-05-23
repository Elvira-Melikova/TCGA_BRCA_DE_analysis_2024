# Устанавливаем необходимые пакеты
install.packages(c("tidyverse", "readxl", "writexl", "data.table", "dplyr", "remotes", "msigdbr", "pheatmap", "RColorBrewer"))

# Устанавливаем пакет dplyr, если он еще не установлен
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

# Устанавливаем DESeq2 через BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", force = TRUE)

# Устанавливаем TCGAbiolinks через BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks", force = TRUE)

# Устанавливаем пакеты из GitHub
remotes::install_github("assaron/r-utils")
remotes::install_github("kevinblighe/EnhancedVolcano")

# Подключаем библиотеки
library(dplyr)
library(tidyverse)
library(readxl)
library(writexl)
library(data.table)
library(org.Hs.eg.db)
library(rUtils) 
library(EnhancedVolcano)
library(msigdbr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(TCGAbiolinks)
library(DESeq2)

# Загружаем данные из GDC
proj <- "TCGA-BRCA"
query <- GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data <- GDCprepare(query)

# Извлекаем матрицу счетов
counts <- assay(data)

# Извлекаем метаданные образцов
metadata <- colData(data)

# Извлекаем типы образцов
sample_type <- metadata$sample_type

# Фильтруем данные по типам образцов
filtered_metadata <- as.data.frame(colData(data)) %>%
  filter(sample_type %in% c("Metastatic", "Primary Tumor"))

# Создаем отдельные наборы данных для "Primary Tumor" и "Metastatic"
primary_tumor <- filtered_metadata %>%
  filter(sample_type == "Primary Tumor")

metastatic <- filtered_metadata %>%
  filter(sample_type == "Metastatic")

# Находим общие ID пациентов
common_patients <- intersect(primary_tumor$patient, metastatic$patient)

# Оставляем только образцы "Primary Tumor" с общими ID пациентов
primary_tumor_common <- primary_tumor %>%
  filter(patient %in% common_patients)

# Фильтруем матрицу счетов по типу образца
filtered_counts <- counts[, sample_type %in% c("Primary Tumor", "Metastatic")]

# Выбираем нужные столбцы по их именам
selected_columns <- c("TCGA-BH-A1ES-01A-11R-A137-07", 
                      "TCGA-BH-A1FE-01A-11R-A13Q-07", 
                      "TCGA-AC-A6IX-01A-12R-A32P-07",
                      "TCGA-E2-A15A-01A-11R-A12D-07",
                      "TCGA-BH-A18V-01A-11R-A12D-07",
                      "TCGA-E2-A15E-01A-11R-A12D-07",
                      "TCGA-E2-A15K-01A-11R-A12P-07",
                      "TCGA-BH-A1FE-06A-11R-A213-07",
                      "TCGA-BH-A1ES-06A-12R-A24H-07",
                      "TCGA-AC-A6IX-06A-11R-A32P-07",
                      "TCGA-BH-A18V-06A-11R-A213-07",
                      "TCGA-E2-A15K-06A-11R-A12P-07",
                      "TCGA-E2-A15A-06A-11R-A12D-07",
                      "TCGA-E2-A15E-06A-11R-A12D-07")

filtered_counts <- filtered_counts[, selected_columns]

# Фильтруем метаданные образцов по типу образца и пациентам
metadata <- colData(data)

# Фильтруем данные по пациентам
patients_to_keep <- c("TCGA-BH-A1ES", "TCGA-BH-A1FE", "TCGA-AC-A6IX", "TCGA-E2-A15A", "TCGA-BH-A18V", "TCGA-E2-A15E", "TCGA-E2-A15K")
filtered_df <- rbind(primary_tumor_common, metastatic)
filtered_df <- filtered_df %>%
  filter(patient %in% patients_to_keep)

# Получаем список ID образцов для сортировки
sorted_sample_ids <- colnames(filtered_counts)

# Сортируем датасеты
sorted_filtered_counts <- filtered_counts[, sorted_sample_ids]
sorted_filtered_df <- filtered_df[sorted_sample_ids, ]

# Проверяем совпадение колонок и строк
all(colnames(sorted_filtered_counts) == rownames(sorted_filtered_df))

# Создаем объект DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = sorted_filtered_counts,
                              colData = sorted_filtered_df,
                              design = ~ tumor_descriptor)

# Анализируем контроль качества и нормализуем объект DESeq
dds <- estimateSizeFactors(dds)

# Извлекаем нормализованные количества чтений
normlzd_dds <- counts(dds, normalized = TRUE)
head(normlzd_dds)

# Преобразуем данные для стабилизации дисперсии
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

# Вычисляем попарные корреляции
vsd_cor <- cor(vsd_mat)

# Строим тепловую карту корреляций
pheatmap(vsd_cor)

# Выполняем PCA для проверки кластеризации образцов
pca_plot_tumor <- plotPCA(vsd, intgroup = "tumor_descriptor") +
  geom_point(size = 5) +
  theme(legend.text = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15))
print(pca_plot_tumor)

pca_plot_patient <- plotPCA(vsd, intgroup = "patient") +
  geom_point(size = 5) +
  theme(legend.text = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15))
print(pca_plot_patient)

# Запускаем конвейер DESeq2
dds$patient <- as.factor(dds$patient)
dds$tumor_descriptor <- as.factor(dds$tumor_descriptor)
design(dds) <- ~ patient + tumor_descriptor
dds <- DESeq(dds)
res <- results(dds)
summary(res)

# Сортируем гены по значению p
resSort <- res[order(res$pvalue),]
head(resSort)

# Аннотируем значимые гены
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
rownames(resSort) <- sub("\\..*$", "", rownames(resSort))
head(resSort, n = 10)
top_genes <- head(rownames(resSort), 30)

geneinfo <- select(org.Hs.eg.db,
                   keys = top_genes,
                   columns = c("ENSEMBL", "SYMBOL", "GENENAME"),
                   keytype = "ENSEMBL")

geneinfo_full <- select(org.Hs.eg.db,
                        keys = rownames(resSort),
                        columns = c("ENSEMBL", "SYMBOL", "GENENAME"),
                        keytype = "ENSEMBL",
                        multiVals = "first")

# Преобразуем ENSEMBL ID в SYMBOL (HUGO)
symbols <- mapIds(org.Hs.eg.db, keys = rownames(resSort), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
symbols <- as.character(symbols)
symbols <- gsub('"', '', symbols)
symbols_clean <- symbols
symbols_clean[is.na(symbols_clean)] <- rownames(resSort)[is.na(symbols_clean)]
rownames(resSort) <- symbols_clean

# Визуализируем результаты значимых генов
rownames(normlzd_dds) <- sub("\\..*$", "", rownames(normlzd_dds))
res_sig <- data.frame(normlzd_dds[geneinfo$ENSEMBL, ])
symbols_30 <- top_genes
symbols_30 <- ifelse(!is.na(geneinfo$SYMBOL), geneinfo$SYMBOL, symbols_30)
rownames(res_sig) <- symbols_30

# Определяем цвета для тепловой карты
heat_colors <- brewer.pal(6, "YlOrRd")

# Создаем объект annotation_col для тепловой карты
annotation_col <- data.frame(tumor_descriptor = sorted_filtered_df$tumor_descriptor)
rownames(annotation_col) <- sorted_filtered_df$barcode
old_row_names <- rownames(annotation_col)
new_row_names <- gsub("-", ".", old_row_names)
rownames(annotation_col) <- new_row_names

# Строим тепловую карту
pheatmap(res_sig,
         color = heat_colors,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         scale = "row",
         annotation_col = annotation_col,
         main = "Heatmap")

# Строим volcano plot
rownames(res) <- sub("\\..*$", "", rownames(res))
volcano_plot <- EnhancedVolcano(resSort,
                                lab = rownames(resSort),
                                x = 'log2FoldChange',
                                y = 'padj',
                                title = 'Volcano Plot',
                                pCutoff = 0.05)
print(volcano_plot)