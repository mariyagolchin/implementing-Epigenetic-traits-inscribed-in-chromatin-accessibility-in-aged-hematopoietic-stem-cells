
setwd("/home/golchinpour/projects/0_scATAC_GSE162662")

library(Signac)
library(Seurat)
library(hdf5r)

# peak_bc_matrix.h5 ===> chromosome -cell matrix
counts <- Read10X_h5(filename = "Data/GSM5723631_Young_HSC_filtered_peak_bc_matrix.h5")
# > head(counts[1:5,1:3])
# 5 x 3 sparse Matrix of class "dgCMatrix"
#                      AAACGAAAGTGATTAG-1 AAACGAAAGTTCAACC-1 AAACGAACAAATTGAG-1
# chr1:3670444-3672753                  .                  .                  .
# chr1:3993282-3994009                  .                  .                  .
# chr1:4173610-4174310                  .                  .                  .
# chr1:4228126-4228612                  .                  .                  .
# chr1:4234171-4234350                  .                  .                  .
# > 

meta <- read.csv(
file = 'Data/GSM5723631_Young_HSC_singlecell.csv.gz',
header = TRUE,
row.names = 1)
head(meta)
#                     total duplicate chimeric unmapped lowmapq mitochondrial
# NO_BARCODE         776900    178773     9289    59731   54195         14043
# AAACGAAAGAAAGCAG-1      1         1        0        0       0             0
# AAACGAAAGAAATACC-1      1         0        0        0       0             0
# AAACGAAAGAACCATA-1      4         0        0        0       0             0
# AAACGAAAGAACCCGA-1      1         0        0        0       0             0
# AAACGAAAGAACGACC-1      1         0        0        1       0             0
#                    passed_filters cell_id is__cell_barcode TSS_fragments
# NO_BARCODE                 460869    None                0             0
# AAACGAAAGAAAGCAG-1              0    None                0             0
# AAACGAAAGAAATACC-1              1    None                0             1
# AAACGAAAGAACCATA-1              4    None                0             1
# AAACGAAAGAACCCGA-1              1    None                0             1
# AAACGAAAGAACGACC-1              0    None                0             0
#                    DNase_sensitive_region_fragments enhancer_region_fragments
# NO_BARCODE                                        0                         0
# AAACGAAAGAAAGCAG-1                                0                         0
# AAACGAAAGAAATACC-1                                1                         0
# AAACGAAAGAACCATA-1                                2                         1
# AAACGAAAGAACCCGA-1                                0                         0
# AAACGAAAGAACGACC-1                                0                         0
#                    promoter_region_fragments on_target_fragments
# NO_BARCODE                                 0                   0
# AAACGAAAGAAAGCAG-1                         0                   0
# AAACGAAAGAAATACC-1                         1                   1
# AAACGAAAGAACCATA-1                         2                   3
# AAACGAAAGAACCCGA-1                         1                   1
# AAACGAAAGAACGACC-1                         0                   0
#                    blacklist_region_fragments peak_region_fragments
# NO_BARCODE                                  0                     0
# AAACGAAAGAAAGCAG-1                          0                     0
# AAACGAAAGAAATACC-1                          0                     1
# AAACGAAAGAACCATA-1                          0                     3
# AAACGAAAGAACCCGA-1                          0                     1
# AAACGAAAGAACGACC-1                          0                     0
#                    peak_region_cutsites
# NO_BARCODE                            0
# AAACGAAAGAAAGCAG-1                    0
# AAACGAAAGAAATACC-1                    2
# AAACGAAAGAACCATA-1                    6
# AAACGAAAGAACCCGA-1                    2
# AAACGAAAGAACGACC-1                    0

setwd("/home/golchinpour/projects/0_scATAC_GSE162662/Data")

  chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = './GSM5723631_Young_HSC_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

> chrom_assay

# ChromatinAssay data with 165871 features for 5028 cells
# Variable features: 0 
# Genome: mm10 
# Annotation present: FALSE 
# Motifs present: FALSE 
# Fragment files: 1 
# > 
chrom_assay$counts  chrom_assay$data   


> head(chrom_assay$data[1:6,1:3]) 
6 x 3 sparse Matrix of class "dgCMatrix"
                     AAACGAAAGTGATTAG-1 AAACGAAAGTTCAACC-1 AAACGAACAAATTGAG-1
chr1-3670444-3672753                  .                  .                  .
chr1-3993282-3994009                  .                  .                  .
chr1-4173610-4174310                  .                  .                  .
chr1-4228126-4228612                  .                  .                  .
chr1-4234171-4234350                  .                  .                  .
chr1-4332385-4332895                  .                  .                  .
> 

> head(chrom_assay$counts[1:4,1:3]) 
4 x 3 sparse Matrix of class "dgCMatrix"
                     AAACGAAAGTGATTAG-1 AAACGAAAGTTCAACC-1 AAACGAACAAATTGAG-1
chr1-3670444-3672753                  .                  .                  .
chr1-3993282-3994009                  .                  .                  .
chr1-4173610-4174310                  .                  .                  .
chr1-4228126-4228612                  .                  .                  .
> 


data <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = meta
)

> data
An object of class Seurat 
165871 features across 5028 samples within 1 assay 
Active assay: peaks (165871 features, 0 variable features)
 2 layers present: counts, data
> 

> head(data)
                      orig.ident nCount_peaks nFeature_peaks total duplicate
AAACGAAAGTGATTAG-1 SeuratProject        20006           8700 17314      2988
AAACGAAAGTTCAACC-1 SeuratProject         8290           3904  7086      1064
AAACGAACAAATTGAG-1 SeuratProject        15697           6854 12583      1790
AAACGAACAAGCAGGT-1 SeuratProject        20326           8886 17930      2668
AAACGAACAAGCCCTG-1 SeuratProject        14083           6324 11334      1859
AAACGAACACACACAT-1 SeuratProject        15845           6794 12901      2368
AAACGAACACTCCTCA-1 SeuratProject        21027           8831 17283      2372
AAACGAAGTTTCCGGG-1 SeuratProject        13747           6119 10569      1083
AAACTCGAGCGAGCTA-1 SeuratProject        18193           7725 15056      2679
AAACTCGAGCGTCAAG-1 SeuratProject        12249           5506  9661      1170
                   chimeric unmapped lowmapq mitochondrial passed_filters
AAACGAAAGTGATTAG-1      234      125    1229           445          12293
AAACGAAAGTTCAACC-1      102       79     484           182           5175
AAACGAACAAATTGAG-1      152      152     842           362           9285
AAACGAACAAGCAGGT-1      223      140    1351           204          13344
AAACGAACAAGCCCTG-1      158      105     649           153           8410
AAACGAACACACACAT-1      166      198     732           176           9261
AAACGAACACTCCTCA-1      274      177    1238           480          12742
AAACGAAGTTTCCGGG-1      164      109     714           198           8301
AAACTCGAGCGAGCTA-1      203      227     916           193          10838
AAACTCGAGCGTCAAG-1      121      103     653           146           7468
                   cell_id is__cell_barcode TSS_fragments
AAACGAAAGTGATTAG-1 _cell_0                1          4868
AAACGAAAGTTCAACC-1 _cell_1                1          2030
AAACGAACAAATTGAG-1 _cell_2                1          4015
AAACGAACAAGCAGGT-1 _cell_3                1          5003
AAACGAACAAGCCCTG-1 _cell_4                1          3532
AAACGAACACACACAT-1 _cell_5                1          3892
AAACGAACACTCCTCA-1 _cell_6                1          5467
AAACGAAGTTTCCGGG-1 _cell_7                1          3609
AAACTCGAGCGAGCTA-1 _cell_8                1          4780
AAACTCGAGCGTCAAG-1 _cell_9                1          3116
                   DNase_sensitive_region_fragments enhancer_region_fragments
AAACGAAAGTGATTAG-1                             7782                      4226
AAACGAAAGTTCAACC-1                             3298                      1745
AAACGAACAAATTGAG-1                             6208                      3077
AAACGAACAAGCAGGT-1                             7986                      4327
AAACGAACAAGCCCTG-1                             5500                      2855
AAACGAACACACACAT-1                             6020                      3163
AAACGAACACTCCTCA-1                             8377                      4208
AAACGAAGTTTCCGGG-1                             5527                      2726
AAACTCGAGCGAGCTA-1                             7300                      3654
AAACTCGAGCGTCAAG-1                             4843                      2502
                   promoter_region_fragments on_target_fragments
AAACGAAAGTGATTAG-1                      4330               10002
AAACGAAAGTTCAACC-1                      1815                4201
AAACGAACAAATTGAG-1                      3664                7793
AAACGAACAAGCAGGT-1                      4354               10416
AAACGAACAAGCCCTG-1                      3263                7015
AAACGAACACACACAT-1                      3544                7827
AAACGAACACTCCTCA-1                      4939               10571
AAACGAAGTTTCCGGG-1                      3326                6912
AAACTCGAGCGAGCTA-1                      4363                9201
AAACTCGAGCGTCAAG-1                      2808                6154
                   blacklist_region_fragments peak_region_fragments
AAACGAAAGTGATTAG-1                         14                 10201
AAACGAAAGTTCAACC-1                         10                  4227
AAACGAACAAATTGAG-1                         45                  7967
AAACGAACAAGCAGGT-1                         30                 10381
AAACGAACAAGCCCTG-1                         21                  7159
AAACGAACACACACAT-1                         44                  8040
AAACGAACACTCCTCA-1                         61                 10698
AAACGAAGTTTCCGGG-1                         29                  6991
AAACTCGAGCGAGCTA-1                         70                  9240
AAACTCGAGCGTCAAG-1                         20                  6234
                   peak_region_cutsites
AAACGAAAGTGATTAG-1                20030
AAACGAAAGTTCAACC-1                 8298
AAACGAACAAATTGAG-1                15707
AAACGAACAAGCAGGT-1                20359
AAACGAACAAGCCCTG-1                14103
AAACGAACACACACAT-1                15856
AAACGAACACTCCTCA-1                21051
AAACGAAGTTTCCGGG-1                13762
AAACTCGAGCGAGCTA-1                18211
AAACTCGAGCGTCAAG-1                12271
> 
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)

> head(annotations)
# GRanges object with 6 ranges and 5 metadata columns:
#                      seqnames          ranges strand |              tx_id
#                         <Rle>       <IRanges>  <Rle> |        <character>
#   ENSMUSE00001236884     chr3 3508030-3508332      + | ENSMUST00000108393
#   ENSMUSE00000676606     chr3 3634150-3634347      + | ENSMUST00000108394
#   ENSMUSE00001345708     chr3 3638059-3638230      + | ENSMUST00000108393
#   ENSMUSE00001345708     chr3 3638059-3638230      + | ENSMUST00000108394
#   ENSMUSE00000149313     chr3 3641223-3641317      + | ENSMUST00000108393
#   ENSMUSE00000149313     chr3 3641223-3641317      + | ENSMUST00000108394
#                        gene_name            gene_id   gene_biotype     type
#                      <character>        <character>    <character> <factor>
#   ENSMUSE00001236884       Hnf4g ENSMUSG00000017688 protein_coding     exon
#   ENSMUSE00000676606       Hnf4g ENSMUSG00000017688 protein_coding     exon
#   ENSMUSE00001345708       Hnf4g ENSMUSG00000017688 protein_coding     exon
#   ENSMUSE00001345708       Hnf4g ENSMUSG00000017688 protein_coding     exon
#   ENSMUSE00000149313       Hnf4g ENSMUSG00000017688 protein_coding     exon
#   ENSMUSE00000149313       Hnf4g ENSMUSG00000017688 protein_coding     exon
#   -------
#   seqinfo: 22 sequences (1 circular) from mm10 genome

Annotation(data) <- annotations

data <- NucleosomeSignal(object = data) #fragment ratio 147-294: <147 
# طول میانگین برای چرخیدن دور هیستون147 
data <- TSSEnrichment(object = data, fast = FALSE)
# بر اساس محل ژن - ترنزکریپشن استارت سایت رو نشون میده

data$blacklist_ratio <- data$blacklist_region_fragments / data$peak_region_fragments

data$pct_reads_in_peaks <- data$peak_region_fragments / data$passed_filters * 100 

vplot1 <- VlnPlot(
  object = data,
  features = c('peak_region_fragments', 'pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'TSS.enrichment'),
  pt.size = 0.1,
  ncol = 5
)
vplot <- vplot1 + 
  theme(
    axis.title = element_text(size = 8),  # Adjust the size for axis titles
    plot.title = element_text(size = 10)  # Adjust the size for the plot title
  )

setwd("/home/golchinpour/projects/0_scATAC_GSE162662")
library(ggplot2)
ggsave("Plots/vplot.png", plot = vplot1, width = 10, height = 6, dpi = 300,bg = "white")


############## Filter outliers #####################
low_prf <- quantile(data[["peak_region_fragments"]]$peak_region_fragments, probs = 0.02)
hig_prf <- quantile(data[["peak_region_fragments"]]$peak_region_fragments, probs = 0.98)
low_prp <- quantile(data[["pct_reads_in_peaks"]]$pct_reads_in_peaks, probs = 0.02)

high_blr <- quantile(data[["blacklist_ratio"]]$blacklist_ratio, probs = 0.98)

hig_ns <- quantile(data[["nucleosome_signal"]]$nucleosome_signal, probs = 0.98)

low_ts <- quantile(data[["TSS.enrichment"]]$TSS.enrichment, probs = 0.02)

print(low_prf)
print(hig_prf)
print(low_prp)
print(high_blr)
print(hig_ns)
print(low_ts)


data <- subset(
  x = data,
  subset = peak_region_fragments > low_prf &
    peak_region_fragments < hig_prf &
    pct_reads_in_peaks > low_prp &
    blacklist_ratio < high_blr &
    nucleosome_signal < hig_ns &
    TSS.enrichment > low_ts
)

> data
An object of class Seurat 
165871 features across 4492 samples within 1 assay 
Active assay: peaks (165871 features, 0 variable features)
 2 layers present: counts, data
> 
دیتای قبل
> data
An object of class Seurat 
165871 features across 5028 samples

# پس یه تعدادی سلول رو فیلتر کرده

# نرمالایز کردن دیتا 
data <- RunTFIDF(data)

data <- FindTopFeatures(data, min.cutoff = 'q0')

data <- RunSVD(data) # analogus to PCA

plotc<- DepthCor(data)
ggsave("Plots/DepthCor.png", plot = plotc, width = 10, height = 6, dpi = 300,bg = "white")


data <- RunUMAP(object = data, reduction = 'lsi', dims = 2:30)
data <- FindNeighbors(object = data, reduction = 'lsi', dims = 2:30)
data <- FindClusters(object = data, verbose = FALSE, algorithm = 3)
dimplot<- DimPlot(object = data, label = TRUE)
ggsave("Plots/dimplot.png", plot = dimplot, width = 10, height = 6, dpi = 300,bg = "white")


############# import both samples ## Multiple samples....########################
function to get both file names and all code
#################################################################################
import_atac <- function(count_path, meta_path, fragment_path){
  counts <- Read10X_h5(filename = count_path)
  
  meta <- read.csv(
  file = meta_path,
  header = TRUE,
  row.names = 1)
  
  
  
    chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = fragment_path,
    min.cells = 10,
    min.features = 200
  )
  
  data <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = meta
  )
  
  Annotation(data) <- annotations
  
  
  data <- NucleosomeSignal(object = data) #fragment ratio 147-294: <147  ---  mononucleosome:nucleosome-free
  
  
  data <- TSSEnrichment(object = data, fast = FALSE)
  
  data$blacklist_ratio <- data$blacklist_region_fragments / data$peak_region_fragments
  
  data$pct_reads_in_peaks <- data$peak_region_fragments / data$passed_filters * 100 
  
  low_prf <- quantile(data[["peak_region_fragments"]]$peak_region_fragments, probs = 0.02)
  hig_prf <- quantile(data[["peak_region_fragments"]]$peak_region_fragments, probs = 0.98)
  low_prp <- quantile(data[["pct_reads_in_peaks"]]$pct_reads_in_peaks, probs = 0.02)
  
  high_blr <- quantile(data[["blacklist_ratio"]]$blacklist_ratio, probs = 0.98)
  
  hig_ns <- quantile(data[["nucleosome_signal"]]$nucleosome_signal, probs = 0.98)
  
  low_ts <- quantile(data[["TSS.enrichment"]]$TSS.enrichment, probs = 0.02)
  
  data <- subset(
    x = data,
    subset = peak_region_fragments > low_prf &
      peak_region_fragments < hig_prf &
      pct_reads_in_peaks > low_prp &
      blacklist_ratio < high_blr &
      nucleosome_signal < hig_ns &
      TSS.enrichment > low_ts
  )
  
  
  
  #data <- RunTFIDF(data)
  #data <- FindTopFeatures(data, min.cutoff = 'q0')
  #data <- RunSVD(data)

  return(data)

}

  ##############################################################
  ##############################################################
  setwd("/home/golchinpour/projects/0_scATAC_GSE162662/Data")
  young <- import_atac("GSM5723631_Young_HSC_filtered_peak_bc_matrix.h5",
         'GSM5723631_Young_HSC_singlecell.csv.gz',
         './GSM5723631_Young_HSC_fragments.tsv.gz')

old <- import_atac("GSM5723632_Aged_HSC_filtered_peak_bc_matrix.h5",
         'GSM5723632_Aged_HSC_singlecell.csv.gz',
         './GSM5723632_Aged_HSC_fragments.tsv.gz')


young$dataset <- "young"
old$dataset <- "old"


data <- merge(young, old)

data <- FindTopFeatures(data, min.cutoff = 'q0')
data <- RunTFIDF(data)
data <- RunSVD(data)
data

data <- RunUMAP(object = data, reduction = 'lsi', dims = 2:30)
data <- FindNeighbors(object = data, reduction = 'lsi', dims = 2:30)

data <- FindClusters(object = data, verbose = FALSE, algorithm = 3, resolution = .4)


umap<- DimPlot(object = data, label = TRUE) 

ggsave("/home/golchinpour/projects/0_scATAC_GSE162662/Plots/umap.png",
  plot = umap, 
  width = 8, 
  height = 6, 
  dpi = 300
)

umap1<- DimPlot(object = data, label = TRUE, group.by = "dataset") + NoLegend()

ggsave("/home/golchinpour/projects/0_scATAC_GSE162662/Plots/umap1.png",
  plot = umap1, 
  width = 8, 
  height = 6, 
  dpi = 300
)

gene.activities <- GeneActivity(data)

> head(gene.activities[1:4,1:3])
# 4 x 3 sparse Matrix of class "dgCMatrix"
#       AAACGAAAGTGATTAG-1_1 AAACGAAAGTTCAACC-1_1 AAACGAACAAATTGAG-1_1
# Hnf4g                    .                    .                    .
# Zfhx4                    1                    .                    .
# Pex2                     .                    .                    .
# UBC                      .                    .                    .
# > 

data[['RNA']] <- CreateAssayObject(counts = gene.activities)

data <- NormalizeData(
  object = data,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(data$nCount_RNA)
)

> head(data)
                        orig.ident nCount_peaks nFeature_peaks total duplicate
AAACGAAAGTGATTAG-1_1 SeuratProject        20006           8678 17314      2988
AAACGAAAGTTCAACC-1_1 SeuratProject         8290           3899  7086      1064
AAACGAACAAATTGAG-1_1 SeuratProject        15697           6830 12583      1790
AAACGAACAAGCAGGT-1_1 SeuratProject        20326           8855 17930      2668
AAACGAACAAGCCCTG-1_1 SeuratProject        14083           6314 11334      1859
AAACGAACACACACAT-1_1 SeuratProject        15845           6767 12901      2368
AAACGAACACTCCTCA-1_1 SeuratProject        21027           8811 17283      2372
AAACGAAGTTTCCGGG-1_1 SeuratProject        13747           6106 10569      1083
AAACTCGAGCGAGCTA-1_1 SeuratProject        18193           7707 15056      2679
AAACTCGAGCGTCAAG-1_1 SeuratProject        12249           5495  9661      1170
                     chimeric unmapped lowmapq mitochondrial passed_filters
AAACGAAAGTGATTAG-1_1      234      125    1229           445          12293
AAACGAAAGTTCAACC-1_1      102       79     484           182           5175
AAACGAACAAATTGAG-1_1      152      152     842           362           9285
AAACGAACAAGCAGGT-1_1      223      140    1351           204          13344
AAACGAACAAGCCCTG-1_1      158      105     649           153           8410
AAACGAACACACACAT-1_1      166      198     732           176           9261
AAACGAACACTCCTCA-1_1      274      177    1238           480          12742
AAACGAAGTTTCCGGG-1_1      164      109     714           198           8301
AAACTCGAGCGAGCTA-1_1      203      227     916           193          10838
AAACTCGAGCGTCAAG-1_1      121      103     653           146           7468
                     cell_id is__cell_barcode TSS_fragments
AAACGAAAGTGATTAG-1_1 _cell_0                1          4868
AAACGAAAGTTCAACC-1_1 _cell_1                1          2030
AAACGAACAAATTGAG-1_1 _cell_2                1          4015
AAACGAACAAGCAGGT-1_1 _cell_3                1          5003
AAACGAACAAGCCCTG-1_1 _cell_4                1          3532
AAACGAACACACACAT-1_1 _cell_5                1          3892
AAACGAACACTCCTCA-1_1 _cell_6                1          5467
AAACGAAGTTTCCGGG-1_1 _cell_7                1          3609
AAACTCGAGCGAGCTA-1_1 _cell_8                1          4780
AAACTCGAGCGTCAAG-1_1 _cell_9                1          3116
                     DNase_sensitive_region_fragments enhancer_region_fragments
AAACGAAAGTGATTAG-1_1                             7782                      4226
AAACGAAAGTTCAACC-1_1                             3298                      1745
AAACGAACAAATTGAG-1_1                             6208                      3077
AAACGAACAAGCAGGT-1_1                             7986                      4327
AAACGAACAAGCCCTG-1_1                             5500                      2855
AAACGAACACACACAT-1_1                             6020                      3163
AAACGAACACTCCTCA-1_1                             8377                      4208
AAACGAAGTTTCCGGG-1_1                             5527                      2726
AAACTCGAGCGAGCTA-1_1                             7300                      3654
AAACTCGAGCGTCAAG-1_1                             4843                      2502
                     promoter_region_fragments on_target_fragments
AAACGAAAGTGATTAG-1_1                      4330               10002
AAACGAAAGTTCAACC-1_1                      1815                4201
AAACGAACAAATTGAG-1_1                      3664                7793
AAACGAACAAGCAGGT-1_1                      4354               10416
AAACGAACAAGCCCTG-1_1                      3263                7015
AAACGAACACACACAT-1_1                      3544                7827
AAACGAACACTCCTCA-1_1                      4939               10571
AAACGAAGTTTCCGGG-1_1                      3326                6912
AAACTCGAGCGAGCTA-1_1                      4363                9201
AAACTCGAGCGTCAAG-1_1                      2808                6154
                     blacklist_region_fragments peak_region_fragments
AAACGAAAGTGATTAG-1_1                         14                 10201
AAACGAAAGTTCAACC-1_1                         10                  4227
AAACGAACAAATTGAG-1_1                         45                  7967
AAACGAACAAGCAGGT-1_1                         30                 10381
AAACGAACAAGCCCTG-1_1                         21                  7159
AAACGAACACACACAT-1_1                         44                  8040
AAACGAACACTCCTCA-1_1                         61                 10698
AAACGAAGTTTCCGGG-1_1                         29                  6991
AAACTCGAGCGAGCTA-1_1                         70                  9240
AAACTCGAGCGTCAAG-1_1                         20                  6234
                     peak_region_cutsites nucleosome_signal
AAACGAAAGTGATTAG-1_1                20030         0.5388420
AAACGAAAGTTCAACC-1_1                 8298         0.5071672
AAACGAACAAATTGAG-1_1                15707         0.4554529
AAACGAACAAGCAGGT-1_1                20359         0.5837989
AAACGAACAAGCCCTG-1_1                14103         0.4113751
AAACGAACACACACAT-1_1                15856         0.4168177
AAACGAACACTCCTCA-1_1                21051         0.4569864
AAACGAAGTTTCCGGG-1_1                13762         0.5380245
AAACTCGAGCGAGCTA-1_1                18211         0.4950305
AAACTCGAGCGTCAAG-1_1                12271         0.3905352
                     nucleosome_percentile TSS.enrichment TSS.percentile
AAACGAAAGTGATTAG-1_1                  0.61       5.768361           0.73
AAACGAAAGTTCAACC-1_1                  0.44       5.725227           0.71
AAACGAACAAATTGAG-1_1                  0.19       6.144159           0.88
AAACGAACAAGCAGGT-1_1                  0.78       5.340361           0.49
AAACGAACAAGCCCTG-1_1                  0.07       5.251553           0.44
AAACGAACACACACAT-1_1                  0.08       4.322344           0.05
AAACGAACACTCCTCA-1_1                  0.20       5.966328           0.82
AAACGAAGTTTCCGGG-1_1                  0.60       4.943969           0.25
AAACTCGAGCGAGCTA-1_1                  0.38       5.817879           0.75
AAACTCGAGCGTCAAG-1_1                  0.03       5.665349           0.67
                     blacklist_ratio pct_reads_in_peaks dataset
AAACGAAAGTGATTAG-1_1     0.001372414           82.98218   young
AAACGAAAGTTCAACC-1_1     0.002365744           81.68116   young
AAACGAACAAATTGAG-1_1     0.005648299           85.80506   young
AAACGAACAAGCAGGT-1_1     0.002889895           77.79526   young
AAACGAACAAGCCCTG-1_1     0.002933371           85.12485   young
AAACGAACACACACAT-1_1     0.005472637           86.81568   young
AAACGAACACTCCTCA-1_1     0.005702000           83.95856   young
AAACGAAGTTTCCGGG-1_1     0.004148191           84.21877   young
AAACTCGAGCGAGCTA-1_1     0.007575758           85.25558   young
AAACTCGAGCGTCAAG-1_1     0.003208213           83.47616   young
                     peaks_snn_res.0.4 seurat_clusters nCount_RNA nFeature_RNA
AAACGAAAGTGATTAG-1_1                 0               0       8806         5907
AAACGAAAGTTCAACC-1_1                 0               0       3681         3057
AAACGAACAAATTGAG-1_1                 0               0       6994         5028
AAACGAACAAGCAGGT-1_1                 0               0       9389         6247
AAACGAACAAGCCCTG-1_1                 0               0       6219         4648
AAACGAACACACACAT-1_1                 0               0       6955         4861
AAACGAACACTCCTCA-1_1                 0               0       9487         6180
AAACGAAGTTTCCGGG-1_1                 0               0       6300         4704
AAACTCGAGCGAGCTA-1_1                 0               0       8254         5632
AAACTCGAGCGTCAAG-1_1                 0               0       5468         4169

data[['RNA']]

fplot<- FeaturePlot(
  object = data,
  features = c('Kit', 'Pecam1', 'Itgam'),
  max.cutoff = 'q95'
)

ggsave("/home/golchinpour/projects/0_scATAC_GSE162662/Plots/FeaturePlot.png",
  plot = fplot, 
  width = 8, 
  height = 6, 
  dpi = 300
)



DefaultAssay(data) <- 'peaks'

da_peaks <- FindMarkers(
  object = data,
  ident.1 = rownames(data[[]][data$dataset == "old",]),
  ident.2 = rownames(data[[]][data$dataset == "young",]),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

> head(da_peaks)
                                 p_val avg_log2FC pct.1 pct.2     p_val_adj
chr2-107231278-107232435 1.094643e-134  11.436486 0.112 0.000 1.941448e-129
chr19-6295925-6310866    7.278839e-133  -0.925082 0.392 0.595 1.290968e-127
chr3-99112422-99115022   5.667198e-128   3.657919 0.165 0.012 1.005129e-122
chr17-58643214-58645092  1.818786e-126   3.793431 0.155 0.010 3.225782e-121
chr3-57509149-57513454   1.996828e-118   3.519500 0.154 0.013 3.541555e-113
chr16-23841670-23845601  3.564753e-109   2.610423 0.194 0.031 6.322410e-104
> 
da_peaks$closest_gene <-ClosestFeature(data, regions = rownames(da_peaks))$gene_name
da_peaks$distance <- ClosestFeature(data, regions = rownames(da_peaks))$distance
da_peaks

covplot<- CoveragePlot(
  object = data,
  region = rownames(da_peaks)[2],
  extend.upstream = 10000,
  extend.downstream = 5000,
  group.by = "dataset"
)

ggsave("/home/golchinpour/projects/0_scATAC_GSE162662/Plots/coverage_plot.png",
  plot = covplot, 
  width = 8, 
  height = 6, 
  dpi = 300
)


plot1 <- VlnPlot(
  object = data,
  features = rownames(da_peaks)[2],
  group.by = "dataset"
)
plot2 <- FeaturePlot(
  object = data,
  features = rownames(da_peaks)[2],
  max.cutoff = 'q95'
)

combineplot<- plot1 | plot2
ggsave("/home/golchinpour/projects/0_scATAC_GSE162662/Plots/combineplot.png",
  plot = combineplot, 
  width = 8, 
  height = 6, 
  dpi = 300
)
