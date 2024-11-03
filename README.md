# Unified Omics Pipelines

```mermaid
graph TD
    %% Main branches for each omics type
    Omics[Omics Pipelines] --> Genomics
    Omics --> Transcriptomics
    Omics --> Proteomics
    Omics --> Metabolomics
    Omics --> Epigenomics
    Omics --> Microbiomics
    Omics --> Phenomics
    Omics --> Lipidomics
    Omics --> Glycomics

    %% Genomics Pipeline
    Genomics --> VG[Variant Calling]
    Genomics --> GA[Genome Assembly]
    Genomics --> AN[Annotation]
    VG --> GT[GATK, BWA, SAMtools]
    GA --> GF[SPAdes, Canu, Flye]
    AN --> GP[Prokka, MAKER, Ensembl VEP]

    %% Transcriptomics Pipeline
    Transcriptomics --> RNASeq[RNA-Seq]
    Transcriptomics --> SCRNASeq[Single-Cell RNA-Seq]
    Transcriptomics --> LRSeq[Long-Read RNA Sequencing]
    RNASeq --> RS[STAR, HISAT2, DESeq2, edgeR]
    SCRNASeq --> SC[Cell Ranger, Seurat, Scanpy]
    LRSeq --> LS[PacBio Iso-Seq, FLAIR, TALON]

    %% Proteomics Pipeline
    Proteomics --> MS[Mass Spec Analysis]
    Proteomics --> PPI[Protein Interaction]
    Proteomics --> SP[Structural Proteomics]
    MS --> MT[MaxQuant, Mascot, Perseus]
    PPI --> PP[STRING, Cytoscape, iRefWeb]
    SP --> ST[Phenix, RELION, Rosetta]

    %% Metabolomics Pipeline
    Metabolomics --> MD[MS & NMR Analysis]
    Metabolomics --> PM[Pathway Mapping]
    Metabolomics --> BD[Biomarker Discovery]
    MD --> MDX[XCMS, MetaboAnalyst, MZmine]
    PM --> PT[KEGG Mapper, Pathway Tools]
    BD --> BT[MetaboAnalyst, SIMCA]

    %% Epigenomics Pipeline
    Epigenomics --> DMA[DNA Methylation]
    Epigenomics --> HMA[Histone Mod Analysis]
    Epigenomics --> ATAC[ATAC-Seq]
    DMA --> DMT[Bismark, methylKit, MethylSig]
    HMA --> HM[MACS2, DiffBind]
    ATAC --> ATA[ATACseqQC, MACS2, Homer]

    %% Microbiomics Pipeline
    Microbiomics --> M16S[16S rRNA Sequencing]
    Microbiomics --> MG[Metagenomics]
    Microbiomics --> MT[Metatranscriptomics]
    M16S --> MQ[QIIME2, mothur, DADA2]
    MG --> MH[MetaPhlAn, HUMAnN, Kraken]
    MT --> MSAM[SAMSA, Metatools, HUMAnN2]

    %% Phenomics Pipeline
    Phenomics --> IA[Image Analysis]
    Phenomics --> QTL[QTL Mapping]
    Phenomics --> HTP[High-Throughput Pheno]
    IA --> IP[ImageJ, PlantCV]
    QTL --> QT[R/qtl, MapQTL]
    HTP --> HT[PhenoBox, PhenoDyn]

    %% Lipidomics Pipeline
    Lipidomics --> LI[Lipid Quantification]
    Lipidomics --> LPA[Pathway Analysis]
    Lipidomics --> LBD[Biomarker Discovery]
    LI --> LIT[LipidSearch, MZmine, LipiDex]
    LPA --> LPT[KEGG, LIPID MAPS]
    LBD --> LB[MetaboAnalyst, SIMCA]

    %% Glycomics Pipeline
    Glycomics --> GI[Glycan Quantification]
    Glycomics --> GP[Glycoproteomics]
    Glycomics --> GPA[Glycan Pathway Analysis]
    GI --> GT[GlycoWorkbench, SimGlycan]
    GP --> GM[GlycoMod, Byonic]
    GPA --> GK[KEGG Glycan, GlyGen]
