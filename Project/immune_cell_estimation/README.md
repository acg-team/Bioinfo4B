# MSI-H vs MSS in TCGA-COAD Dataset

## Project Overview

This project investigates the immunogenic differences between microsatellite instability-high (MSI-H) and microsatellite stable (MSS) colorectal cancers (CRC) using data from the TCGA-COAD dataset. Studies have shown that MSI-H CRCs are more immunogenic and have higher tumor-infiltrating lymphocytes (TILs) compared to MSS CRCs. Our objective is to analyze immune cell presence in these groups and validate findings with existing deep learning studies and gene expression data.

## Data Sources

We use two main data sources for this analysis:
1. **Deep Learning Study Results**: We refer to the findings from the article [1] which provides insights into immune cell analysis across various datasets.
2. **TCGA-COAD Gene Expression Data**: We are going to use gene expression data from TCGA-COAD to study immune cell infiltration.

## Tasks

### Task 1: Analyzing immune cells in the TCGA-COAD dataset from slide data

- **Download and Load Data**: Download the slide data and load the clinical information.
- **Filter Data**: Filter the dataset to include only COAD samples.
- **Data Visualization**: Create plots to visualize the distribution of TIL% in MSI-H and MSS samples.

[Task 1 Notebook](https://colab.research.google.com/drive/1ugA9hmgpY9u68tTyInTY4-2uRCFUVVeC?usp=sharing)

### Task 2: Calculating Immunoscore for TCGA-COAD Gene Expression Data

- **Load Gene Expression Data**: Read the gene expression data from the TCGA-COAD dataset.
- **Calculate Immunoscore**: Use the tidyestimate R package [3] to calculate immunoscores based on the gene expression data.
- **Statistical Testing**: Perform statistical tests to determine if any observed differences are statistically significant.

[Task 2 Notebook](https://colab.research.google.com/drive/1L_rQjx3evLs8WPzQdUqgukEkSfmVSBga?usp=sharing)

## References
1. Saltz, Joel, et al. "Spatial organization and molecular correlation of tumor-infiltrating lymphocytes using deep learning on pathology images." Cell reports 23.1 (2018): 181-193.
2. TCGA Data Portal: [https://portal.gdc.cancer.gov/](https://portal.gdc.cancer.gov/)
3. Yoshihara, K., Shahmoradgoli, M., Martínez, E. et al. "Inferring tumour purity and stromal and immune cell admixture from expression data," Nature Communications 4, 2612 (2013)

