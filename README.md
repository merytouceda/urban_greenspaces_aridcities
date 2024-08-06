NOTE: :wrench: This repository is under construction

# Soil metagenomics of urban greenspaces in an arid city
Repository for the publication: "Differences in the genomic potential of soil bacterial and phage communities between urban greenspaces and natural arid soils". 


#### :file_folder: This repository
This repository harbors the scripts for analysis and visualization, and final data products from the study of soil bacteria and viruses in urban greenspaces and natural soils around the arid city of Tucson (Arizona, USA).  

Each folder contains scripts and data, and a README explaining its contents. The structure is briefly as follows: 

##### Folders and contents: 
- Bacteria: This folder contains the scripts and data for the analysis and vizualization of bacterial community structure, taxonomy, functional traits, and genetic potential. Bacterial community structure is studied on the Kraken2+Bracken annotations of metagenomic reads. Gene annotation includes the following databases, and thus there will be an R script for the analysis of the annotations with each: 
        - KEGG (used for funtional trait calculation, such as sugar-aminoacid preference)
        - NCyc (nitrogen cycling genes) 
        - CAZY (carbon cycling genes) 
        - CARD (anbtibiotic resistance genes)
        - BACmet (heavy metal resistance genes).
  
- Virus: this folder contains the scripts and data for the analysis and visualization of viral community structure, taxonomy, life history strategies, and AMG (auxiliary metabolic gene) annotation.
  

#### :mount_fuji: Field sampling design
The study included __24 samples__, 2 land uses (urban, natural), 3 vegetation types (grassland, forest, shrubland), 3 replicates per site. Here is the breakdown of the sites: 
   
| landuse        | site              | Vegetation      |
| -------------  | ----------------  |---------------- |
| urban          | Reid Park (RP)    | Grassland       |
| urban          | Himmel Park (HP)  | Grassland       |
| urban          | Rose Canyon (RC)  | Forest          |
| natural        | Sabino Canyon (SC)| Shrubland       |
| natural        | SRER-11B (11B)    | Grassland       |
| natural        | SRER-8 (8)        | Grassland       |
| natural        | SRER-Ex45 (Ex45)  | Grassland       |
| natural        | SRER-UAB (UAB)    | Grassland       |

 *SRER stands for Santa Rita Experimental Range
 Replicates are taken in a transect at 3,6, and 9 meters for each site.

   
#### :bookmark_tabs: Other data locations
Raw data for this project can be found at the [NCBI website](https://www.ncbi.nlm.nih.gov/),  project number: PRJNA1143147. Additionally, final data products can be found as this [Zenodo project](https://zenodo.org/records/13152735).


#### :newspaper: Publication
The publication associated with this work is still under review. 
