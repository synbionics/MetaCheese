# metacheese
Pipeline to analyze cheese samples starting from compressed FASTQ files (`.fq.gz`).  
The program runs the main steps of a metagenomic analysis and saves results in numbered folders—one for each stage of the process.

---

## Workflow

![Workflow overview](docs/img/workflow.png)


---

## Main features
The pipeline automatically executes these steps, in order:

- ### Pre-processing 🟠
    **1. Extraction and adapter trimming**  
         *Tool: AdapterRemoval*  
         Extracts FASTQ files from archives (if present), then removes adapter and low-quality sequences, producing “clean” reads.
      
     **2. Host genome filtering**  
         *Tool: Bowtie2*  
         Aligns reads to the host species genome (e.g., *Bos taurus*) and discards mapped reads, keeping only microbial sequences.

- ### A) Microbiome profiling 🔴
     **3. Taxonomic profiling**  
        *Tool: MetaPhlAn*  
        Estimates the sample’s taxonomic composition, reporting detected species and relative abundance.
  
     **4. Functional profiling**  
        *Tool: HUMAnN*  
        Analyzes the biological/metabolic functions potentially present in the microbiome.

- ### B) Assembly and binning 🔵
     **5. Metagenome assembly**  
        *Tool: SPAdes*  
        Reconstructs contiguous sequences (contigs) by assembling reads.
  
     **6. Contig filtering**  
        *Tool: custom scripts*  
        Selects contigs above length/coverage thresholds.
    
     **7. Index creation for binning**  
        *Tool: Bowtie2*  
        Builds indices for contig binning.
    
     **8. Read mapping and coverage calculation**  
        *Tool: Bowtie2 + custom scripts*  
        Maps reads back and computes coverage for each contig.
    
     **9. Metagenomic binning**  
        *Tool: MetaBAT2*  
        Groups contigs into MAGs (Metagenome-Assembled Genomes).
    
     **10. MAG quality assessment**  
         *Tool: CheckM, CheckM2*  
         Evaluates completeness and contamination of MAGs.
    
     **11. Filtering and annotation of high-quality MAGs**  
         *Tool: custom scripts + TORMES*  
         Selects the best MAGs and annotates them (genes, pathways, resistance, etc.).
    
     **12. Final analysis and reporting**  
         *Tool: custom scripts, R*  
         Produces summary tables and plots.

---

## Project structure

    metacheese/
    ├── config/                   # Global configuration files
    │   └── Config.yml            # Main pipeline parameters
    │
    ├── data/                     # Support data and references
    │   ├── calculate_diversity.R # Diversity analysis in R
    │   └── gene/                 # Reference genomes (e.g., Bos taurus)
    │
    ├── docs/                     # Technical docs and figures
    │
    ├── input/                    # Input data organized by sample
    │   ├── campione_prova1/      # Example/test
    │   └── PDO/                  # Real dataset
    │       ├── D_PR_01A_L2_1.fq.gz
    │       └── ...
    │
    ├── output/                   # Output organized by run (date + code)
    │   ├── 20250728_PDMIDR/
    │   │   ├── 01_AdapterRemoval/
    │   │   ├── 03_bowtie2_output/
    │   │   └── config.yml
    │   └── ...
    │
    ├── scripts/                  # Pipeline scripts and templates
    │   ├── run_pipeline.sh       # Master script
    │   ├── pipeline/             # Individual step scripts
    │   ├── templates/            # Script templates
    │   └── utils/                # Utilities (e.g., build_bowtie2_index.sh, delete.sh)
    │
    ├── Dockerfile                # Reproducible Docker environment
    ├── docker-compose.yml        # Advanced setup (volumes, container)
    └── README.md                 # This file

---

## Requirements

- **Docker** (recommended version ≥ 20.10)  
- **docker-compose** (for container and volume management)  
- Recommended resources: **≥ 64 GB RAM** and multi-core CPU (some steps are heavy).

No extra tools required: **everything is already included in the container**.

---

## Docker installation

### Linux (Ubuntu/Debian)
    sudo apt-get update
    sudo apt-get install -y docker.io docker-compose
    sudo systemctl enable docker --now

### Windows/Mac
Download *Docker Desktop* and follow the official instructions.

### Verify installation
    docker --version
    docker compose version
    docker run hello-world

If you see “Hello from Docker!”, the installation is OK.

---

## Environment setup

1) Go to the project folder (cloned or copied locally):
    
        cd metacheese

2) Build the Docker image:
    
        docker build -t metacheese .

3) Check the image exists:
    
        docker images

4) Recommended start with docker-compose (volumes already mapped):
    
        docker compose up -d
        docker ps                  # Running containers
        docker exec -it metacheese-container bash

   - Local folders remain synced with the container.
   - You can monitor progress directly under `output/`.

5) (Alternative) Manual container start:

        docker run -it --rm           -v $PWD/scripts:/main/scripts           -v $PWD/config:/main/config           -v $PWD/data:/main/data          -v $PWD/input:/main/input           -v $PWD/output:/main/output           metacheese /bin/bash

6) Quick management

   - Stop the container:

         docker stop metacheese-container

   - List containers/images and remove:

         docker ps -a
         docker images
         docker rm metacheese-container
         docker rmi metacheese:latest

---

## Host genome preparation (Bowtie2 index)

For host filtering (eg. *Bos taurus*) you need **Bowtie2 index**.

### Steps

1) **Download the genome** (genomic FASTA, `.fna/.fa/.fasta`) from NCBI/Ensembl to your machine.

2) **Create the species folder** and **move the FASTA** inside it (run these in the project folder “metacheese”):

    mkdir -p data/gene/Bos_taurus
    mv /path/to/your/file/GCA_002263795.4_ARS-UCD2.0_genomic.fna data/gene/Bos_taurus/
    # optional (for clarity): rename the file
    mv data/gene/Bos_taurus/GCA_002263795.4_ARS-UCD2.0_genomic.fna data/gene/Bos_taurus/genome.fna

   > **Any filename** is fine as long as it’s inside `data/gene/Bos_taurus/` with `.fa/.fna/.fasta`.

3) **Build the index** (run **inside the container**):

    bash scripts/build_bowtie2_index.sh

   When prompted, enter:

    Bos_taurus

### Expected result  

The script creates index files with **prefix** `data/gene/Bos_taurus/Bos_taurus`, as expected by `Config.yml`.

    data/gene/Bos_taurus/
    ├── genome.fna                                # (or your .fna/.fa/.fasta)
    ├── Bos_taurus.1.bt2
    ├── Bos_taurus.2.bt2
    ├── Bos_taurus.3.bt2
    ├── Bos_taurus.4.bt2
    ├── Bos_taurus.rev.1.bt2
    └── Bos_taurus.rev.2.bt2
---

## Quickstart example

1) **Prepare input data**  
   Create a folder under `input/` (e.g., `input/PDO/`) and place the `*.fq.gz` files to analyze.

2) **(Optional) Configure resources**  
   Edit `config/Config.yml` to set threads/RAM and other step parameters.
     FORMAT OF STEP DEFINITIONS: [stepXX]="<TEMPLATE_FILE> <OUTPUT_FILE> <PREFIX> <ph1=YAMLkey1> <ph2=YAMLkey2> ..."
    	- TEMPLATE_FILE  = name of the template script (inside scripts/templates/)
     	- OUTPUT_FILE    = path where the generated script will be written
    	- PREFIX         = prefix used in placeholders (usually the step number, e.g. 03 or 00-01)
     	- ph=YAMLkey     = mapping between:
                        	• ph  = placeholder name used in the template (@<PREFIX>_<ph>@)
                         	• YAMLkey = key name inside config.yml (under section stepXX)
    
     HOW TO CHANGE NAME OR VALUE:
     	- If you only want to change the value → edit config.yml
     	- If you want to rename the placeholder in the template → also change the left part (ph) in STEPS
     	- If you want to rename the key in config.yml → also change the right part (YAMLkey) in STEPS
    
     Example of renaming:
     	  Template:   @03_vartemplate@
     	  STEPS:      vartemplate=newvarconfig
    	  Config.yml: step03.newvarconfig: /new/path
     	→ Result: @03_vartemplate@ becomes /new/path

3) **Run the pipeline**

        bash scripts/run_pipeline.sh

   - Choose `1` (new run)  
   - Enter the **input folder name** (e.g.,s. `PDO`)  
   - Enter a **descriptive code** (e.g., `PDMIDR`)  
   - This will create `output/<date>_<code>` (e.g., `output/20250728_PDMIDR`).

4) **During the run**  
   Steps are executed in order and results are saved in subfolders under `output/<date>_<code>/`.

5) **Resume an existing run**

        bash scripts/run_pipeline.sh

   - Choose `2` (continue run)  
   - Enter the full output folder name (e.g., `20250728_PDMIDR`)  
   - Specify the step to restart from (e.g., `05` or `04b-last`).

---

## Outputs and results

After each run, you’ll find a new folder under `output/` named `data_codice` (e.g., 20250728_PDMIDR).
Each phase has its own subfolder containing the files produced by that step.

Typical structure:

    output/20250728_PDMIDR/
    ├── 01_AdapterRemoval/       # Adapter-trimmed/cleaned reads
    ├── 03_bowtie2_output/       # Microbial reads + mapping summary
    ├── 04_metaphlan_output/     # Taxonomic abundance tables
    ├── 04b_humann_output/       # Functional profiling (HUMAnN)
    ├── 05_spades_output/        # Assembled contigs (.fasta)
    ├── 06_contig_filter/        # Contigs filtered by quality/length
    ├── 07_Bowtie_Index/         # Bowtie2 indices for binning
    ├── 08_mapping_coverage/     # Contig coverage
    ├── 10_metabat_depth/        # Depth for binning
    ├── 11_metabat_MAG/          # MAGs (FASTA)
    ├── 12_checkm/               # Quality reports (CheckM)
    ├── 13_checkm2/              # Quality reports (CheckM2)
    ├── 14_MAGs_high_quality/    # Selected MAGs + metadata
    ├── 15_tormes_MAGs/          # Final annotations (TORMES)
    └── config.yml               # Parameters used for the run

> Tip: ensure each folder contains up-to-date files; empty folders may indicate errors in previous steps.

---

## Selective output cleanup

To free space or rerun specific steps without deleting everything, use:

    bash scripts/clean_output_folders.sh

**What it does:** 
asks which `output/<data_codice>` folder to clean and, for each defined subfolder, deletes contents while **preserving** any files/subfolders listed.

### Configuration (at the top of the script)

- `TO_DELETE`: list of subfolders to process (relative to `output/<data_codice>/`).
Comment out those you don’t want to touch.

- `PRESERVE`: map (associative array) of exceptions to keep for each folder.
Comma-separated values; supports **wildcards** and **subfolders**:
  - Preconfigured example:
    - `PRESERVE["04_metaphlan_output"]="diversity,merged_abundance_table.txt,*.txt"`
    - `PRESERVE["05_spades_output"]="contigs/filtered"`

If the value is empty (`""`), **nothing is preserved** and the entire folder is removed.

### Typical example

1. Run the script:
    
        bash scripts/clean_output_folders.sh

2. Choose the output folder (e.g., `20250728_PDMIDR`) and confirm.  
   The script will delete contents of the folders listed in `TO_DELETE`, keeping whatever is defined in `PRESERVE`.

> **Warning:** deletions are **permanent**. Back up results you need to keep.

---

## Essential troubleshooting

- **Docker won’t start / permissions**  
  Add your user to the docker group and restart the session::  
      sudo usermod -aG docker $USER

- **`yq` not found**  
  It’s already installed inside the container. If you run `run_pipeline.sh` outside the container, you may not have it in your PATH.

- **Insufficient resources (RAM/CPU)**  
  Reduce dataset size or increase resources. Some steps (SPAdes, MetaBAT) are particularly heavy.

- **Missing Bowtie2 index**  
  Run first:  
      bash scripts/build_bowtie2_index.sh

---

## Citation

**AdapterRemoval v2**
   Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 12;9(1):88 <http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2>

**MetaPhlAn**
   https://doi.org/10.1038/s41587-023-01688-w
   Aitor Blanco-Miguez, Francesco Beghini, Fabio Cumbo, Lauren J. McIver, Kelsey N. Thompson, Moreno Zolfo, Paolo Manghi, Leonard Dubois, Kun D. Huang, Andrew Maltez Thomas, Gianmarco Piccinno, Elisa Piperni, Michal Punčochář, Mireia Valles-Colomer, Adrian Tett, Francesca Giordano, Richard Davies, Jonathan Wolf, Sarah E. Berry, Tim D. Spector, Eric A. Franzosa, Edoardo Pasolli, Francesco Asnicar, Curtis Huttenhower, Nicola Segata. Nature Biotechnology (2023)

**HUMAnN**
   Francesco Beghini1 ,Lauren J McIver2 ,Aitor Blanco-Mìguez1 ,Leonard Dubois1 ,Francesco Asnicar1 ,Sagun Maharjan2,3 ,Ana Mailyan2,3 ,Andrew Maltez Thomas1 ,Paolo Manghi1 ,Mireia Valles-Colomer1 ,George Weingart2,3 ,Yancong Zhang2,3 ,Moreno Zolfo1 ,Curtis Huttenhower2,3 ,Eric A Franzosa2,3 ,Nicola Segata1,4
   https://doi.org/10.7554/eLife.65088

   1 Department CIBIO, University of Trento, Italy
   2 Harvard T. H. Chan School of Public Health, Boston, MA, USA
   3 The Broad Institute of MIT and Harvard, Cambridge, MA, USA
   4 IEO, European Institute of Oncology IRCCS, Milan, Italy

**SPAdes**
   https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpbi.102

**SAMtools**
   Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li
   GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008*

**Tormes**
   Narciso M. Quijada, David Rodríguez-Lázaro, Jose María Eiros e Marta Hernández (2019). TORMES: an automated pipeline for whole bacterial genome analysis. Bioinformatica , 35(21), 4207–4212, https://doi.org/10.1093/bioinformatics/btz220

## Credits & License

Author(s): Dorin / Davide  


