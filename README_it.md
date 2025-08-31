# metacheese
Pipeline per analizzare campioni di formaggio a partire da file FASTQ compressi (`.fq.gz`).  
Il programma esegue tutti i passaggi principali dell’analisi metagenomica e salva i risultati in cartelle numerate, una per ogni fase del processo.

---

## Workflow

![Workflow overview](docs/img/workflow.png)


---

## Funzionalità principali
La pipeline esegue automaticamente questi passaggi, in ordine:

- ### Pre-processing 🟠
    **1. Estrazione e rimozione adattatori**  
     *Tool: AdapterRemoval*  
     Estrae i file FASTQ dagli archivi (se presenti), poi rimuove le sequenze di adattatori e di bassa qualità, ottenendo reads “puliti”.
  
    **2. Filtraggio del genoma ospite**  
     *Tool: Bowtie2*  
     Allinea i reads al genoma della specie ospite (es. *Bos taurus*) e scarta quelli che si mappano, lasciando solo le sequenze microbiche.
- ### A) Profilazione del microbioma 🔴
    **3. Profilazione tassonomica**  
       *Tool: MetaPhlAn*  
       Stima la composizione tassonomica del campione, indicando specie presenti e abbondanza relativa.
    
    **4. Profilazione funzionale**  
       *Tool: HUMAnN*  
       Analizza le funzioni biologiche/metaboliche potenzialmente presenti nel microbioma.

- ### B) Assemblaggio e binning 🔵
    **5. Assemblaggio del metagenoma**  
       *Tool: SPAdes*  
       Ricostruisce sequenze contigue (contig) assemblando i reads.
    
    **6. Filtraggio dei contig**  
       *Tool: script custom*  
       Seleziona solo i contig sopra soglie di lunghezza/copertura.
    
    **7. Creazione indici per binning**  
       *Tool: Bowtie2*  
       Costruisce indici per il binning dei contig.
    
    **8. Mappatura e calcolo copertura**  
       *Tool: Bowtie2 + script custom*  
       Mappa nuovamente i reads e calcola la copertura di ciascun contig.
    
    **9. Binning metagenomico**  
       *Tool: MetaBAT2*  
       Raggruppa i contig in MAGs (genomi metagenomici assemblati).
    
    **10. Valutazione qualità MAGs**  
        *Tool: CheckM, CheckM2*  
        Valuta completezza e contaminazione dei MAGs.
    
    **11. Filtraggio e annotazione MAGs di alta qualità**  
        *Tool: script custom + TORMES*  
        Seleziona i MAGs migliori e li annota (geni, pathway, resistenze, ecc).
    
    **12. Analisi finale e report**  
        *Tool: script custom, R*  
        Produce tabelle riassuntive e grafici dei risultati.

---

## Struttura del progetto

    metacheese/
    ├── config/                   # File di configurazione globali della pipeline
    │   └── Config.yml            # Parametri principali della pipeline
    │
    ├── data/                     # Dati di supporto e reference
    │   ├── calculate_diversity.R # Analisi di diversità in R
    │   └── gene/                 # Genomi reference (es. Bos taurus)
    │
    ├── docs/                     # Documentazione tecnica e grafici
    │
    ├── input/                    # Dati di input organizzati per campione
    │   ├── campione_prova1/      # Esempio/test
    │   └── PDO/                  # Dataset reale
    │       ├── D_PR_01A_L2_1.fq.gz
    │       └── ...
    │
    ├── output/                   # Output ordinato per run (data + codice)
    │   ├── 20250728_PDMIDR/
    │   │   ├── 01_AdapterRemoval/
    │   │   ├── 03_bowtie2_output/
    │   │   └── config.yml
    │   └── ...
    │
    ├── scripts/                  # Script e template della pipeline
    │   ├── run_pipeline.sh       # Script master
    │   ├── pipeline/             # Script dei singoli step
    │   ├── templates/            # Template degli script
    │   └── utils/                # Utility (es. build_bowtie2_index.sh, delete.sh)
    │
    ├── Dockerfile                # Ambiente Docker riproducibile
    ├── docker-compose.yml        # Setup avanzato (volumi, container)
    └── README.md                 # Questo file

---

## Requisiti

- **Docker** (consigliata versione ≥ 20.10)  
- **docker-compose** (per gestione container e volumi)  
- Risorse consigliate: **≥64 GB RAM** e CPU multi-core (alcuni step sono molto pesanti).  

Non serve installare tool aggiuntivi: **tutto è già incluso nel container**.

---

## Installazione Docker

### Linux (Ubuntu/Debian)
    sudo apt-get update
    sudo apt-get install -y docker.io docker-compose
    sudo systemctl enable docker --now

### Windows/Mac
Scarica *Docker Desktop* e segui le istruzioni dal sito ufficiale.

### Verifica installazione
    docker --version
    docker compose version
    docker run hello-world

Se compare il messaggio “Hello from Docker!”, l’installazione è OK.

---

## Preparazione dell’ambiente

1) Posizionati nella cartella del progetto (clonata o copiata localmente):
    
        cd metacheese

2) Costruisci l’immagine Docker:
    
        docker build -t metacheese .

3) Verifica che l’immagine esista:
    
        docker images

4) Avvio consigliato con docker-compose (volumi già mappati):
    
        docker compose up -d
        docker ps                  # Container in esecuzione
        docker exec -it metacheese-container bash

   - Le cartelle locali restano sincronizzate con il container.  
   - Puoi monitorare lo stato direttamente in `output/`.

5) (Alternativa) Avvio manuale del container:

        docker run -it --rm           -v $PWD/scripts:/main/scripts           -v $PWD/config:/main/config           -v $PWD/data:/main/data          -v $PWD/input:/main/input           -v $PWD/output:/main/output           metacheese /bin/bash

6) Gestione rapida

   - Ferma il container:

         docker stop metacheese-container

   - Elenca container/immagini e rimozione:

         docker ps -a
         docker images
         docker rm metacheese-container
         docker rmi metacheese:latest

---

## Preparazione del genoma ospite (Bowtie2 index)

Per il filtraggio dell’ospite (es. *Bos taurus*) serve un **indice Bowtie2**.

### Passi

1) **Scarica il genoma** (file FASTA genomico, estensione `.fna/.fa/.fasta`) da NCBI/Ensembl in locale.

2) **Prepara la cartella specie** e **sposta il FASTA** dentro (eseguire comandi dentro la cartella del progetto "metacheese"):

    mkdir -p data/gene/Bos_taurus
    mv /percorso/del/tuo/file/GCA_002263795.4_ARS-UCD2.0_genomic.fna data/gene/Bos_taurus/
    # opzionale (solo per chiarezza): rinomina il file
    mv data/gene/Bos_taurus/GCA_002263795.4_ARS-UCD2.0_genomic.fna data/gene/Bos_taurus/genome.fna

   > Va bene **qualsiasi nome** file, purché sia dentro `data/gene/Bos_taurus/` e con estensione `.fa/.fna/.fasta`.

3) **Costruisci l’indice** (lanciare **dentro al container**):

    bash scripts/build_bowtie2_index.sh

   Quando richiesto, inserisci:

    Bos_taurus

### Risultato atteso

Lo script crea i file indice con **prefisso** `data/gene/Bos_taurus/Bos_taurus`, come si aspetta il `Config.yml`.

    data/gene/Bos_taurus/
    ├── genome.fna                                # (o il tuo .fna/.fa/.fasta)
    ├── Bos_taurus.1.bt2
    ├── Bos_taurus.2.bt2
    ├── Bos_taurus.3.bt2
    ├── Bos_taurus.4.bt2
    ├── Bos_taurus.rev.1.bt2
    └── Bos_taurus.rev.2.bt2
---

## Esempio rapido d’uso (Quickstart)

1) **Prepara i dati di input**  
   Crea una cartella dentro `input/` (es. `input/PDO/`) e inserisci i file `*.fq.gz` da analizzare.

2) **(Opzionale) Configura le risorse**  
   Modifica `config/Config.yml` per impostare thread/RAM e altri parametri di alcuni step.

3) **Avvia la pipeline**

        bash scripts/run_pipeline.sh

   - Scegli `1` (nuova esecuzione)  
   - Inserisci il **nome della cartella input** (es. `PDO`)  
   - Inserisci un **codice descrittivo** (es. `PDMIDR`)  
   - Verrà creata `output/<data>_<codice>` (es. `output/20250728_PDMIDR`).

4) **Durante la run**  
   Gli step vengono eseguiti in ordine e i risultati salvati in sottocartelle dentro `output/<data_codice>/`.

5) **Riprendere una run esistente**

        bash scripts/run_pipeline.sh

   - Scegli `2` (continua esecuzione)  
   - Inserisci il nome completo della cartella output (es. `20250728_PDMIDR`)  
   - Indica lo step da cui ripartire (es. `05` o `04b-last`).

---

## Output e risultati

Alla fine di ogni esecuzione trovi una nuova cartella in `output/` con nome `data_codice` (es. `20250728_PDMIDR`).  
Ogni fase ha la propria sottocartella con i file generati da quello step.

Struttura tipica:

    output/20250728_PDMIDR/
    ├── 01_AdapterRemoval/       # Reads filtrati e puliti dagli adattatori
    ├── 03_bowtie2_output/       # Reads microbici + summary mapping
    ├── 04_metaphlan_output/     # Tabelle di abbondanza tassonomica
    ├── 04b_humann_output/       # Profilazione funzionale (HUMAnN)
    ├── 05_spades_output/        # Contig assemblati (.fasta)
    ├── 06_contig_filter/        # Contig filtrati per qualità/lunghezza
    ├── 07_Bowtie_Index/         # Indici Bowtie2 per binning
    ├── 08_mapping_coverage/     # Coverage dei contig
    ├── 10_metabat_depth/        # Profondità per binning
    ├── 11_metabat_MAG/          # MAGs (FASTA)
    ├── 12_checkm/               # Report qualità (CheckM)
    ├── 13_checkm2/              # Report qualità (CheckM2)
    ├── 14_MAGs_high_quality/    # MAGs selezionati + metadata
    ├── 15_tormes_MAGs/          # Annotazioni finali (TORMES)
    └── config.yml               # Parametri usati per la run

> Suggerimento: controlla che ogni cartella contenga file aggiornati; cartelle vuote possono indicare errori in step precedenti.

---

## Pulizia selettiva dell’output

Per fare spazio o ripetere alcuni step senza perdere tutto, usa:

    bash scripts/clean_output_folders.sh

**Cosa fa:** chiede quale cartella `output/<data_codice>` ripulire e, per ciascuna sottocartella definita, elimina i contenuti **preservando** eventuali file/sottocartelle indicati in una lista.

### Configurazione (in testa allo script)

- `TO_DELETE`: elenco delle sottocartelle da processare (relative a `output/<data_codice>/`).  
  Puoi commentare quelle che non vuoi toccare.

- `PRESERVE`: mappa (associative array) delle eccezioni da tenere per ciascuna cartella.  
  Valori separati da virgola, supporta **wildcard** e **sottocartelle**:
  - Esempio preconfigurato:
    - `PRESERVE["04_metaphlan_output"]="diversity,merged_abundance_table.txt,*.txt"`
    - `PRESERVE["05_spades_output"]="contigs/filtered"`

Se il valore è vuoto (`""`), **non viene preservato nulla** e l’intera cartella viene eliminata.

### Esempio tipico

1. Lancia lo script:
    
        bash scripts/clean_output_folders.sh

2. Scegli la cartella di output (es. `20250728_PDMIDR`) e conferma.  
   Lo script eliminerà i contenuti delle cartelle elencate in `TO_DELETE`, tenendo quanto definito in `PRESERVE`.

> **Attenzione:** le eliminazioni sono **definitive**. Fai un backup se devi conservare i risultati.

---

## Troubleshooting essenziale

- **Docker non parte / permessi**  
  Aggiungi l’utente al gruppo docker e riavvia la sessione:  
      sudo usermod -aG docker $USER

- **`yq` non trovato**  
  All’interno del container è già installato. Se lanci `run_pipeline.sh` fuori dal container, potresti non averlo nel PATH.

- **Risorse insufficienti (RAM/CPU)**  
  Riduci la dimensione del dataset o aumenta le risorse. Alcuni step (SPAdes, MetaBAT) sono particolarmente pesanti.

- **Indice Bowtie2 mancante**  
  Esegui prima:  
      bash scripts/build_bowtie2_index.sh

---

## Citazioni

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
   Narciso M. Quijada, David Rodríguez-Lázaro, Jose María Eiros e Marta Hernández (2019). TORMES: una pipeline automatizzata per l'analisi dell'intero genoma batterico. Bioinformatica , 35(21), 4207–4212, https://doi.org/10.1093/bioinformatics/btz220

## Crediti e licenza

Autore/i: Dorin / Peraz  

