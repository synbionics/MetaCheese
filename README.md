## Obiettivo
Eseguire la pipeline end-to-end in modo riproducibile e raccogliere velocemente i risultati chiave (tassonomia, funzione, assemblaggio, MAG).

## Prerequisiti (prima di lanciare)
* Nextflow installato (>=23.x) e Docker attivo in locale/WSL2; sul nodo di cluster usare il profilo dedicato.
* File di input nella cartella `input/<sample>/` con nomi coerenti.
* Spazio su disco per `work/` e `output/`; permessi di scrittura.

## Esecuzione

### Locale/WSL2
```bash
nextflow run main.nf \
  --input_sample <sample> \
  --codice <codice> \
  -resume
```
###  Nodo di cluster ( senza scheduler )
```bash
nextflow run main.nf \
 --input_sample < sample > \
 --codice <codice> \
 -profile hpc - interactive \
 -resume
```

## Controlli
Per fare dei controlli accurati si possono consultare i file `report.html`,
`timeline.html` e `.nextflow.log`, `trace.csv`; in quest’ultimo verificare che
lo stato di nuessun processo sia `status=FAILED`. Controllare che tutte le
cartelle di output esistano.

## Output
`output/<data>_<codice>/04_metaphlan_output/` → profili tassonomici per
campione.

`output/<data>_<codice>/04b_humann_output/` → genefamilies.tsv, pathabundance.tsv
e merge.

`output/<data>_<codice>/15_tormes_output/` → report/annotazioni

## Ripresa
Usare sempre `-resume` per riavviare una run interrotto. Se si desidera un
esperimento “nuovo”, cambiare `–codice` (o la data) così da ottenere una
nuova cartella `output/<data>_<codice>/`.