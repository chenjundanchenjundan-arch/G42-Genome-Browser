# G42 Genome Browser



## Overview



G42 Genome Browser is a web-based tool for searching and visualizing gene information from the G42 watermelon genome annotation.



The application provides an integrated interface for gene structure visualization, sequence retrieval, functional annotation queries, and BLAST-based sequence similarity searches.



## Features



- Gene ID search

- Gene and transcript information retrieval

- Gene structure visualization, including exon, CDS, UTR, and intron regions

- Transcript sequence and protein sequence retrieval

- Functional annotation search

- Genomic region search

- Batch gene ID queries

- BLASTN, BLASTP, BLASTX, TBLASTN, and TBLASTX searches

- BLAST result annotation and gene identification



## Technology



- Python

- Streamlit

- Pandas

- Biopython

- Plotly

- NCBI BLAST+



## Data



The application uses the G42 genome assembly, gene annotation, CDS sequences, protein sequences, and functional annotation files.



Large reference datasets and BLAST databases are not included in this source-code repository.



## Installation



Install the required Python packages:



```bash

pip install -r requirements.txt

```



## Usage



Run the application:



```bash

streamlit run app.py

```



**Note:** The current application is configured for a local Windows environment. Paths, BLAST executables, reference datasets, and database configuration must be adapted before deployment to a public Linux server.



## License



To be determined.

