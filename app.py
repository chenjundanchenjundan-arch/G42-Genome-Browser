
from pathlib import Path
from collections import defaultdict
import gzip
import streamlit as st
import pandas as pd
import plotly.graph_objects as go

from Bio import SeqIO
from Bio.Seq import Seq
from pyfaidx import Fasta
from urllib.parse import unquote


# ============================================================
# 1. Basic configuration
# ============================================================

BASE = Path(r"F:\worker\G42\G42_final\G42")

GFF_FILE = BASE / "G42.gff3.gz"
CDS_FILE = BASE / "G42.cds.fa.gz"
PROTEIN_FILE = BASE / "G42.protein.fa.gz"
GENOME_FILE = BASE / "G42.genome.fasta"
ANNOTATION_FILE = BASE / "G42_functional_annotations.tsv"

st.set_page_config(
    page_title="G42 Genome Browser",
    page_icon="🧬",
    layout="wide"
)


# ============================================================
# 2. Parsing functions
# ============================================================

def parse_attributes(attribute_string):
    """Parse GFF3 attributes."""
    result = {}

    for item in attribute_string.strip().split(";"):
        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            result[unquote(key)] = unquote(value)

    return result


@st.cache_resource(show_spinner="Loading G42 genome data...")
def load_data():

    # --------------------------------------------------------
    # Load genome FASTA
    # --------------------------------------------------------

    genome = Fasta(
        str(GENOME_FILE),
        as_raw=True,
        sequence_always_upper=True
    )

    # --------------------------------------------------------
    # Load CDS FASTA
    # --------------------------------------------------------

    cds_dict = {}

    with gzip.open(CDS_FILE, "rt") as f:
        for record in SeqIO.parse(f, "fasta"):
            cds_dict[record.id] = str(record.seq)

    # --------------------------------------------------------
    # Load protein FASTA
    # --------------------------------------------------------

    protein_dict = {}

    with gzip.open(PROTEIN_FILE, "rt") as f:
        for record in SeqIO.parse(f, "fasta"):
            protein_dict[record.id] = str(record.seq)

    # --------------------------------------------------------
    # Load functional annotations
    # --------------------------------------------------------

    annotation_df = pd.read_csv(
        ANNOTATION_FILE,
        sep="\t",
        dtype=str
    ).fillna("")

    annotation_dict = {}

    for _, row in annotation_df.iterrows():
        annotation_dict[row["ID"]] = row["Functional_Description"]

    # --------------------------------------------------------
    # Parse GFF3
    # --------------------------------------------------------

    genes = {}
    transcripts = {}

    with gzip.open(GFF_FILE, "rt") as f:

        for line in f:

            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            chrom = fields[0]
            feature = fields[2]

            start = int(fields[3])
            end = int(fields[4])

            score = fields[5]
            strand = fields[6]
            phase = fields[7]

            attrs = parse_attributes(fields[8])

            # Gene records
            if feature == "gene":

                gene_id = attrs.get("ID")

                if gene_id:
                    genes[gene_id] = {
                        "gene_id": gene_id,
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "features": []
                    }

            # Transcript records
            elif feature in ("mRNA", "transcript"):

                transcript_id = attrs.get("ID")
                gene_id = attrs.get("Parent")

                if transcript_id and gene_id:

                    transcripts[transcript_id] = {
                        "transcript_id": transcript_id,
                        "gene_id": gene_id,
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "features": []
                    }

            # Child features
            else:

                parent = attrs.get("Parent", "")

                if not parent:
                    continue

                for transcript_id in parent.split(","):

                    if transcript_id in transcripts:

                        transcripts[transcript_id]["features"].append({
                            "type": feature,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "phase": phase,
                            "score": score,
                            "attributes": attrs
                        })

    # --------------------------------------------------------
    # Build gene -> transcript relationship
    # --------------------------------------------------------

    gene_transcripts = defaultdict(list)

    for transcript_id, transcript in transcripts.items():

        gene_id = transcript["gene_id"]

        gene_transcripts[gene_id].append(transcript_id)

        if gene_id not in genes:
            genes[gene_id] = {
                "gene_id": gene_id,
                "chrom": transcript["chrom"],
                "start": transcript["start"],
                "end": transcript["end"],
                "strand": transcript["strand"],
                "features": []
            }

    return (
        genome,
        genes,
        transcripts,
        dict(gene_transcripts),
        cds_dict,
        protein_dict,
        annotation_dict
    )


# ============================================================
# 3. Sequence formatting
# ============================================================

def format_fasta(sequence, sequence_id, width=60):

    lines = [f">{sequence_id}"]

    for i in range(0, len(sequence), width):
        lines.append(sequence[i:i + width])

    return "\n".join(lines)


# ============================================================
# 4. Gene structure visualization
# ============================================================


def build_transcript_comparison(
    transcript_ids,
    transcripts,
    cds_dict,
    protein_dict,
    annotation_dict
):
    """
    Compare all transcripts belonging to one gene.
    """

    rows = []

    for tid in transcript_ids:

        transcript = transcripts[tid]
        features = transcript["features"]

        exon_count = sum(
            1 for x in features if x["type"] == "exon"
        )

        intron_count = sum(
            1 for x in features if x["type"] == "intron"
        )

        cds_length = sum(
            x["end"] - x["start"] + 1
            for x in features
            if x["type"] == "CDS"
        )

        cds_sequence = cds_dict.get(tid, "")
        protein_sequence = protein_dict.get(tid, "")

        annotation = annotation_dict.get(
            tid,
            "Unknown"
        )

        rows.append({
            "Transcript ID": tid,
            "Exon count": exon_count,
            "Intron count": intron_count,
            "CDS length (bp)": cds_length,
            "CDS FASTA length (bp)": len(cds_sequence),
            "Protein length (aa)": len(protein_sequence),
            "Functional annotation": annotation
        })

    return pd.DataFrame(rows)


def draw_chromosome_location(
    genome,
    genes,
    selected_gene_id
):
    """
    Display all genes on the selected chromosome
    and highlight the selected gene.
    """

    selected = genes[selected_gene_id]

    chrom = selected["chrom"]

    chromosome_length = len(genome[chrom])

    x_positions = []
    gene_ids = []
    starts = []
    ends = []

    for gene_id, gene in genes.items():

        if gene["chrom"] != chrom:
            continue

        midpoint = (gene["start"] + gene["end"]) / 2

        x_positions.append(midpoint)
        gene_ids.append(gene_id)
        starts.append(gene["start"])
        ends.append(gene["end"])

    fig = go.Figure()

    # All genes on this chromosome
    fig.add_trace(
        go.Scatter(
            x=x_positions,
            y=[0] * len(x_positions),
            mode="markers",
            marker=dict(
                color="#94A3B8",
                size=5,
                opacity=0.65
            ),
            text=gene_ids,
            customdata=list(zip(starts, ends)),
            hovertemplate=(
                "Gene: %{text}<br>"
                "Start: %{customdata[0]:,}<br>"
                "End: %{customdata[1]:,}"
                "<extra></extra>"
            ),
            name="Other genes"
        )
    )

    # Selected gene
    midpoint = (
        selected["start"] + selected["end"]
    ) / 2

    fig.add_trace(
        go.Scatter(
            x=[midpoint],
            y=[0],
            mode="markers",
            marker=dict(
                color="#DC2626",
                size=15,
                symbol="diamond"
            ),
            text=[selected_gene_id],
            hovertemplate=(
                "Selected gene: %{text}<br>"
                "Start: " + f"{selected['start']:,}" +
                "<br>End: " + f"{selected['end']:,}" +
                "<extra></extra>"
            ),
            name="Selected gene"
        )
    )

    # Chromosome backbone
    fig.add_shape(
        type="line",
        x0=0,
        x1=chromosome_length,
        y0=0,
        y1=0,
        line=dict(
            color="#64748B",
            width=3
        ),
        layer="below"
    )

    # Highlight selected gene interval
    fig.add_shape(
        type="rect",
        x0=selected["start"],
        x1=selected["end"],
        y0=-0.18,
        y1=0.18,
        fillcolor="rgba(220,38,38,0.18)",
        line=dict(
            color="#DC2626",
            width=2
        )
    )

    fig.update_layout(
        title=f"{chrom} — Genome-wide Gene Distribution",
        xaxis_title="Chromosome coordinate (bp)",
        xaxis=dict(
            range=[0, chromosome_length],
            tickformat=",",
            showgrid=True
        ),
        yaxis=dict(
            visible=False,
            range=[-0.5, 0.5]
        ),
        height=250,
        template="plotly_white",
        margin=dict(l=20, r=20, t=50, b=50),
        legend=dict(
            orientation="h",
            y=-0.25
        )
    )

    return fig




def draw_gene_structure(transcript):

    features = transcript["features"]

    start = transcript["start"]
    end = transcript["end"]
    strand = transcript["strand"]

    # Track positions
    tracks = {
        "Exon": 3,
        "CDS": 2,
        "UTR": 1,
        "Intron": 0
    }

    colors = {
        "Exon": "#16A34A",
        "CDS": "#2563EB",
        "UTR": "#F59E0B",
        "Intron": "#64748B"
    }

    fig = go.Figure()

    # --------------------------------------------------------
    # 1. Draw transcript coordinate backbone
    # --------------------------------------------------------

    for name, y in tracks.items():

        fig.add_shape(
            type="line",
            x0=start,
            x1=end,
            y0=y,
            y1=y,
            line=dict(
                color="#CBD5E1",
                width=1,
                dash="dot"
            ),
            layer="below"
        )

    # --------------------------------------------------------
    # 2. Draw feature tracks
    # --------------------------------------------------------

    feature_mapping = {
        "exon": "Exon",
        "CDS": "CDS",
        "five_prime_UTR": "UTR",
        "three_prime_UTR": "UTR",
        "UTR": "UTR",
        "intron": "Intron"
    }

    displayed = set()

    for feature in features:

        ftype = feature["type"]

        if ftype not in feature_mapping:
            continue

        track = feature_mapping[ftype]

        y = tracks[track]

        x0 = feature["start"]
        x1 = feature["end"]

        if track == "Intron":

            fig.add_shape(
                type="line",
                x0=x0,
                x1=x1,
                y0=y,
                y1=y,
                line=dict(
                    color=colors[track],
                    width=4
                )
            )

        else:

            height = 0.28 if track == "CDS" else 0.38

            fig.add_shape(
                type="rect",
                x0=x0,
                x1=x1,
                y0=y - height,
                y1=y + height,
                fillcolor=colors[track],
                line=dict(width=0),
                opacity=0.95
            )

        displayed.add(track)

    # --------------------------------------------------------
    # 3. Add legend
    # --------------------------------------------------------

    for name in ["Exon", "CDS", "UTR", "Intron"]:

        if name not in displayed:
            continue

        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(
                    color=colors[name],
                    size=10,
                    symbol="square"
                ),
                name=name
            )
        )

    # --------------------------------------------------------
    # 4. Add strand direction
    # --------------------------------------------------------

    if strand == "+":

        arrow_x = end
        arrow_text = "5′ → 3′"

    else:

        arrow_x = start
        arrow_text = "3′ ← 5′"

    fig.add_annotation(
        x=arrow_x,
        y=3.65,
        text=arrow_text,
        showarrow=False,
        font=dict(
            size=12,
            color="#DC2626"
        ),
        xanchor="right" if strand == "+" else "left"
    )

    # --------------------------------------------------------
    # 5. Layout
    # --------------------------------------------------------

    fig.update_layout(
        title=dict(
            text=f"Transcript Structure: {transcript['transcript_id']}",
            x=0.02
        ),

        xaxis=dict(
            title="Genomic coordinate (bp)",
            range=[start, end],
            showgrid=True,
            gridcolor="#E5E7EB",
            tickformat=","
        ),

        yaxis=dict(
            tickmode="array",
            tickvals=[3, 2, 1, 0],
            ticktext=["Exon", "CDS", "UTR", "Intron"],
            range=[-0.6, 4.0],
            fixedrange=True
        ),

        height=390,

        template="plotly_white",

        margin=dict(
            l=90,
            r=30,
            t=70,
            b=65
        ),

        legend=dict(
            orientation="h",
            y=-0.25,
            x=0
        ),

        hovermode="x"
    )

    return fig


# ============================================================
# 5. Main website
# ============================================================

st.title("🧬 G42 Genome Browser")

st.markdown(
    "### Watermelon G42 Genome and Gene Information Portal"
)

st.caption(
    "Gene structure | Functional annotation | DNA | CDS | Protein"
)

st.divider()


# Load data
try:

    (
        genome,
        genes,
        transcripts,
        gene_transcripts,
        cds_dict,
        protein_dict,
        annotation_dict
    ) = load_data()

except Exception as e:

    st.error(f"Data loading failed: {e}")
    st.stop()


# Sidebar
st.sidebar.title("G42 Database")

st.sidebar.metric("Genes", f"{len(genes):,}")
st.sidebar.metric("Transcripts", f"{len(transcripts):,}")
st.sidebar.metric("CDS sequences", f"{len(cds_dict):,}")
st.sidebar.metric("Protein sequences", f"{len(protein_dict):,}")

st.sidebar.divider()

st.sidebar.info(
    "Enter a G42 gene ID to retrieve its genomic information."
)


# Search
gene_id_input = st.text_input(
    "Enter Gene ID",
    placeholder="Example: ClG42_01g0000100"
)

search_button = st.button(
    "Search Gene",
    type="primary"
)

if search_button:

    st.session_state["searched_gene"] = gene_id_input.strip()


# --------------------------------------------------------
# Functional annotation search
# --------------------------------------------------------

st.divider()

st.subheader("Search by Functional Annotation")

function_query = st.text_input(
    "Enter functional keyword",
    placeholder="Example: kinase, transporter, transcription factor"
)

if st.button("Search Functional Annotation"):

    if not function_query.strip():

        st.warning("Please enter a keyword.")

    else:

        query = function_query.strip().lower()

        matched_genes = {}

        for transcript_id, description in annotation_dict.items():

            if query in description.lower():

                if transcript_id in transcripts:

                    gid = transcripts[transcript_id]["gene_id"]

                    if gid not in matched_genes:
                        matched_genes[gid] = []

                    matched_genes[gid].append(
                        (transcript_id, description)
                    )

        st.session_state["function_search_results"] = matched_genes
        st.session_state["function_search_keyword"] = function_query.strip()


# Display results
matched_genes = st.session_state.get(
    "function_search_results",
    {}
)

if matched_genes:

    keyword = st.session_state.get(
        "function_search_keyword",
        ""
    )

    st.success(
        f"Found {len(matched_genes)} matching genes for: {keyword}"
    )

    result_rows = []

    for gid, items in matched_genes.items():

        for tid, description in items:

            result_rows.append({
                "Gene ID": gid,
                "Transcript ID": tid,
                "Functional Description": description
            })

    result_df = pd.DataFrame(result_rows)

    st.dataframe(
        result_df,
        width="stretch",
        hide_index=True
    )

    available_gene_ids = sorted(matched_genes.keys())

    selected_function_gene = st.selectbox(
        "Select a gene from search results",
        available_gene_ids,
        key="function_gene_choice"
    )

    if st.button(
        "Open Selected Gene",
        key="open_function_gene"
    ):

        st.session_state["searched_gene"] = selected_function_gene

        st.rerun()

elif st.session_state.get("function_search_keyword"):

    st.warning("No matching genes found.")


# ============================================================
# Genome Region Search
# ============================================================

st.divider()
st.header("Genome Region Search")

st.write(
    "Search all genes overlapping a chromosome interval."
)

chromosome_names = sorted(list(genome.keys()))

with st.form("region_search_form"):

    selected_chrom = st.selectbox(
        "Select chromosome",
        chromosome_names,
        key="region_chrom"
    )

    chromosome_length = len(genome[selected_chrom])

    col1, col2 = st.columns(2)

    with col1:
        region_start = st.number_input(
            "Start position (bp)",
            min_value=1,
            max_value=chromosome_length,
            value=1,
            step=1000
        )

    with col2:
        region_end = st.number_input(
            "End position (bp)",
            min_value=1,
            max_value=chromosome_length,
            value=min(200000, chromosome_length),
            step=1000
        )

    region_submit = st.form_submit_button(
        "Search Region",
        type="primary"
    )


if region_submit:

    if region_start > region_end:

        st.session_state["region_results"] = None
        st.error("Start position must be smaller than or equal to end position.")

    else:

        region_rows = []

        for gid, g in genes.items():

            if g["chrom"] != selected_chrom:
                continue

            # Overlap condition:
            # gene start <= query end AND gene end >= query start

            if (
                g["start"] <= region_end
                and g["end"] >= region_start
            ):

                tids = gene_transcripts.get(gid, [])

                region_rows.append({
                    "Gene ID": gid,
                    "Chromosome": g["chrom"],
                    "Start": g["start"],
                    "End": g["end"],
                    "Length (bp)": g["end"] - g["start"] + 1,
                    "Strand": g["strand"],
                    "Transcript count": len(tids)
                })

        region_df = pd.DataFrame(region_rows)

        if not region_df.empty:
            region_df = region_df.sort_values(
                "Start"
            ).reset_index(drop=True)

        st.session_state["region_results"] = region_df

        st.session_state["region_query"] = (
            selected_chrom,
            int(region_start),
            int(region_end)
        )


region_df = st.session_state.get("region_results")

if region_df is not None:

    chrom, qstart, qend = st.session_state["region_query"]

    st.subheader(
        f"Region Results: {chrom}:{qstart:,}-{qend:,}"
    )

    st.success(
        f"Found {len(region_df):,} overlapping genes."
    )

    st.dataframe(
        region_df,
        width="stretch",
        hide_index=True
    )

    st.download_button(
        "Download Region Gene List",
        data=region_df.to_csv(
            sep="\t",
            index=False
        ),
        file_name=f"{chrom}_{qstart}_{qend}_genes.tsv",
        mime="text/tab-separated-values",
        key="download_region_results"
    )


# ============================================================
# Batch Gene ID Search
# ============================================================

st.divider()
st.header("Batch Gene ID Search")

st.write(
    "Enter multiple G42 Gene IDs, one ID per line. "
    "Commas and semicolons are also supported."
)

batch_input = st.text_area(
    "Gene ID list",
    placeholder=(
        "ClG42_01g0000100\n"
        "ClG42_01g0000200\n"
        "ClG42_01g0000300"
    ),
    height=130,
    key="batch_gene_input"
)

if st.button(
    "Run Batch Search",
    type="primary",
    key="run_batch_search"
):

    raw_ids = (
        batch_input
        .replace(",", "\n")
        .replace(";", "\n")
    )

    input_ids = list(dict.fromkeys(raw_ids.split()))

    if not input_ids:

        st.warning("Please enter at least one Gene ID.")

        st.session_state["batch_results"] = None

    elif len(input_ids) > 500:

        st.warning(
            "A maximum of 500 Gene IDs is supported per search. "
            "Please split your list into smaller batches."
        )

        st.session_state["batch_results"] = None

    else:

        batch_rows = []
        missing_ids = []

        for gid in input_ids:

            if gid not in genes:

                missing_ids.append(gid)
                continue

            g = genes[gid]

            tids = sorted(
                gene_transcripts.get(gid, [])
            )

            descriptions = []

            for tid in tids:

                desc = annotation_dict.get(tid, "")

                if desc and desc not in descriptions:
                    descriptions.append(desc)

            batch_rows.append({

                "Gene ID": gid,

                "Chromosome": g["chrom"],

                "Start": g["start"],

                "End": g["end"],

                "Length (bp)": (
                    g["end"] - g["start"] + 1
                ),

                "Strand": g["strand"],

                "Transcript count": len(tids),

                "Transcript IDs": ";".join(tids),

                "Functional annotation": " | ".join(
                    descriptions
                )
            })

        batch_df = pd.DataFrame(batch_rows)

        st.session_state["batch_results"] = batch_df
        st.session_state["batch_missing"] = missing_ids
        st.session_state["batch_input_count"] = len(input_ids)


batch_df = st.session_state.get("batch_results")

if batch_df is not None:

    missing_ids = st.session_state.get(
        "batch_missing",
        []
    )

    total_input = st.session_state.get(
        "batch_input_count",
        0
    )

    c1, c2, c3 = st.columns(3)

    c1.metric("Input IDs", total_input)

    c2.metric("Matched genes", len(batch_df))

    c3.metric("Not found", len(missing_ids))

    st.dataframe(
        batch_df,
        width="stretch",
        hide_index=True
    )

    st.download_button(
        "Download Batch Results",
        data=batch_df.to_csv(
            sep="\t",
            index=False
        ),
        file_name="G42_batch_gene_search.tsv",
        mime="text/tab-separated-values",
        key="download_batch_results"
    )

    if missing_ids:

        with st.expander("View unmatched Gene IDs"):

            st.write(missing_ids)

gene_id = st.session_state.get("searched_gene", "")


if gene_id:

    if gene_id not in genes:

        st.error(
            f"Gene ID not found: {gene_id}"
        )

        st.stop()

    gene = genes[gene_id]

    transcript_ids = sorted(
        gene_transcripts.get(gene_id, [])
    )

    st.success(f"Gene found: {gene_id}")

    # --------------------------------------------------------
    # Gene basic information
    # --------------------------------------------------------

    st.header("1. Gene Information")

    gene_length = gene["end"] - gene["start"] + 1

    c1, c2, c3, c4 = st.columns(4)

    c1.metric("Chromosome", gene["chrom"])
    c2.metric("Gene length", f"{gene_length:,} bp")
    c3.metric("Start", f"{gene['start']:,}")
    c4.metric("End", f"{gene['end']:,}")

    st.write(
        f"**Strand:** {'Positive (+)' if gene['strand'] == '+' else 'Negative (-)'}"
    )

    st.write(
        f"**Number of transcripts:** {len(transcript_ids)}"
    )

    
    # --------------------------------------------------------
    # Chromosome location
    # --------------------------------------------------------

    st.header("2. Chromosome Location")

    st.write(
        f"Selected gene: {gene_id}"
    )

    st.write(
        f"Chromosome: {gene['chrom']} | "
        f"Position: {gene['start']:,}–{gene['end']:,} bp"
    )

    chromosome_fig = draw_chromosome_location(
        genome,
        genes,
        gene_id
    )

    st.plotly_chart(
        chromosome_fig,
        width="stretch"
    )


    # --------------------------------------------------------
    # Transcript selection
    # --------------------------------------------------------

    st.header("2. Transcript Information")

    if not transcript_ids:

        st.warning("No transcript records found for this gene.")
        st.stop()

    selected_transcript_id = st.selectbox(
        "Select transcript",
        transcript_ids
    )

    transcript = transcripts[selected_transcript_id]

    # --------------------------------------------------------
    # Transcript comparison
    # --------------------------------------------------------

    
    st.subheader("2.1 Transcript Comparison")

    comparison_df = build_transcript_comparison(
        transcript_ids,
        transcripts,
        cds_dict,
        protein_dict,
        annotation_dict
    )

    st.dataframe(
        comparison_df,
        width="stretch",
        hide_index=True
    )

    st.download_button(
        "Download Transcript Comparison",
        data=comparison_df.to_csv(
            sep="\t",
            index=False
        ),
        file_name=f"{gene_id}.transcript_comparison.tsv",
        mime="text/tab-separated-values"
    )


    # --------------------------------------------------------
    # Gene structure
    # --------------------------------------------------------

    st.header("3. Gene Structure")

    fig = draw_gene_structure(transcript)

    st.plotly_chart(
        fig,
        width="stretch"
    )

    feature_counts = pd.Series(
        [x["type"] for x in transcript["features"]]
    ).value_counts()

    st.write("**Feature statistics**")

    if len(feature_counts):
        st.dataframe(
            feature_counts.rename("Count"),
            width="stretch"
        )

    # --------------------------------------------------------
    # Functional annotation
    # --------------------------------------------------------

    st.header("4. Functional Annotation")

    annotation = annotation_dict.get(
        selected_transcript_id,
        "No functional annotation available."
    )

    st.info(annotation)

    # --------------------------------------------------------
    # Sequence information
    # --------------------------------------------------------

    st.header("5. Sequence Information")

    tab1, tab2, tab3 = st.tabs([
        "Genomic DNA",
        "CDS",
        "Protein"
    ])

    # Genomic DNA
    with tab1:

        chrom = transcript["chrom"]
        start = transcript["start"]
        end = transcript["end"]
        strand = transcript["strand"]

        genomic_sequence = str(
            genome[chrom][start - 1:end]
        ).upper()

        if strand == "-":
            genomic_sequence = str(
                Seq(genomic_sequence).reverse_complement()
            )

        st.write(
            f"Genomic DNA length: {len(genomic_sequence):,} bp"
        )

        st.text_area(
            "Genomic DNA sequence",
            value=genomic_sequence,
            height=200
        )

        st.download_button(
            "Download Genomic DNA FASTA",
            data=format_fasta(
                genomic_sequence,
                selected_transcript_id + "_genomic"
            ),
            file_name=f"{selected_transcript_id}.genomic.fa",
            mime="text/plain"
        )

    # CDS
    with tab2:

        cds_sequence = cds_dict.get(
            selected_transcript_id,
            ""
        )

        if cds_sequence:

            st.write(
                f"CDS length: {len(cds_sequence):,} bp"
            )

            st.text_area(
                "CDS sequence",
                value=cds_sequence,
                height=200
            )

            st.download_button(
                "Download CDS FASTA",
                data=format_fasta(
                    cds_sequence,
                    selected_transcript_id
                ),
                file_name=f"{selected_transcript_id}.cds.fa",
                mime="text/plain"
            )

        else:
            st.warning("CDS sequence not found.")

    # Protein
    with tab3:

        protein_sequence = protein_dict.get(
            selected_transcript_id,
            ""
        )

        if protein_sequence:

            st.write(
                f"Protein length: {len(protein_sequence):,} aa"
            )

            st.text_area(
                "Protein sequence",
                value=protein_sequence,
                height=200
            )

            st.download_button(
                "Download Protein FASTA",
                data=format_fasta(
                    protein_sequence,
                    selected_transcript_id
                ),
                file_name=f"{selected_transcript_id}.protein.fa",
                mime="text/plain"
            )

        else:
            st.warning("Protein sequence not found.")

    # --------------------------------------------------------
    # Transcript feature table
    # --------------------------------------------------------

    st.header("6. Transcript Feature Details")

    feature_rows = []

    for item in transcript["features"]:

        feature_rows.append({
            "Feature": item["type"],
            "Start": item["start"],
            "End": item["end"],
            "Length (bp)": item["end"] - item["start"] + 1,
            "Strand": item["strand"],
            "Phase": item["phase"]
        })

    feature_df = pd.DataFrame(feature_rows)

    if not feature_df.empty:

        feature_df = feature_df.sort_values(
            ["Start", "End"]
        )

        st.dataframe(
            feature_df,
            width="stretch",
            hide_index=True
        )

        st.download_button(
            "Download Feature Table",
            data=feature_df.to_csv(
                sep="\t",
                index=False
            ),
            file_name=f"{selected_transcript_id}.features.tsv",
            mime="text/tab-separated-values"
        )

    st.divider()

    st.caption(
        "G42 Genome Browser | Local development version"
    )

else:

    st.info(
        "Please enter a Gene ID above to begin searching."
    )

    st.subheader("Database Overview")

    c1, c2, c3 = st.columns(3)

    c1.metric("Genes", f"{len(genes):,}")
    c2.metric("Transcripts", f"{len(transcripts):,}")
    c3.metric("Chromosomes", len(genome.keys()))

    st.markdown(
        """
        **Example Gene ID**

        `ClG42_01g0000100`

        Enter the Gene ID in the search box to explore its structure,
        annotation and sequences.
        """
    )


