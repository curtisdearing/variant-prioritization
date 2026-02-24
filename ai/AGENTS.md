# Variant Prioritization System (U4U Prototype)

## Project Overview

The **Variant Prioritization System** is a clinical variant annotation and prioritization tool designed to help users interpret genetic data. It serves as a prototype for the U4U product, allowing clients to input their genetic information and receive useful health insights.

### Key Features
*   **Multi-format Input:** Supports VCF, CSV, and TXT files (including raw 23andMe data).
*   **Automated Annotation:** Integrates with Ensembl VEP, ClinVar, and gnomAD APIs to annotate variants.
*   **Prioritization Logic:** Scores variants based on pathogenicity (ClinVar), functional consequence, and population frequency.
*   **Tiered Classification:** Categorizes variants into Critical, High, Medium, and Low tiers.
*   **N8N Integration:** Optional hook to send variants to an N8N automation pipeline for further analysis.

### Tech Stack
*   **Language:** Python 3.10+
*   **UI Framework:** Streamlit
*   **Data Processing:** Pandas, Pysam (VCF parsing)
*   **External APIs:** Ensembl REST API (VEP, Variation)
*   **Package Management:** `uv` (recommended), `pip`

## Architecture

The application is primarily contained within `main.py`, which handles the Streamlit UI, file parsing, API interactions, and scoring logic.

1.  **Input Parsing:**
    *   **VCF:** Parsed using `pysam`.
    *   **TXT/23andMe:** Custom parsers detect format and extract rsIDs/genotypes.
    *   **CSV:** Standard CSV parsing with expected columns.
2.  **Resolution & Annotation:**
    *   Variants with only rsIDs are resolved to genomic coordinates using the Ensembl API.
    *   All variants are annotated via the Ensembl VEP REST API.
3.  **Scoring & Prioritization:**
    *   **ClinVar:** Pathogenic variants get the highest scores.
    *   **Consequence:** Loss-of-Function (LoF) variants score higher than missense or silent ones.
    *   **Frequency:** Rare variants (low gnomAD AF) score higher; common ones score lower.
4.  **Output:**
    *   Interactive Streamlit dashboard.
    *   Downloadable CSV report.

## Building and Running

### Prerequisites
*   Python 3.10+
*   `uv` (recommended) or `pip`

### Setup & Execution

**Using `uv` (Recommended):**
```bash
# Install dependencies and run
./run.sh
# OR manually:
uv run streamlit run main.py
```

**Using `pip`:**
```bash
pip install -r requirements.txt
streamlit run main.py
```

The application will typically be available at `http://localhost:8501`.

## Key Files

*   `main.py`: The core application file containing the UI layout, business logic, parsers, and API clients.
*   `pyproject.toml`: Project configuration and dependencies (managed by `uv`).
*   `run.sh`: Shell script to easily run the application using `uv`.
*   `README.md`: General project documentation (Note: May reference FastAPI, but the active implementation is Streamlit).
*   `uploads/`: Temporary directory for uploaded files during processing.
*   `data/`: Directory for static resources like rsID filter lists. (Note: create this directory and add filter files like `acmg81_rsids.txt` to enable filtering features).

## Development Conventions

*   **Code Structure:** Monolithic script (`main.py`) for rapid prototyping. Future refactoring should split logic into modules (e.g., `parsers.py`, `api.py`, `scoring.py`).
*   **Error Handling:** The UI uses `st.error` and `st.warning` to communicate issues to the user gracefully.
*   **Dependencies:** `pysam` is an optional dependency for VCF support; the app guards against its absence.
*   **Configuration:** N8N integration is configured via environment variables (`N8N_WEBHOOK_URL`) or directly in the code.
