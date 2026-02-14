"""
Variant Prioritization Tool — Streamlit UI
==========================================
A user-friendly web interface for uploading variant files (VCF, TXT, CSV),
annotating them with Ensembl VEP, ClinVar, and gnomAD data, and prioritizing
them by clinical significance.

Includes clearly marked integration points for an N8N automation pipeline.

Run with:  streamlit run main.py
"""

import os
import uuid
import shutil
import csv
import re

import streamlit as st
import requests
import pandas as pd

# pysam is only needed for VCF parsing — guard the import so the rest of the
# app still works even if it's not installed (e.g. during UI development).
try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False

UPLOAD_DIR = "uploads"
os.makedirs(UPLOAD_DIR, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════
#  N8N PIPELINE INTEGRATION
# ═══════════════════════════════════════════════════════════════════════════
#
#  Your N8N pipeline accepts the NUMERIC part of an rsID (e.g. rs12345 → 12345)
#  and returns text output.  The helper below handles the call.
#
#  ──────────  HOW TO CONNECT YOUR N8N PIPELINE  ──────────
#
#  1. In N8N, create a Webhook node as the trigger for your pipeline.
#  2. Copy the webhook URL (e.g. https://your-n8n-instance.com/webhook/abc123)
#  3. Paste it into N8N_WEBHOOK_URL below (or set the env-var).
#  4. The function sends a POST request with JSON body:
#         { "rsid_number": "12345" }
#     Make sure your N8N workflow's Webhook node accepts POST + JSON.
#  5. The pipeline should return a JSON response.  We read the text from
#     the field specified by N8N_RESPONSE_KEY (default "output").
#
#  ──────────  CONFIGURATION  ──────────

# TODO: Replace with your actual N8N webhook URL, or set the
#       N8N_WEBHOOK_URL environment variable.
N8N_WEBHOOK_URL = os.environ.get("N8N_WEBHOOK_URL", "")

# TODO: If your N8N pipeline returns the text in a different JSON field,
#       change this key.  For example, if the response is {"result": "..."},
#       set this to "result".
N8N_RESPONSE_KEY = os.environ.get("N8N_RESPONSE_KEY", "output")

# TODO: Adjust the timeout (in seconds) for your N8N pipeline.
N8N_TIMEOUT_SECONDS = int(os.environ.get("N8N_TIMEOUT_SECONDS", "120"))


def extract_rsid_number(rsid: str) -> str | None:
    """
    Extracts the numeric portion from an rsID string.
    Examples:
        "rs12345"  →  "12345"
        "RS429358" →  "429358"
        "12345"    →  "12345"   (already numeric)
    Returns None if no number can be found.
    """
    if rsid is None:
        return None
    match = re.search(r"(\d+)", rsid)
    return match.group(1) if match else None


def call_n8n_pipeline(rsid: str) -> str:
    """
    Sends the numeric part of *rsid* to the N8N webhook and returns
    the text response from the pipeline.

    ──── Integration Notes ────
    • Input sent to N8N:   { "rsid_number": "<numeric_id>" }
    • Expected response:   { "<N8N_RESPONSE_KEY>": "<text output>" }

    If the webhook URL is not configured, or the call fails, a
    descriptive placeholder message is returned instead.
    """
    rsid_number = extract_rsid_number(rsid)
    if rsid_number is None:
        return "⚠️ Could not extract a numeric ID from this rsID."

    if not N8N_WEBHOOK_URL:
        return (
            "🔧 **N8N webhook not configured.**\n\n"
            "Set `N8N_WEBHOOK_URL` in `main.py` or as an environment variable "
            "to enable the pipeline."
        )

    # ── Build the request payload ──
    # TODO: Adjust the payload structure if your N8N pipeline expects
    #       different fields or additional data.
    payload = {
        "rsid_number": rsid_number,
    }

    # TODO: Add authentication headers if your N8N instance requires them.
    #       For example:
    #       headers = {
    #           "Content-Type": "application/json",
    #           "Authorization": "Bearer <YOUR_TOKEN>",
    #       }
    headers = {
        "Content-Type": "application/json",
    }

    try:
        response = requests.post(
            N8N_WEBHOOK_URL,
            json=payload,
            headers=headers,
            timeout=N8N_TIMEOUT_SECONDS,
        )
        response.raise_for_status()

        data = response.json()

        # TODO: Parse the response.  By default we look for the key
        #       defined by N8N_RESPONSE_KEY.  Change this logic if your
        #       pipeline returns a different structure.
        if isinstance(data, dict):
            return data.get(N8N_RESPONSE_KEY, str(data))
        elif isinstance(data, list) and len(data) > 0:
            # Some N8N webhook responses are wrapped in a list
            first = data[0]
            if isinstance(first, dict):
                return first.get(N8N_RESPONSE_KEY, str(first))
            return str(first)
        else:
            return str(data)

    except requests.exceptions.ConnectionError:
        return "❌ Could not connect to N8N. Is the workflow active?"
    except requests.exceptions.Timeout:
        return f"⏱️ N8N pipeline timed out after {N8N_TIMEOUT_SECONDS}s."
    except requests.exceptions.HTTPError as e:
        return f"❌ N8N returned an error: {e}"
    except Exception as e:
        return f"❌ Unexpected error calling N8N: {e}"


# ═══════════════════════════════════════════════════════════════════════════
#  FILE PARSERS
# ═══════════════════════════════════════════════════════════════════════════

def parse_vcf(file_path):
    """Parse a VCF file using pysam and return a list of variant dicts."""
    if not PYSAM_AVAILABLE:
        raise ImportError(
            "pysam is required to parse VCF files. "
            "Install it with: pip install pysam"
        )
    variants = []
    vcf = pysam.VariantFile(file_path)
    for record in vcf.fetch():
        for alt in record.alts:
            variants.append({
                "chrom": record.chrom,
                "pos": record.pos,
                "ref": record.ref,
                "alt": alt,
                "rsid": record.id,
                "variant_type": "coordinate",
            })
    return variants


def _is_23andme_format(file_path):
    """Detect whether a .txt file is in 23andMe raw-data format.

    23andMe files have a long comment header (lines starting with #)
    followed by a column-header line and tab-separated data:

        # This data file generated by 23andMe at: ...
        # ...
        # rsid	chromosome	position	genotype
        rs548049170	1	69869	TT

    We scan past the comment block looking for either the specific
    column header or a tab-separated data line.
    """
    with open(file_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue

            # The 23andMe column header
            if stripped.lower().startswith("# rsid"):
                return True

            # Skip other comment/preamble lines and keep looking
            if stripped.startswith("#"):
                continue

            # First non-comment line — check if it's tab-separated data
            parts = stripped.split("\t")
            if len(parts) >= 4 and parts[0].lower().startswith("rs"):
                return True

            # Non-comment line that doesn't look like 23andMe data
            return False
    return False


def parse_rsid_txt(file_path):
    """Parse a plain-text file with one rsID per line."""
    variants = []
    with open(file_path, "r") as f:
        for line in f:
            rsid = line.strip()
            if rsid:
                variants.append({
                    "chrom": None,
                    "pos": None,
                    "ref": None,
                    "alt": None,
                    "rsid": rsid,
                    "variant_type": "rsid_only",
                })
    return variants


def parse_23andme(file_path):
    """Parse a 23andMe raw data file (tab-separated).

    Format:
        # rsid  chromosome  position  genotype
        rs548049170  1  69869  TT
        rs9283150    1  565508 AA

    Since 23andMe files contain the subject's *genotype* (e.g. "TT")
    but not the reference allele, we treat every variant as rsID-only
    and let Ensembl resolve the full ref/alt coordinates.
    """
    variants = []
    with open(file_path, "r") as f:
        for line in f:
            stripped = line.strip()
            # Skip blanks and comment/header lines
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split("\t")
            if len(parts) < 4:
                continue
            rsid, chrom, pos, genotype = parts[0], parts[1], parts[2], parts[3]
            # Skip internal 23andMe IDs (e.g. "i7001348") — only process rs IDs
            if not rsid.lower().startswith("rs"):
                continue
            variants.append({
                "chrom": chrom,
                "pos": int(pos) if pos else None,
                "ref": None,
                "alt": None,
                "rsid": rsid,
                "genotype": genotype,
                "variant_type": "rsid_only",  # resolved via Ensembl
            })
    return variants


def parse_csv(file_path):
    """Parse a CSV file with optional columns: chrom, pos, ref, alt, rsid."""
    variants = []
    with open(file_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            variant_type = (
                "coordinate" if row.get("chrom") and row.get("pos") else "rsid_only"
            )
            variants.append({
                "chrom": row.get("chrom"),
                "pos": int(row.get("pos")) if row.get("pos") else None,
                "ref": row.get("ref"),
                "alt": row.get("alt"),
                "rsid": row.get("rsid"),
                "variant_type": variant_type,
            })
    return variants


# ═══════════════════════════════════════════════════════════════════════════
#  RSID RESOLUTION  (Ensembl REST API)
# ═══════════════════════════════════════════════════════════════════════════

def resolve_rsid(rsid):
    """Convert an rsID to genomic coordinates using the Ensembl REST API."""
    response = requests.get(
        f"https://rest.ensembl.org/variation/human/{rsid}",
        headers={"Content-Type": "application/json"},
    )
    if response.status_code != 200:
        return None

    data = response.json()

    mappings = data.get("mappings", [])
    if not mappings:
        return None

    mapping = mappings[0]

    allele_string = mapping.get("allele_string", "")
    alleles = allele_string.split("/")
    if len(alleles) < 2:
        return None

    ref = alleles[0]
    alts = alleles[1:]

    resolved = []
    for alt in alts:
        resolved.append({
            "chrom": mapping["seq_region_name"],
            "pos": mapping["start"],
            "ref": ref,
            "alt": alt,
            "rsid": rsid,
            "variant_type": "coordinate",
        })
    return resolved


# ═══════════════════════════════════════════════════════════════════════════
#  CLINVAR & GNOMAD LOOKUP  (extracted from VEP colocated_variants)
# ═══════════════════════════════════════════════════════════════════════════

def get_clinvar_significance(vep_result):
    """Extract ClinVar clinical significance from a VEP result."""
    colocated = vep_result.get("colocated_variants", [])
    for variant in colocated:
        clin_sig = variant.get("clin_sig", [])
        if clin_sig:
            if "pathogenic" in clin_sig:
                return "pathogenic"
            elif "likely_pathogenic" in clin_sig:
                return "likely_pathogenic"
            elif "uncertain_significance" in clin_sig:
                return "uncertain_significance"
            elif "likely_benign" in clin_sig:
                return "likely_benign"
            elif "benign" in clin_sig:
                return "benign"
    return None


def get_gnomad_frequency(vep_result):
    """Extract gnomAD allele frequency from a VEP result."""
    colocated = vep_result.get("colocated_variants", [])
    for variant in colocated:
        frequencies = variant.get("frequencies", {})
        for allele, freq_data in frequencies.items():
            gnomadg_af = freq_data.get("gnomadg")
            if gnomadg_af is not None:
                return gnomadg_af
            gnomade_af = freq_data.get("gnomade")
            if gnomade_af is not None:
                return gnomade_af
    return None


# ═══════════════════════════════════════════════════════════════════════════
#  VEP ANNOTATION  (Ensembl REST API)
# ═══════════════════════════════════════════════════════════════════════════

def run_vep_api(variants, progress_callback=None):
    """Annotate a list of coordinate-resolved variants via the VEP REST API."""
    annotated = []
    for i, v in enumerate(variants):
        chrom = str(v["chrom"]).replace("chr", "")
        region_string = f"{chrom}:{v['pos']}-{v['pos']}:1/{v['alt']}"
        response = requests.post(
            "https://rest.ensembl.org/vep/human/region",
            headers={"Content-Type": "application/json"},
            json={"variants": [region_string]},
        )
        if response.status_code != 200:
            raise Exception(
                f"VEP API failed for variant {region_string}: {response.text}"
            )
        result = response.json()
        if result:
            # Carry the original rsID through so we can use it for N8N
            result[0]["_original_rsid"] = v.get("rsid")
            annotated.append(result[0])

        if progress_callback:
            progress_callback(i + 1, len(variants))

    return annotated


# ═══════════════════════════════════════════════════════════════════════════
#  VARIANT PRIORITIZATION  (scoring & tiering)
# ═══════════════════════════════════════════════════════════════════════════

HIGH_IMPACT = {
    "stop_gained", "frameshift_variant", "splice_donor_variant",
    "splice_acceptor_variant", "start_lost", "stop_lost",
}
MODERATE_IMPACT = {
    "missense_variant", "inframe_deletion", "inframe_insertion",
}
LOW_IMPACT = {
    "synonymous_variant", "intron_variant", "intergenic_variant",
}


def prioritize_variant(vep_result):
    """
    Score a variant based on clinical significance:
      1. ClinVar classification  (highest priority)
      2. Functional consequence  (LoF > missense > silent)
      3. Population frequency    (rare ≈ more likely pathogenic)
    Returns a dict with score, tier, reasons, and annotation data.
    """
    score = 0
    reasons = []

    # Extract consequence terms and gene symbols
    consequences = set()
    genes = set()
    for tc in vep_result.get("transcript_consequences", []):
        consequences.update(tc.get("consequence_terms", []))
        if tc.get("gene_symbol"):
            genes.add(tc["gene_symbol"])

    clinvar = get_clinvar_significance(vep_result)
    gnomad_af = get_gnomad_frequency(vep_result)

    # ---------- 1. ClinVar ----------
    if clinvar == "pathogenic":
        score += 1000
        reasons.append("⚠️ ClinVar: PATHOGENIC")
    elif clinvar == "likely_pathogenic":
        score += 500
        reasons.append("⚠️ ClinVar: Likely Pathogenic")
    elif clinvar == "benign":
        score = 1
        reasons.append("✅ ClinVar: Benign")
        return _build_result(vep_result, genes, clinvar, gnomad_af, score, "low", reasons)
    elif clinvar == "likely_benign":
        score = 5
        reasons.append("✅ ClinVar: Likely Benign")
    elif clinvar == "uncertain_significance":
        score += 50
        reasons.append("❓ ClinVar: VUS (uncertain)")

    # ---------- 2. Functional consequence ----------
    if consequences & HIGH_IMPACT:
        score += 100
        reasons.append(f"High impact: {', '.join(consequences & HIGH_IMPACT)}")
    elif consequences & MODERATE_IMPACT:
        score += 50
        reasons.append(f"Moderate impact: {', '.join(consequences & MODERATE_IMPACT)}")
    elif consequences & LOW_IMPACT:
        score += 5
        reasons.append(f"Low impact: {', '.join(consequences & LOW_IMPACT)}")
    else:
        score += 1
        reasons.append("No coding consequence")

    # ---------- 3. Population frequency ----------
    if gnomad_af is not None:
        if gnomad_af == 0:
            score += 30
            reasons.append("Absent in gnomAD (AF=0)")
        elif gnomad_af < 0.0001:
            score += 20
            reasons.append(f"Ultra-rare (AF={gnomad_af:.6f})")
        elif gnomad_af < 0.001:
            score += 10
            reasons.append(f"Very rare (AF={gnomad_af:.4f})")
        elif gnomad_af < 0.01:
            score += 5
            reasons.append(f"Rare (AF={gnomad_af:.3f})")
        else:
            score -= 20
            reasons.append(f"Common variant (AF={gnomad_af:.2%})")
    else:
        reasons.append("No frequency data")

    # ---------- 4. Gene presence ----------
    if genes:
        reasons.append(f"Gene(s): {', '.join(sorted(genes))}")
    else:
        score -= 10
        reasons.append("Intergenic (no gene)")

    # ---------- Tiering ----------
    if score >= 500:
        tier = "critical"
    elif score >= 100:
        tier = "high"
    elif score >= 30:
        tier = "medium"
    else:
        tier = "low"

    return _build_result(vep_result, genes, clinvar, gnomad_af, score, tier, reasons)


def _build_result(vep_result, genes, clinvar, gnomad_af, score, tier, reasons):
    """Helper to build a standardised result dict."""
    return {
        "variant_id": vep_result.get("id"),
        "rsid": vep_result.get("_original_rsid"),
        "location": f"{vep_result['seq_region_name']}:{vep_result['start']}",
        "consequence": vep_result.get("most_severe_consequence", "unknown"),
        "genes": list(genes),
        "clinvar": clinvar,
        "gnomad_af": gnomad_af,
        "score": score,
        "tier": tier,
        "reasons": reasons,
    }


# ═══════════════════════════════════════════════════════════════════════════
#  STREAMLIT UI
# ═══════════════════════════════════════════════════════════════════════════

TIER_COLORS = {
    "critical": "🔴",
    "high": "🟠",
    "medium": "🟡",
    "low": "🟢",
}

TIER_CSS = {
    "critical": "#ef4444",
    "high": "#f97316",
    "medium": "#eab308",
    "low": "#22c55e",
}


def render_methodology_sidebar():
    """Render an explanatory sidebar describing the scoring methodology."""
    with st.sidebar:
        st.header("ℹ️ How It Works")
        st.markdown(
            """
            **This tool prioritizes genetic variants** by combining
            multiple evidence sources into a single actionable score.

            ---

            #### 📂 Accepted File Formats
            - **VCF** — standard variant call format
            - **TXT** — one rsID per line (e.g. `rs7412`)
            - **TXT (23andMe)** — raw data download from 23andMe
            - **CSV** — columns: `chrom`, `pos`, `ref`, `alt`, `rsid`

            ---

            #### 🧬 Scoring Components

            | Source | Weight |
            |--------|--------|
            | ClinVar pathogenic | +1 000 |
            | ClinVar likely path. | +500 |
            | ClinVar VUS | +50 |
            | High-impact (LoF) | +100 |
            | Moderate-impact | +50 |
            | Rare in gnomAD | +10 – +30 |
            | Common in gnomAD | −20 |

            ---

            #### 🏷️ Tier Definitions

            | Tier | Score | Meaning |
            |------|-------|---------|
            | 🔴 Critical | ≥ 500 | ClinVar pathogenic |
            | 🟠 High | ≥ 100 | Likely LoF / likely path. |
            | 🟡 Medium | ≥ 30 | Missense + rare, or VUS |
            | 🟢 Low | < 30 | Common / benign / silent |
            """
        )

        st.divider()

        # ── N8N toggle ──
        st.subheader("🔗 N8N Pipeline")
        n8n_enabled = st.toggle(
            "Enable N8N analysis",
            value=False,
            help="When enabled, each variant with an rsID will be sent to "
                 "your N8N pipeline for additional analysis.",
        )
        if n8n_enabled and not N8N_WEBHOOK_URL:
            st.warning(
                "N8N is enabled but no webhook URL is configured. "
                "Set `N8N_WEBHOOK_URL` in main.py or as an env variable."
            )
        return n8n_enabled


def render_header():
    """Render the page header with an explanation of the tool."""
    st.markdown(
        """
        <h1 style="text-align:center; margin-bottom:0;">
            🧬 Variant Prioritization Tool
        </h1>
        <p style="text-align:center; color:#888; margin-top:4px;">
            Upload your variant file and get clinically prioritized results in seconds.
        </p>
        """,
        unsafe_allow_html=True,
    )
    st.divider()

    with st.expander("📖 **What does this tool do?**", expanded=False):
        st.markdown(
            """
            This tool helps you **identify which genetic variants matter most**
            from a clinical perspective.

            1. **Upload** a VCF, TXT, or CSV file containing your variants.
            2. Variants with only rsIDs are automatically **resolved** to
               genomic coordinates via the Ensembl REST API.
            3. Each variant is **annotated** using the Ensembl VEP (Variant
               Effect Predictor), which provides functional consequence,
               ClinVar classifications, and gnomAD population frequencies.
            4. A **prioritization score** is computed and each variant is
               assigned a tier — from 🔴 Critical to 🟢 Low.
            5. Optionally, each variant can be sent through an **N8N
               automation pipeline** for additional custom analysis.
            """
        )


def render_metrics(prioritized):
    """Render summary metric cards."""
    tiers = {"critical": 0, "high": 0, "medium": 0, "low": 0}
    for v in prioritized:
        tiers[v["tier"]] += 1

    cols = st.columns(5)
    cols[0].metric("Total Variants", len(prioritized))
    cols[1].metric("🔴 Critical", tiers["critical"])
    cols[2].metric("🟠 High", tiers["high"])
    cols[3].metric("🟡 Medium", tiers["medium"])
    cols[4].metric("🟢 Low", tiers["low"])


def render_variant_card(variant, n8n_enabled):
    """Render a single variant as an expandable card."""
    tier_emoji = TIER_COLORS.get(variant["tier"], "⚪")
    tier_color = TIER_CSS.get(variant["tier"], "#888")

    label = (
        f"{tier_emoji} **{variant['variant_id'] or variant['location']}** "
        f"— Score {variant['score']}  "
        f"(`{variant['tier'].upper()}`)"
    )

    with st.expander(label, expanded=(variant["tier"] in ("critical", "high"))):
        c1, c2, c3 = st.columns(3)
        c1.markdown(f"**Location:** `{variant['location']}`")
        c2.markdown(f"**Consequence:** `{variant['consequence']}`")
        c3.markdown(f"**Genes:** {', '.join(variant['genes']) if variant['genes'] else '—'}")

        c4, c5, c6 = st.columns(3)
        c4.markdown(f"**ClinVar:** {variant['clinvar'] or '—'}")
        c5.markdown(
            f"**gnomAD AF:** "
            f"{variant['gnomad_af']:.6f}" if variant["gnomad_af"] is not None else "**gnomAD AF:** —"
        )
        c6.markdown(f"**rsID:** {variant['rsid'] or '—'}")

        st.markdown("**Scoring Reasons:**")
        for r in variant["reasons"]:
            st.markdown(f"- {r}")

        # ── N8N Analysis ──
        if variant.get("rsid"):
            st.divider()
            if n8n_enabled:
                with st.spinner("Running N8N pipeline…"):
                    n8n_result = call_n8n_pipeline(variant["rsid"])
                st.markdown("**🔗 N8N Pipeline Analysis:**")
                st.info(n8n_result)
            else:
                st.caption(
                    "💡 Enable the N8N toggle in the sidebar to run "
                    "additional analysis on this variant."
                )


def render_download(prioritized):
    """Offer a CSV download of all results."""
    rows = []
    for v in prioritized:
        rows.append({
            "Variant ID": v["variant_id"],
            "rsID": v["rsid"],
            "Location": v["location"],
            "Consequence": v["consequence"],
            "Genes": ", ".join(v["genes"]) if v["genes"] else "",
            "ClinVar": v["clinvar"] or "",
            "gnomAD AF": v["gnomad_af"] if v["gnomad_af"] is not None else "",
            "Score": v["score"],
            "Tier": v["tier"],
            "Reasons": " | ".join(v["reasons"]),
        })
    df = pd.DataFrame(rows)
    st.download_button(
        label="📥 Download results as CSV",
        data=df.to_csv(index=False),
        file_name="variant_prioritization_results.csv",
        mime="text/csv",
    )


# ═══════════════════════════════════════════════════════════════════════════
#  MAIN  (Streamlit page)
# ═══════════════════════════════════════════════════════════════════════════

def main():
    st.set_page_config(
        page_title="Variant Prioritization",
        page_icon="🧬",
        layout="wide",
    )

    n8n_enabled = render_methodology_sidebar()
    render_header()

    # ── File uploader ──
    uploaded_file = st.file_uploader(
        "Upload a variant file",
        type=["vcf", "txt", "csv"],
        help="VCF, plain-text rsID list, 23andMe raw data (.txt), or CSV with variant columns.",
    )

    if uploaded_file is None:
        st.info("👆 Upload a file to get started.")
        return

    # ── Save uploaded file to disk (parsers expect a path) ──
    job_id = str(uuid.uuid4())
    job_dir = os.path.join(UPLOAD_DIR, job_id)
    os.makedirs(job_dir, exist_ok=True)
    file_path = os.path.join(job_dir, uploaded_file.name)
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())

    filename = uploaded_file.name.lower()

    # ── Parse ──
    with st.spinner("Parsing file…"):
        try:
            if filename.endswith(".vcf"):
                variants = parse_vcf(file_path)
            elif filename.endswith(".txt"):
                # Auto-detect 23andMe format vs plain rsID list
                if _is_23andme_format(file_path):
                    variants = parse_23andme(file_path)
                else:
                    variants = parse_rsid_txt(file_path)
            elif filename.endswith(".csv"):
                variants = parse_csv(file_path)
            else:
                st.error("Unsupported file format.")
                return
        except Exception as e:
            st.error(f"Failed to parse file: {e}")
            return

    st.success(f"Parsed **{len(variants)}** variant(s) from `{uploaded_file.name}`.")

    # ── Separate coordinate vs rsID-only variants ──
    coordinate_variants = [v for v in variants if v["variant_type"] == "coordinate"]
    rsid_only_variants = [v for v in variants if v["variant_type"] == "rsid_only"]

    # ── Resolve rsIDs ──
    if rsid_only_variants:
        st.markdown(f"Resolving **{len(rsid_only_variants)}** rsID(s) to genomic coordinates…")
        progress_bar = st.progress(0)
        resolved_count = 0
        failed_rsids = []

        for i, rsid_variant in enumerate(rsid_only_variants):
            rsid = rsid_variant["rsid"]
            resolved = resolve_rsid(rsid)
            if resolved:
                coordinate_variants.extend(resolved)
                resolved_count += 1
            else:
                failed_rsids.append(rsid)
            progress_bar.progress((i + 1) / len(rsid_only_variants))

        if resolved_count:
            st.success(f"Resolved {resolved_count} rsID(s).")
        if failed_rsids:
            st.warning(f"Could not resolve {len(failed_rsids)} rsID(s): {', '.join(failed_rsids[:10])}")

    if not coordinate_variants:
        st.warning("No coordinate variants available for annotation.")
        return

    # ── VEP annotation ──
    st.markdown(f"Annotating **{len(coordinate_variants)}** variant(s) with Ensembl VEP…")
    progress_bar = st.progress(0)

    def vep_progress(current, total):
        progress_bar.progress(current / total)

    try:
        vep_results = run_vep_api(coordinate_variants, progress_callback=vep_progress)
    except Exception as e:
        st.error(f"VEP annotation failed: {e}")
        return

    # ── Prioritize ──
    with st.spinner("Scoring and prioritizing…"):
        prioritized = [prioritize_variant(v) for v in vep_results]
        prioritized.sort(key=lambda x: x["score"], reverse=True)

    # ── Display results ──
    st.divider()
    st.subheader("📊 Results")
    render_metrics(prioritized)
    st.divider()

    # Download button
    render_download(prioritized)

    st.divider()

    # Variant cards
    st.subheader("🔎 Variant Details")
    for v in prioritized:
        render_variant_card(v, n8n_enabled)

    # Cleanup temp files
    shutil.rmtree(job_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
