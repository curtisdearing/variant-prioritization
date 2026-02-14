from fastapi import FastAPI, UploadFile, File
import shutil
import os
import uuid
import pysam
import csv
import requests

app = FastAPI()

UPLOAD_DIR = "uploads"
os.makedirs(UPLOAD_DIR, exist_ok=True)


# -----------------------------
# PARSERS (unchanged)
# -----------------------------

def parse_vcf(file_path):
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
                "variant_type": "coordinate"
            })
    return variants


def parse_rsid_txt(file_path):
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
                    "variant_type": "rsid_only"
                })
    return variants


def parse_csv(file_path):
    variants = []
    with open(file_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            variant_type = "coordinate" if row.get("chrom") and row.get("pos") else "rsid_only"
            variants.append({
                "chrom": row.get("chrom"),
                "pos": int(row.get("pos")) if row.get("pos") else None,
                "ref": row.get("ref"),
                "alt": row.get("alt"),
                "rsid": row.get("rsid"),
                "variant_type": variant_type
            })
    return variants


# -----------------------------
# RSID RESOLUTION
# -----------------------------

def resolve_rsid(rsid):
    """Convert rsID to genomic coordinates using Ensembl API"""
    response = requests.get(
        f"https://rest.ensembl.org/variation/human/{rsid}",
        headers={"Content-Type": "application/json"}
    )
    
    if response.status_code != 200:
        return None
    
    data = response.json()
    
    # Get GRCh38 mappings
    mappings = data.get('mappings', [])
    if not mappings:
        return None
    
    # Use first mapping (primary assembly)
    mapping = mappings[0]
    
    # Parse alleles
    allele_string = mapping.get('allele_string', '')
    alleles = allele_string.split('/')
    if len(alleles) < 2:
        return None
    
    ref = alleles[0]
    alts = alleles[1:]  # Can be multiple alt alleles
    
    resolved = []
    for alt in alts:
        resolved.append({
            'chrom': mapping['seq_region_name'],
            'pos': mapping['start'],
            'ref': ref,
            'alt': alt,
            'rsid': rsid,
            'variant_type': 'coordinate'
        })
    
    return resolved


# -----------------------------
# CLINVAR & GNOMAD LOOKUP
# -----------------------------

def get_clinvar_significance(vep_result):
    """Extract ClinVar clinical significance from VEP result"""
    # VEP includes ClinVar data in colocated_variants
    colocated = vep_result.get('colocated_variants', [])
    
    for variant in colocated:
        clin_sig = variant.get('clin_sig', [])
        if clin_sig:
            # Return most severe classification
            if 'pathogenic' in clin_sig:
                return 'pathogenic'
            elif 'likely_pathogenic' in clin_sig:
                return 'likely_pathogenic'
            elif 'uncertain_significance' in clin_sig:
                return 'uncertain_significance'
            elif 'likely_benign' in clin_sig:
                return 'likely_benign'
            elif 'benign' in clin_sig:
                return 'benign'
    
    return None


def get_gnomad_frequency(vep_result):
    """Extract gnomAD allele frequency from VEP result"""
    # VEP includes gnomAD frequencies in colocated_variants
    colocated = vep_result.get('colocated_variants', [])
    
    for variant in colocated:
        frequencies = variant.get('frequencies', {})
        
        # Get the variant's allele
        for allele, freq_data in frequencies.items():
            # gnomAD v3 genome frequency
            gnomadg_af = freq_data.get('gnomadg')
            if gnomadg_af is not None:
                return gnomadg_af
            
            # Fallback to gnomAD exome
            gnomade_af = freq_data.get('gnomade')
            if gnomade_af is not None:
                return gnomade_af
    
    return None


# -----------------------------
# VEP API
# -----------------------------

def run_vep_api(variants):
    annotated = []
    for v in variants:
        chrom = str(v['chrom']).replace('chr', '')
        region_string = f"{chrom}:{v['pos']}-{v['pos']}:1/{v['alt']}"
        response = requests.post(
            "https://rest.ensembl.org/vep/human/region",
            headers={"Content-Type": "application/json"},
            json={"variants": [region_string]}
        )
        if response.status_code != 200:
            raise Exception(f"VEP API failed for variant {region_string}: {response.text}")
        result = response.json()
        if result:
            annotated.append(result[0])
    return annotated


# -----------------------------
# PRIORITIZATION WITH CLINVAR & GNOMAD
# -----------------------------

HIGH_IMPACT = {'stop_gained', 'frameshift_variant', 'splice_donor_variant', 
               'splice_acceptor_variant', 'start_lost', 'stop_lost'}
MODERATE_IMPACT = {'missense_variant', 'inframe_deletion', 'inframe_insertion'}
LOW_IMPACT = {'synonymous_variant', 'intron_variant', 'intergenic_variant'}


def prioritize_variant(vep_result):
    """
    Score variant based on what actually matters clinically:
    1. ClinVar classification (highest priority)
    2. Functional consequence (loss of function > missense > silent)
    3. Population frequency (rare = more likely pathogenic)
    """
    score = 0
    reasons = []
    
    # Extract data
    consequences = set()
    genes = set()
    
    for tc in vep_result.get('transcript_consequences', []):
        consequences.update(tc.get('consequence_terms', []))
        if tc.get('gene_symbol'):
            genes.add(tc['gene_symbol'])
    
    clinvar = get_clinvar_significance(vep_result)
    gnomad_af = get_gnomad_frequency(vep_result)
    
    # ============================================
    # SCORING LOGIC (what matters most to least)
    # ============================================
    
    # 1. CLINVAR - Overrides everything (if present)
    if clinvar == 'pathogenic':
        score += 1000
        reasons.append("⚠️ ClinVar: PATHOGENIC")
    elif clinvar == 'likely_pathogenic':
        score += 500
        reasons.append("⚠️ ClinVar: Likely Pathogenic")
    elif clinvar == 'benign':
        score = 1  # Known benign - ignore everything else
        reasons.append("✓ ClinVar: Benign")
        return {
            'variant_id': vep_result.get('id'),
            'location': f"{vep_result['seq_region_name']}:{vep_result['start']}",
            'consequence': vep_result.get('most_severe_consequence', 'unknown'),
            'genes': list(genes),
            'clinvar': clinvar,
            'gnomad_af': gnomad_af,
            'score': score,
            'tier': 'low',
            'reasons': reasons
        }
    elif clinvar == 'likely_benign':
        score = 5
        reasons.append("✓ ClinVar: Likely Benign")
    elif clinvar == 'uncertain_significance':
        score += 50
        reasons.append("? ClinVar: VUS (uncertain)")
    
    # 2. FUNCTIONAL CONSEQUENCE
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
    
    # 3. POPULATION FREQUENCY (rarity = pathogenicity signal)
    if gnomad_af is not None:
        if gnomad_af == 0:
            score += 30
            reasons.append(f"Absent in gnomAD (AF=0)")
        elif gnomad_af < 0.0001:  # <0.01%
            score += 20
            reasons.append(f"Ultra-rare (AF={gnomad_af:.6f})")
        elif gnomad_af < 0.001:  # <0.1%
            score += 10
            reasons.append(f"Very rare (AF={gnomad_af:.4f})")
        elif gnomad_af < 0.01:  # <1%
            score += 5
            reasons.append(f"Rare (AF={gnomad_af:.3f})")
        else:  # Common
            score -= 20
            reasons.append(f"Common variant (AF={gnomad_af:.2%})")
    else:
        reasons.append("No frequency data")
    
    # 4. GENE PRESENCE
    if genes:
        reasons.append(f"Gene(s): {', '.join(sorted(genes))}")
    else:
        score -= 10
        reasons.append("Intergenic (no gene)")
    
    # ============================================
    # TIERING (based on clinical actionability)
    # ============================================
    
    if score >= 500:
        tier = 'critical'  # ClinVar pathogenic
    elif score >= 100:
        tier = 'high'  # LoF or likely pathogenic
    elif score >= 30:
        tier = 'medium'  # Missense + rare, or VUS
    else:
        tier = 'low'  # Common, silent, or benign
    
    return {
        'variant_id': vep_result.get('id'),
        'location': f"{vep_result['seq_region_name']}:{vep_result['start']}",
        'consequence': vep_result.get('most_severe_consequence', 'unknown'),
        'genes': list(genes),
        'clinvar': clinvar,
        'gnomad_af': gnomad_af,
        'score': score,
        'tier': tier,
        'reasons': reasons
    }


# -----------------------------
# MAIN ROUTE (updated)
# -----------------------------

@app.post("/upload")
async def upload_file(file: UploadFile = File(...)):
    job_id = str(uuid.uuid4())
    job_dir = os.path.join(UPLOAD_DIR, job_id)
    os.makedirs(job_dir, exist_ok=True)
    
    file_path = os.path.join(job_dir, file.filename)
    with open(file_path, "wb") as buffer:
        shutil.copyfileobj(file.file, buffer)
    
    filename = file.filename.lower()
    
    try:
        if filename.endswith(".vcf"):
            variants = parse_vcf(file_path)
        elif filename.endswith(".txt"):
            variants = parse_rsid_txt(file_path)
        elif filename.endswith(".csv"):
            variants = parse_csv(file_path)
        else:
            return {"error": "Unsupported file format"}
    except Exception as e:
        return {"error": str(e)}
    
    # Separate coordinate and rsid-only variants
    coordinate_variants = [v for v in variants if v["variant_type"] == "coordinate"]
    rsid_only_variants = [v for v in variants if v["variant_type"] == "rsid_only"]
    
    # Resolve rsIDs to coordinates
    resolution_stats = {'resolved': 0, 'failed': 0}
    for rsid_variant in rsid_only_variants:
        rsid = rsid_variant['rsid']
        resolved = resolve_rsid(rsid)
        if resolved:
            coordinate_variants.extend(resolved)
            resolution_stats['resolved'] += 1
        else:
            resolution_stats['failed'] += 1
    
    if coordinate_variants:
        try:
            vep_results = run_vep_api(coordinate_variants)
            
            # PRIORITIZE
            prioritized = [prioritize_variant(v) for v in vep_results]
            prioritized.sort(key=lambda x: x['score'], reverse=True)
            
            # Count tiers
            tiers = {'critical': 0, 'high': 0, 'medium': 0, 'low': 0}
            for v in prioritized:
                tiers[v['tier']] += 1
            
            response = {
                "job_id": job_id,
                "total": len(prioritized),
                "tiers": tiers,
                "top_10": prioritized[:10]
            }
            
            # Add resolution stats if any rsIDs were processed
            if rsid_only_variants:
                response['rsid_resolution'] = resolution_stats
            
            return response
            
        except Exception as e:
            return {"error": f"VEP/Prioritization failed: {str(e)}"}
    
    return {
        "job_id": job_id,
        "message": "No coordinate variants found",
        "variants_preview": variants[:10]
    }


@app.get("/")
def read_root():
    return {"message": "Backend is working"}
