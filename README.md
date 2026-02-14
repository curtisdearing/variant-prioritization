# Variant Prioritization System

A clinical variant annotation and prioritization tool that integrates VEP, ClinVar, and gnomAD data.

## Features

- **Multi-format input**: VCF, CSV, TXT (rsID lists)
- **rsID resolution**: Automatically converts rsIDs to genomic coordinates
- **Clinical scoring**: Prioritizes variants based on:
  - ClinVar pathogenicity
  - Functional consequence (LoF > missense > silent)
  - Population frequency (gnomAD)
  - Gene context

## Tiers

- **CRITICAL**: ClinVar pathogenic/likely pathogenic (score ≥500)
- **HIGH**: Novel loss-of-function, ultra-rare deleterious (score 100-499)
- **MEDIUM**: VUS, rare missense (score 30-99)
- **LOW**: Benign, common, silent variants (score <30)

## Installation
```bash
pip install fastapi uvicorn pysam requests --break-system-packages
```

## Usage
```bash
uvicorn main:app --reload
```

Then navigate to `http://localhost:8000/docs` and upload a variant file.

## Example Input (rsID list)
```
rs80356868
rs121913254
rs113488022
```

## Example Output
```json
{
  "total": 3,
  "tiers": {
    "critical": 3,
    "high": 0,
    "medium": 0,
    "low": 0
  },
  "top_10": [
    {
      "variant_id": "17:43047661-43047661:1/A",
      "consequence": "stop_gained",
      "genes": ["BRCA1"],
      "clinvar": "pathogenic",
      "gnomad_af": null,
      "score": 1100,
      "tier": "critical",
      "reasons": [
        "⚠️ ClinVar: PATHOGENIC",
        "High impact: stop_gained",
        "Gene(s): BRCA1"
      ]
    }
  ]
}
```

## API Endpoints

- `GET /` - Health check
- `POST /upload` - Upload variant file (VCF, CSV, or TXT)

## Requirements

- Python 3.8+
- FastAPI
- pysam
- requests
