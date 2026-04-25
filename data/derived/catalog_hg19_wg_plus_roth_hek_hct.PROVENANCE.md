# catalog_hg19_wg_plus_roth_hek_hct.jsonl — provenance

Append-only augmented denominator for the Roth HEK293T/HCT116 System B
pre-registration. The file is intentionally not committed because whole-genome
catalog JSONL files are gitignored and several GB on disk. This provenance file
and its SHA256 companion freeze the exact local artifact used for the tag-A
pre-registration.

## Parent Catalog

- File: `data/derived/catalog_hg19_wg.jsonl`
- SHA256: `d20661c5d5fc0c42491d9a94ef9485f69dbf071ea9a69f43a887d6f82b1357dc`
- Records: 19,787,820
- Provenance: `data/derived/catalog_hg19_wg.PROVENANCE.md`

The parent catalog is a whole-genome ThermoCas9 candidate catalog filtered to
HM450 probe windows (within 500 bp of any HM450 probe). Roth VEGFA T3 and T9
are real hg19 ThermoCas9 sites, but their nearest HM450 probes are outside
that 500 bp filter. Therefore the pre-registered catalog rule selects an
augmented denominator.

## Added Records

Two deterministic `CandidateSite` records are appended to the parent JSONL:

| target | candidate_id | hg19 critical C, 1-based | strand | PAM | family | catalog_resolution |
|---|---|---:|---|---|---|---|
| VEGFA T9 | `chr6:43742149-:NNNNCCA` | chr6:43,742,150 | - | `CCATCCA` | `NNNNCCA` | augmented_addition |
| VEGFA T3 | `chr6:43742293-:NNNNCGA` | chr6:43,742,294 | - | `CTCACGA` | `NNNNCGA` | augmented_addition |

The additions were generated from `data/raw/hg19/chr6.fa` using the same
`CandidateSite` schema and local-context convention as `thermocas.catalog`.
Each added record was validated against the anchored 23 bp spacer plus 7 bp
PAM from `data/positives/positives_roth_hek_hct_v0.tsv`:

| target | anchored sequence |
|---|---|
| VEGFA T9 | `ATCCCCTGAGAGGACAGGGAACCCCATCCA` |
| VEGFA T3 | `AGAGACCTGAACAGCGGAGAGTCCTCACGA` |

The small uncommitted local helper file
`data/derived/catalog_hg19_wg_plus_roth_hek_hct.additions.jsonl` had SHA256
`bd79cbd13ee58ce0f99a8da8b8ffde8f3dde7df0ec99086f77df0ac1005b522c` at
generation time. Its contents are fully recoverable from the table above plus
the hg19 FASTA.

## Build Commands

```
uv run python - <<'PY'
import json
from pathlib import Path
from thermocas.catalog import _local_context
from thermocas.models import CandidateSite, Strand

records = [
    {
        "candidate_id": "chr6:43742149-:NNNNCCA",
        "chrom": "chr6",
        "critical_c_pos": 43742149,
        "strand": Strand.MINUS,
        "pam": "CCATCCA",
        "pam_family": "NNNNCCA",
        "is_cpg_pam": False,
    },
    {
        "candidate_id": "chr6:43742293-:NNNNCGA",
        "chrom": "chr6",
        "critical_c_pos": 43742293,
        "strand": Strand.MINUS,
        "pam": "CTCACGA",
        "pam_family": "NNNNCGA",
        "is_cpg_pam": True,
    },
]

with open("data/raw/hg19/chr6.fa") as fh:
    fh.readline()
    seq = "".join(line.strip().upper() for line in fh)

out = Path("data/derived/catalog_hg19_wg_plus_roth_hek_hct.additions.jsonl")
with out.open("w") as f:
    for r in records:
        site = CandidateSite(
            candidate_id=r["candidate_id"],
            chrom=r["chrom"],
            critical_c_pos=r["critical_c_pos"],
            strand=r["strand"],
            pam=r["pam"],
            pam_family=r["pam_family"],
            is_cpg_pam=r["is_cpg_pam"],
            local_seq_100bp=_local_context(seq, len(seq), r["critical_c_pos"], r["strand"]),
            nearest_gene=None,
            regulatory_context=None,
        )
        f.write(json.dumps(site.model_dump(mode="json"), separators=(",", ":")) + "\n")
PY

cp data/derived/catalog_hg19_wg.jsonl \
  data/derived/catalog_hg19_wg_plus_roth_hek_hct.jsonl
cat data/derived/catalog_hg19_wg_plus_roth_hek_hct.additions.jsonl >> \
  data/derived/catalog_hg19_wg_plus_roth_hek_hct.jsonl
```

## Output Checks

```
wc -l data/derived/catalog_hg19_wg.jsonl \
      data/derived/catalog_hg19_wg_plus_roth_hek_hct.jsonl \
      data/derived/catalog_hg19_wg_plus_roth_hek_hct.additions.jsonl
```

returned:

```
19787820 data/derived/catalog_hg19_wg.jsonl
19787822 data/derived/catalog_hg19_wg_plus_roth_hek_hct.jsonl
       2 data/derived/catalog_hg19_wg_plus_roth_hek_hct.additions.jsonl
```

The augmented catalog SHA256 is:

`3ea33d1a5f294f3875ee6dcad0289ac49826c64d0dac564a0c33512cdac585fe`

## Scope

This augmented denominator is frozen only for the Roth HEK293T/HCT116 System B
pre-registration. It must not be substituted for the original frozen WG
catalog in earlier System A benchmarks unless a separate pre-registration
explicitly does so.
