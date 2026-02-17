# Imprinting Database Sources (GRCh38)

This repository file was generated from HumanICR JBrowse track JSON endpoints:

- Root browser: `https://jb2.humanicr.org/?data=data/json/humanicr`
- Track list: `https://jb2.humanicr.org/data/json/humanicr/trackList.json`
- Tracks used:
  - `known_ICR_25`
  - `putative_ICRs_35_65_final_N1488`
  - `ICRs_2_or_more_method_val_N1348`
  - `ICRs_2_or_more_method_N4239`
- Gamete methylation tracks used to infer `methylated_parent`:
  - `Oocyte_100`
  - `Sperm_with_sue_100`

Build script:

```bash
python tools/build_imprinting_database.py
```

Output file:

- `src/methphaser/data/imprinting_regions.hg38.bed` (5501 rows)
- Parent-direction labels in current build:
  - `paternal`: 1093
  - `maternal`: 481
  - `unknown`: 3927

Licensing note:

- The HumanICR website includes a commercial-licensing notice. Verify downstream
  use requirements for your organization before redistribution/commercial use.
