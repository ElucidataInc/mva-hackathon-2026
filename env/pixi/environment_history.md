
## Environment verification 2026-09-30T03:58:57Z

Log: /home/ubuntu/mva-hackathon-2026/tools/pixi/logs/verify_20260930T035853Z_NadCJt/verification.txt
Failed checks: 0. These are executable and path checks, not end-to-end analysis validation.

## Pixi relocation 2026-09-30T03:58:57Z

Workspace: /home/ubuntu/mva-hackathon-2026/tools/pixi
Evidence: /home/ubuntu/mva-hackathon-2026/tools/pixi/logs/relocate_20260930T035851Z_bVO5UX
Old environments removed after verification.

## Environment verification 2026-09-30T03:59:12Z

Log: /home/ubuntu/mva-hackathon-2026/tools/pixi/logs/verify_20260930T035910Z_dWNOM2/verification.txt
Failed checks: 0. These are executable and path checks, not end-to-end analysis validation.

### Verification and cleanup
Recorded UTC: 2026-09-30T04:22:05Z
Verification source: logs/verify_20260930T035910Z_dWNOM2/verification.txt
Failed checks: 0

Verification report SHA256:
e19db164868ce35dc99a88d9827118798607bb4d33ba34cd74a4e2cd63042e5d  logs/verify_20260930T035910Z_dWNOM2/verification.txt

Current configuration SHA256:
112db3f1314dd223eb0f4cb41e78d9733979937c68b1e7bce58b178bf4c6e105  pixi.toml
6db503f8425e981a9a6a09f9a9a5bd54459ed26d1a038bbc5049d69f1943a8a3  pixi.lock

## Environments added 2026-09-30T07:47:11Z

New environments; qc, phenotype, annotate and rank are unchanged in pixi.lock.
- align: bwa 0.7.19, samtools 1.24, mosdepth 0.3.14, whatshap 2.8 (bioconda)
- sv: delly 2.6.0, samtools 1.24 (bioconda)
- manta: python 2.7 runtime for the official Manta 1.6.0 binary in tools/manta-1.6.0.centos6_x86_64
- gatk: imported from tools/gatk-4.7.0.0/gatkcondaenv.yml (official GATK 4.7.0.0 release); only the pip line was rewritten to "gatkpythonpackages @ file://..." for pixi. See logs/add_gatk_env_*/.
Official downloads (checksums verified at download): GATK 4.7.0.0 zip (SHA-256 from GitHub release), Temurin JDK 17.0.20.1 (SHA-256 from Adoptium API), Manta 1.6.0 binary (no checksum published).

Current configuration SHA256:
a1c44264d20de96f714c0ae67179f1794af95b67048be3a08535aa858654cc86  pixi.toml
2a28b988553027e7da6286fa897bdcbcc5be336baf340f11dc0d3b92302d3baf  pixi.lock

## Environment added 2026-09-30T16:54:21Z

- popgen: vcftools 0.1.17, plink 1.90b7.7, plink2 2.0.0a.6.9 (bioconda). Existing environments unchanged in pixi.lock.

Current configuration SHA256:
b9a8b9f24515960fe4066e9827328a2c115e0c303a559e3df36385508f6171a8  pixi.toml
327de82f87f7dc0f00a13ef15d286de5c9543d16672759747bc1888f6f4cf4e8  pixi.lock
