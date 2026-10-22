#!/usr/bin/env bash
set -euo pipefail
umask 077
root="$HOME/mva-hackathon-2026"
cd "$root"
pixi_bin="$root/tools/pixi/0.81.0/bin/pixi"
[[ -x "$pixi_bin" && -f "$root/tools/pixi/pixi.lock" ]] || { echo 'Complete installation first'; exit 2; }
mkdir -p tools/pixi/logs
run=$(mktemp -d "$root/tools/pixi/logs/verify_$(date -u +%Y%m%dT%H%M%SZ)_XXXXXX")
cp -- "${BASH_SOURCE[0]}" "$run/verify.sh"
exec > >(tee "$run/verification.txt") 2>&1
failures=0
check() {
  printf '\nCHECK:'; printf ' %q' "$@"; printf '\n'
  if "$@"; then echo 'PASS'; else rc=$?; echo "FAIL: exit $rc"; failures=$((failures+1)); fi
}
check "$pixi_bin" --version
check "$pixi_bin" run --manifest-path "$root/tools/pixi/pixi.toml" --locked -e qc bcftools --version
check "$pixi_bin" run --manifest-path "$root/tools/pixi/pixi.toml" --locked -e qc samtools --version
check "$pixi_bin" run --manifest-path "$root/tools/pixi/pixi.toml" --locked -e qc fastqc --version
check "$pixi_bin" run --manifest-path "$root/tools/pixi/pixi.toml" --locked -e qc multiqc --version
check "$pixi_bin" run --manifest-path "$root/tools/pixi/pixi.toml" --locked -e phenotype python -c 'import pandas, requests, networkx, obonet, pyhpo; from importlib.metadata import version; print({p: version(p) for p in ["pandas", "requests", "networkx", "obonet", "pyhpo"]})'
check "$pixi_bin" run --manifest-path "$root/tools/pixi/pixi.toml" --locked -e annotate vep --help
check "$pixi_bin" run --manifest-path "$root/tools/pixi/pixi.toml" --locked -e rank java -version
check "$pixi_bin" run --manifest-path "$root/tools/pixi/pixi.toml" --locked -e rank java -jar "$root/tools/exomiser-cli-15.1.0/exomiser-cli-15.1.0.jar" --help
check test -d references/vep/116/homo_sapiens/116_GRCh38
check test -s references/vep/116/Homo_sapiens.GRCh38.dna.primary_assembly.fa
check test -s references/vep/116/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
check test -d references/exomiser/2602/2602_hg38
check test -d references/exomiser/2602/2602_phenotype
check sha256sum tools/pixi/pixi.toml tools/pixi/pixi.lock
check df -h .
printf '\nFailed checks: %s\nVerification output: %s/verification.txt\n' "$failures" "$run"
printf '\n## Environment verification %s\n\nLog: %s/verification.txt\nFailed checks: %s. These are executable and path checks, not end-to-end analysis validation.\n' "$(date -u +%FT%TZ)" "$run" "$failures" >> tools/pixi/environment_history.md
(( failures == 0 ))
