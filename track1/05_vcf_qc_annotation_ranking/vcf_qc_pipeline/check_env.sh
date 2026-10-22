#!/usr/bin/env bash
# Report what is available on this machine for the VCF QC pipeline:
# system resources, environment managers, conda envs, and which tools / VEP caches /
# SnpEff databases exist. Read-only: it installs and changes nothing.
#   bash check_env.sh            (or: bash check_env.sh > env_report.txt)

TOOLS="python3 bcftools bgzip tabix vcftools plot-vcfstats gatk somalier peddy vep vep_install snpEff snpeff hap.py rtg java"

hr() { printf '\n==== %s ====\n' "$1"; }
have() { command -v "$1" >/dev/null 2>&1; }

hr "System"
echo "host:   $(hostname)"
echo "os:     $( (. /etc/os-release 2>/dev/null && echo "$PRETTY_NAME") || uname -sr)"
echo "cpus:   $(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null)"
if have free; then echo "memory: $(free -h | awk '/^Mem:/{print $2" total, "$7" available"}')"; fi
echo "disk (home, current dir, /tmp):"
df -h "$HOME" "$PWD" /tmp 2>/dev/null | awk 'NR==1 || !seen[$0]++' | sed 's/^/  /'

hr "Environment managers"
for m in conda mamba micromamba module docker singularity apptainer; do
    if have "$m"; then echo "  [x] $m  ($(command -v "$m"))"; else echo "  [ ] $m"; fi
done
echo "  active conda env: ${CONDA_DEFAULT_ENV:-none}"
echo "  active virtualenv: ${VIRTUAL_ENV:-none}"

hr "Tools on current PATH"
for t in $TOOLS; do
    if have "$t"; then printf '  [x] %-14s %s\n' "$t" "$(command -v "$t")"; else printf '  [ ] %s\n' "$t"; fi
done
if have python3; then
    python3 - <<'EOF'
import importlib
for m in ("numpy", "matplotlib", "pandas", "hail"):
    try:
        mod = importlib.import_module(m)
        print(f"  [x] python: {m} {getattr(mod, '__version__', '')}")
    except Exception:
        print(f"  [ ] python: {m}")
EOF
fi

CONDA_BIN=""
for c in conda mamba micromamba; do have "$c" && { CONDA_BIN="$c"; break; }; done
if [ -z "$CONDA_BIN" ]; then
    for p in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3" "$HOME/mambaforge" /opt/conda /opt/anaconda3 /opt/miniconda3; do
        [ -x "$p/bin/conda" ] && { CONDA_BIN="$p/bin/conda"; break; }
    done
fi
if [ -n "$CONDA_BIN" ]; then
    hr "Conda environments ($CONDA_BIN)"
    "$CONDA_BIN" env list 2>/dev/null | grep -v '^#' | grep -v '^$' | sed 's/^/  /'
    hr "Tools inside each conda environment"
    "$CONDA_BIN" env list 2>/dev/null | grep -v '^#' | awk 'NF{print $NF}' | while read -r envp; do
        [ -d "$envp/bin" ] || continue
        found=""
        for t in $TOOLS; do [ -x "$envp/bin/$t" ] && found="$found $t"; done
        printf '  %-45s %s\n' "$envp" "${found:- (none of the pipeline tools)}"
    done
fi

if have module; then
    hr "Environment modules matching pipeline tools"
    module avail 2>&1 | tr -s ' ' '\n' | grep -iE 'bcftools|htslib|samtools|vcftools|gatk|vep|snpeff|somalier|peddy|rtg|hap' | sort -u | sed 's/^/  /'
fi

hr "VEP caches"
for d in "$HOME/.vep" /data/vep /opt/vep /usr/share/vep ${VEP_CACHE:-}; do
    [ -d "$d" ] || continue
    echo "  $d:"
    ls -d "$d"/*/*_GRCh* 2>/dev/null | sed 's/^/    /' || echo "    (no species/version folders)"
done
[ -d "$HOME/.vep" ] || echo "  no ~/.vep (set VEP_CACHE=/path before running to check another location)"

hr "SnpEff databases"
for s in snpEff snpeff; do
    have "$s" || continue
    base=$(dirname "$(readlink -f "$(command -v "$s")" 2>/dev/null || command -v "$s")")
    for d in "$base"/../share/snpeff*/data "$base"/data "$HOME/snpEff/data"; do
        [ -d "$d" ] && { echo "  $d:"; ls "$d" | sed 's/^/    /'; }
    done
    break
done
have snpEff || have snpeff || echo "  snpEff not on PATH"

hr "Summary for the VCF QC pipeline"
ok() { have "$1" && echo "yes" || echo "no"; }
echo "  core (python3+numpy+matplotlib): $(python3 -c 'import numpy, matplotlib' 2>/dev/null && echo yes || echo no)"
echo "  bcftools: $(ok bcftools)   vcftools: $(ok vcftools)   gatk: $(ok gatk)"
echo "  vep: $(ok vep)   snpEff: $( (have snpEff || have snpeff) && echo yes || echo no)   somalier: $(ok somalier)   peddy: $(ok peddy)"
echo "  (tools found only inside a conda env need 'conda activate <env>' first)"
