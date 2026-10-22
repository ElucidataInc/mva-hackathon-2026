set -euo pipefail
PROJECT_ROOT=/home/ubuntu/mva-hackathon-2026
PIXI=$PROJECT_ROOT/tools/pixi/0.81.0/bin/pixi
PIXI_MANIFEST=$PROJECT_ROOT/tools/pixi/pixi.toml
REF=$PROJECT_ROOT/references/genome/hs38d1_maskedGRC_v2/GRCh38_no_alt_plus_hs38d1_maskedGRC_exclusions_v2.fna
GATK=$PROJECT_ROOT/tools/gatk-4.7.0.0/gatk
export JAVA_HOME=$PROJECT_ROOT/tools/jdk-17.0.20.1+1
export PATH=$JAVA_HOME/bin:$PATH

run() { local e=$1; shift; "$PIXI" run --manifest-path "$PIXI_MANIFEST" --frozen --no-install -e "$e" "$@"; }
fq() { echo "$FQ_DIR"/*_"$1"_"$2"_001.fastq.gz; }

FQ_DIR=$PROJECT_ROOT/data/fastq
VCF=$(echo "$PROJECT_ROOT"/data/variants/*.vcf.gz)
[[ -f $VCF ]] || { echo "expected one VCF in data/variants" >&2; exit 1; }
SAMPLE=$(run qc bcftools query -l "$VCF")
[[ $SAMPLE =~ ^[^[:space:]]+$ ]] || { echo "expected one sample in $VCF" >&2; exit 1; }
LANES=$(ls "$FQ_DIR" | sed -n 's/.*_\(L[0-9]\{3\}\)_R1_001\.fastq\.gz$/\1/p')
THREADS=$(nproc)
SCRIPTS=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
TRACK1=$(dirname "$(dirname "$SCRIPTS")")
RESULTS=$(dirname "$SCRIPTS")/results
MVA_GENES=$(awk 'NR > 1 && NF {print $1}' "$SCRIPTS/marker_genes.tsv" | paste -sd'|' -)
A=$RESULTS/alignment
QC=$RESULTS/qc
RG_TSV=$A/read_groups.tsv
FINAL_BAM=$A/$SAMPLE.md.bam

RUN_DIR=$(mkdir -p "${1:?run folder required}" && cd "$1" && pwd)
mkdir -p "$A/tmp" "$QC"
