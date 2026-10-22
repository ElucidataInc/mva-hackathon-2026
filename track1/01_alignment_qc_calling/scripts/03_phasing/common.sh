source "$(dirname "${BASH_SOURCE[0]}")/../01_alignment/common.sh"
case "${2:-}" in
    "") IN=$VCF; PH=$RESULTS/phasing ;;
    provided_hc) IN=$TRACK1/05_vcf_qc_annotation_ranking/results/hc_out/final.high_confidence.vcf.gz; PH=$RESULTS/phasing/provided_hc ;;
    deepvariant_hc) IN=$RESULTS/annotation/deepvariant/final.high_confidence.vcf.gz; PH=$RESULTS/phasing/deepvariant_hc ;;
    *) echo "unknown call set: $2" >&2; exit 1 ;;
esac
CHROMS=$(grep -E '^chr([0-9]+|X)\s' "$REF.fai" | cut -f1)
mkdir -p "$PH/chrom"
