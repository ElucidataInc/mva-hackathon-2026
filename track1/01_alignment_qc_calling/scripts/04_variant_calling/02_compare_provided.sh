#!/usr/bin/env bash
source "$(dirname "$0")/../01_alignment/common.sh"
eval "$("$PIXI" shell-hook --manifest-path "$PIXI_MANIFEST" --frozen --no-install -e qc)"
DV=$RESULTS/deepvariant
C=$DV/compare
mkdir -p "$C"

prep() {
    bcftools view -f PASS -Ou "$1" | bcftools norm -m -any -f "$REF" -Ou |
        bcftools view -e 'ALT="*"' -Oz -o "$2"
    bcftools index -t "$2"
}

awk '{c = $1; sub(/^chr/, "", c); print c "\t" $1}' "$REF.fai" > "$C/chr_map.txt"
bcftools annotate --rename-chrs "$C/chr_map.txt" -Ou "$VCF" | prep - "$C/provided.norm.vcf.gz"
prep "$DV/$SAMPLE.deepvariant.vcf.gz" "$C/deepvariant.norm.vcf.gz"
bcftools isec -p "$C/isec" "$C/provided.norm.vcf.gz" "$C/deepvariant.norm.vcf.gz"

{
    printf 'set\tSNV\tindel\n'
    for s in 0000:provided_only 0001:deepvariant_only 0002:shared; do
        f=$C/isec/${s%%:*}.vcf
        printf '%s\t%s\t%s\n' "${s#*:}" "$(bcftools view -H -v snps "$f" | wc -l)" "$(bcftools view -H -v indels "$f" | wc -l)"
    done
    paste <(bcftools query -f '[%GT]\n' "$C/isec/0002.vcf") <(bcftools query -f '[%GT]\n' "$C/isec/0003.vcf") |
        awk '{for (i = 1; i <= 2; i++) {split($i, a, /[\/|]/); $i = a[1] <= a[2] ? a[1] "/" a[2] : a[2] "/" a[1]}; n++; m += $1 == $2}
             END {printf "shared_gt_concordance\t%.4f\n", m / n}'
    bcftools view -f PASS -m2 -M2 -v snps -i 'GT="het"' "$DV/$SAMPLE.deepvariant.vcf.gz" | bcftools query -f '[%VAF]\n' |
        sort -g | awk '{v[NR] = $1; b += $1 >= 0.3 && $1 <= 0.7} END {printf "dv_het_snv_vaf_median\t%.3f\ndv_het_snv_vaf_0.3_0.7\t%.4f\n", v[int(NR / 2) + 1], b / NR}'
} > "$C/comparison_summary.tsv"
