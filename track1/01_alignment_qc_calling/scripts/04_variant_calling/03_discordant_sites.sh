#!/usr/bin/env bash
source "$(dirname "$0")/../01_alignment/common.sh"
eval "$("$PIXI" shell-hook --manifest-path "$PIXI_MANIFEST" --frozen --no-install -e qc)"
C=$RESULTS/deepvariant/compare
I=$C/isec

{
    for s in 0000:provided_only 0001:deepvariant_only; do
        bcftools query -f '%CHROM\n' "$I/${s%%:*}.vcf" |
            awk -v s="${s#*:}" -v OFS='\t' '{c[$1 ~ /^chr([0-9]+|X|Y)$/ ? "chr1-22XY" : "unplaced_decoy"]++} END {for (k in c) print s, "contigs", k, c[k]}'
    done
    bcftools view -H -T "$I/0000.vcf" "$RESULTS/deepvariant/$SAMPLE.deepvariant.vcf.gz" | cut -f7 | sort | uniq -c |
        awk -v OFS='\t' '{print "provided_only", "deepvariant_filter", $2, $1}'
    bcftools annotate --rename-chrs "$C/chr_map.txt" -Ou "$VCF" | bcftools view -H -T "$I/0001.vcf" | cut -f7 | sort | uniq -c |
        awk -v OFS='\t' '{print "deepvariant_only", "provided_filter", $2, $1}'
} > "$C/discordant_summary.tsv"
