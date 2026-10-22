#!/usr/bin/env bash
# Module 05 – sample QC: vcftools (missingness, het/F, HWE, Ts/Tv, depth), bcftools +guess-ploidy,
# +check-ploidy, roh, +mendelian2, somalier (relatedness/sex/ancestry), peddy, then sample_qc.py
# (KING kinship, PCA, sex check, contamination proxies, ROH summary, outlier flags).

mod_sampleqc() {
    local d="${OUTDIR}/05_sample_qc"; mkdir -p "$d"
    local ns; ns=$(n_samples)
    if [[ "${ns:-0}" -eq 0 ]]; then mnote "sites-only VCF: no sample QC"; return 3; fi

    # ---------------- vcftools
    if have vcftools; then
        mkdir -p "$d/vcftools"
        local vin=(--gzvcf "$VCF_IN"); [[ "$VCF_IN" != *.gz ]] && vin=(--vcf "$VCF_IN")
        for opt in --missing-indv --missing-site --het --hardy --TsTv-summary --depth --site-mean-depth; do
            vcftools "${vin[@]}" "$opt" --out "$d/vcftools/vcftools" >>"$d/vcftools/vcftools.log" 2>&1 || \
                mnote "vcftools $opt failed"
        done
    else
        mnote "vcftools not found"
    fi

    if have bcftools; then
        local X Y g; X=$(contig_name X); Y=$(contig_name Y); g=$(genome_tag)
        # ---------------- sex: +guess-ploidy (GT likelihood of haploid vs diploid X)
        local gp=(-t GT)
        bcf_header | grep -q '^##FORMAT=<ID=PL' && gp=(-t PL)
        if [[ -n "$g" ]]; then gp+=(--genome "$g"); else gp+=(-r "$X"); fi
        if bcftools +guess-ploidy ${gp[@]+"${gp[@]}"} "$VCF_IN" >"$d/guess_ploidy.raw" 2>"$d/guess_ploidy.err"; then
            # default output is "<sample> <M|F>" (newer versions may add a '# [1]...' header line)
            { printf 'sample\tpredicted_sex\n'; awk -F'\t' '!/^#/{ if ($1=="SEX") print $2"\t"$3; else print $1"\t"$2 }' "$d/guess_ploidy.raw"; } \
                >"$d/guess_ploidy.tsv"
        else
            mnote "+guess-ploidy failed"
        fi
        bcftools +check-ploidy "$VCF_IN" >"$d/check_ploidy.txt" 2>&1 || true

        # ---------------- runs of homozygosity (genotype-based)
        local afopt=(--AF-dflt 0.4)
        [[ "$ns" -ge 10 ]] && afopt=(--estimate-AF -)
        bcftools roh --threads "$THREADS" -G30 ${afopt[@]+"${afopt[@]}"} -O r -o "$d/roh.txt" "$VCF_IN" >/dev/null 2>"$d/roh.err" || \
            mnote "bcftools roh failed – fallback ROH scan used"

        # ---------------- Mendelian errors
        if [[ -n "$PED" ]]; then
            local rules=()
            case "$(genome_tag)" in hg38|b38) rules=(--rules GRCh38) ;; hg19|b37) rules=(--rules GRCh37) ;; esac
            bcftools +mendelian2 "$VCF_IN" -P "$PED" -m c ${rules[@]+"${rules[@]}"} >"$d/mendelian2.txt" 2>"$d/mendelian2.err" || {
                awk '$3!="0" && $4!="0"{print $4","$3","$2}' "$PED" >"$d/trios.txt"
                bcftools +mendelian "$VCF_IN" -T "$d/trios.txt" -c >"$d/mendelian2.txt" 2>>"$d/mendelian2.err" || \
                    mnote "+mendelian2/+mendelian failed"; }
        fi
    else
        mnote "bcftools not found: no guess-ploidy/roh/mendelian2"
    fi

    # ---------------- somalier: relatedness, sex and ancestry projected onto 1000G
    if have somalier && [[ -n "$SOMALIER_SITES" && -n "$REF" ]]; then
        mkdir -p "$d/somalier/extracted"
        somalier extract -d "$d/somalier/extracted" --sites "$SOMALIER_SITES" -f "$REF" "$VCF_IN" >"$d/somalier/extract.log" 2>&1 && \
        somalier relate ${PED:+--ped "$PED"} -o "$d/somalier/somalier" "$d/somalier/extracted/"*.somalier >"$d/somalier/relate.log" 2>&1 || \
            mnote "somalier relate failed"
        if [[ -n "$SOMALIER_LABELS" && -n "$SOMALIER_REF_DIR" ]]; then
            somalier ancestry --labels "$SOMALIER_LABELS" "$SOMALIER_REF_DIR"/*.somalier ++ "$d/somalier/extracted/"*.somalier \
                -o "$d/somalier/somalier-ancestry" >"$d/somalier/ancestry.log" 2>&1 || mnote "somalier ancestry failed"
        fi
    else
        mnote "somalier skipped (needs somalier, --somalier-sites, --ref)"
    fi

    # ---------------- peddy: sex / het / pedigree checks
    if have peddy && [[ -n "$PED" ]]; then
        mkdir -p "$d/peddy"
        local site_opt=(); [[ "$(genome_tag)" == *38 ]] && site_opt=(--sites hg38)
        peddy -p "$THREADS" ${site_opt[@]+"${site_opt[@]}"} --prefix "$d/peddy/peddy" "$VCF_IN" "$PED" >"$d/peddy/peddy.log" 2>&1 || \
            mnote "peddy failed"
    else
        mnote "peddy skipped (needs peddy + --ped)"
    fi

    # ---------------- Python sample QC on profiler outputs
    [[ -f "${OUTDIR}/04_profile/profile.json" ]] || { mnote "profile missing"; return 1; }
    local sexsrc="${SEX_FILE:-$PED}"
    "$PYTHON" "${SCRIPTS}/sample_qc.py" --profile-dir "${OUTDIR}/04_profile" --out "$d" \
        ${sexsrc:+--sex-file "$sexsrc"} $( [[ -s "$d/roh.txt" ]] && echo --roh "$d/roh.txt" )
}
