#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
mkdir -p "$SVREF" && cd "$SVREF"
for u in https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
         https://zenodo.org/records/11642574/files/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz \
         https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.5.tar.gz; do
    [[ -s ${u##*/} ]] || curl -fL -C - -O "$u"
done
[[ -d Annotations_Human ]] || { gzip -t Annotations_Human_3.5.tar.gz && tar xzf Annotations_Human_3.5.tar.gz; }
