#!/usr/bin/bash
set -eux

wget 'https://data.broadinstitute.org/gsea-msigdb/gsea/software/desktop/4.2/GSEA_4.2.1.zip' -O "GSEA_4.2.1.zip"

unzip -u "GSEA_4.2.1.zip"

cp -r "GSEA_4.2.1" "$CONDA_PREFIX/share"

ln -fs "$CONDA_PREFIX/share/GSEA_4.2.1/gsea-cli.sh" "$CONDA_PREFIX/bin"

rm -rf "GSEA_4.2.1" "GSEA_4.2.1.zip"
