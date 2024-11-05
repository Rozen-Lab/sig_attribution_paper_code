#!/bin/bash
PROJECT_DIR=/home/e0012078/packages/sig_attribution_paper_code
BASH_SCRIPT_DIR=analysis/code/SBS/msa_default/bash/SBS/syn
ABS_BASH_SCRIPT_DIR="$PROJECT_DIR"/$BASH_SCRIPT_DIR
cd $ABS_BASH_SCRIPT_DIR

bash run_Breast-AdenoCA_non_msi_SBS.sh
bash run_ColoRect-AdenoCA_msi_SBS.sh
bash run_ColoRect-AdenoCA_non_msi_SBS.sh
bash run_Eso-AdenoCA_non_msi_SBS.sh
bash run_Kidney-RCC_msi_SBS.sh
bash run_Kidney-RCC_non_msi_SBS.sh
bash run_Liver-HCC_msi_SBS.sh
bash run_Liver-HCC_non_msi_SBS.sh
bash run_Lung-AdenoCA_non_msi_SBS.sh
bash run_Ovary-AdenoCA_msi_SBS.sh
bash run_Ovary-AdenoCA_non_msi_SBS.sh
bash run_Skin-Melanoma_msi_SBS.sh
bash run_Skin-Melanoma_non_msi_SBS.sh
bash run_Stomach-AdenoCA_msi_SBS.sh
bash run_Stomach-AdenoCA_non_msi_SBS.sh
