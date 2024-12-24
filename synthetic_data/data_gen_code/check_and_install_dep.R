if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

if (!requireNamespace("magrittr", quietly = TRUE)) {
  install.packages("magrittr")
}

if (!requireNamespace("ICAMS", quietly = TRUE)) {
  remotes::install_github("steverozen/ICAMS", ref = "v3.0.8-branch")
}

if (!requireNamespace("cosmicsig", quietly = TRUE)) {
  remotes::install_github(repo = "Rozen-Lab/cosmicsig", ref = "v1.2.0-branch")
}

if (!requireNamespace("PCAWG7", quietly = TRUE)) {
  remotes::install_github(repo = "steverozen/PCAWG7", ref = "v0.1.3-branch")
}

if (!requireNamespace("SynSigGen", quietly = TRUE) ||
    packageVersion("SynSigGen") < "1.2.2") {
  remotes::install_github(repo = "steverozen/SynSigGen",ref = "1.2.2-branch")
}

if (!requireNamespace("mSigAct", quietly = TRUE)) {
  remotes::install_github(repo = "steverozen/mSigAct", ref = "v3.0.1-branch")
}

if (!requireNamespace("MutationalPatterns", quietly = TRUE)) {
  install.packages(pkgs = "https://bioconductor.org/packages/3.16/bioc/src/contrib/MutationalPatterns_3.8.1.tar.gz", 
                   repos = NULL)
}

# Make sure to have Rtools installed on your Windows machine before proceeding
# with the installation of NNLM, which is a dependency package for
# signature.tools.lib. See https://cran.r-project.org/bin/windows/Rtools/ for
# more details

# To accommodate the download duration of the four BSgenome dependency packages
# required for signature.tools.lib—BSgenome.Cfamiliaris.UCSC.canFam3,
# BSgenome.Mmusculus.UCSC.mm10, BSgenome.Hsapiens.1000genomes.hs37d5, and
# BSgenome.Hsapiens.UCSC.hg38—which takes more than 60 seconds, you may need to
# adjust the download timeout option in R.
if (getOption('timeout') == 60) {
  options(timeout = 1000)
}
if (!requireNamespace("signature.tools.lib", quietly = TRUE)) {
  remotes::install_github(repo = "Nik-Zainal-Group/signature.tools.lib", ref = "v2.3.0")
}


