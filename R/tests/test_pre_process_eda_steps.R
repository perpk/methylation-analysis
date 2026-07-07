if (!require(testthat)) {
    install.packages("testthat")
}
if (!require(minfiData)) {
    BiocManager::install("minfiData")
}

library(testthat)
library(withr)
library(minfiData)
library(minfi)

test_that("remove_sex_chromosome_probes returns the expected values", {
    source("R/remove_sex_chromosome_probes.R")
    source("R/results_container.R")

    m_set <- preprocessRaw(rgSet = minfiData::RGsetEx[1:100, 1:2])
    m_set <- mapToGenome(m_set)

    processed_dir <- withr::local_tempdir()
    project_context <- list(paths = list(processed = processed_dir))

    res <- remove_sex_chromosome_probes(
        context = project_context,
        methyl_set = m_set
    )

    expect_true("mismatch_container" %in% names(res))
    expect_s4_class(res$mismatch_container, "ResultsContainer")
    expect_s4_class(res$mismatch_container@object, "GenomicMethylSet")
    expect_true(all(!(seqnames(res$mismatch_container@object) %in% c("chrX", "chrY"))))
    expect_equal(nrow(res$mismatch_container@object), nrow(m_set) - sum(seqnames(m_set) %in% c("chrX", "chrY")))
    expect_equal(ncol(res$mismatch_container@object), ncol(m_set))
    expect_equal(colnames(res$mismatch_container@object), colnames(m_set))
})

test_that("remove_cross_reactive_probes returns the expected values", {
    if (!("devtools" %in% installed.packages())) {
        install.packages("devtools")
    }
    devtools::install_github("markgene/maxprobes")
    source("R/remove_cross_reactive_probes.R")
    source("R/results_container.R")
    source("R/progress_mgr.R")

    m_set <- preprocessRaw(rgSet = minfiData::RGsetEx[1:100, 1:2])
    m_set <- mapToGenome(m_set)

    processed_dir <- withr::local_tempdir()
    project_context <- list(
        paths = list(processed = processed_dir),
        platform = "450k"
    )

    res <- remove_cross_reactive_probes(
        context = project_context,
        methyl_set = m_set
    )

    print(nrow(res$methyl_set_cross_reactive_clean_container@object))
    print(nrow(m_set))

    expect_true("methyl_set_cross_reactive_clean_container" %in% names(res))
    expect_s4_class(res$methyl_set_cross_reactive_clean_container, "ResultsContainer")
    expect_s4_class(res$methyl_set_cross_reactive_clean_container@object, "GenomicMethylSet")
    expect_true(all(!(rownames(res$methyl_set_cross_reactive_clean_container@object) %in% c("cg00000029", "cg00000108"))))
    expect_equal(nrow(res$methyl_set_cross_reactive_clean_container@object), nrow(m_set) - 5)
    expect_equal(ncol(res$methyl_set_cross_reactive_clean_container@object), ncol(m_set))
    expect_equal(colnames(res$methyl_set_cross_reactive_clean_container@object), colnames(m_set))
})
