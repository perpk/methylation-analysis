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

test_that("extract_methyl_set returns the expected values", {
    source("R/remove_sex_chromosome_probes.R")

    m_set <- preprocessRaw(rgSet = minfiData::RGsetEx[1:100, 1:2])

    processed_dir <- withr::local_tempdir()
    project_context <- list(paths = list(processed = processed_dir))

    res <- remove_sex_chromosome_probes(
        context = project_context,
        methyl_set = m_set
    )

    expect_true("mismatch_container" %in% names(res))
    expect_s4_class(res$mismatch_container, "ResultsContainer")
    expect_true(file.exists(res$mismatch_container@filename))
    expect_s4_class(res$mismatch_container@object, "MethylSet")
    expect_true(all(!(seqnames(res$mismatch_container@object) %in% c("chrX", "chrY"))))
    expect_equal(nrow(res$mismatch_container@object), nrow(m_set) - sum(seqnames(m_set) %in% c("chrX", "chrY")))
    expect_equal(ncol(res$mismatch_container@object), ncol(m_set))
    expect_equal(colnames(res$mismatch_container@object), colnames(m_set))
    expect_equal(rownames(res$mismatch_container@object), rownames(m_set)[!(seqnames(m_set) %in% c("chrX", "chrY"))])
})
