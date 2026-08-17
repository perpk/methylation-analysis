rm(list = ls())
gc(full = TRUE)

library(tidyverse)

setClass(
    "DataSources",
    slots = list(
        targets = "ANY",
        m_values = "ANY",
        pca_df = "ANY"
    )
)

data_sources <- list(
    sgpd = new("DataSources",
        targets = "/tmp/test",
        m_values = "/tmp/test",
        pca_df = "/tmp/test"
    ),
)

m_values_sgpd <- readRDS(data_sources$sgpd@m_values)
targets_sgpd <- readRDS(data_sources$sgpd@targets)
pca_df <- readRDS(data_sources$sgpd@pca_df)
