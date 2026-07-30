
<!-- README.md is generated from README.Rmd. Please edit that file and run
     devtools::build_readme() (or rmarkdown::render("README.Rmd")) -->

# ReproMeta

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

***ReproMeta*** provides a unified R interface for assessing replication
success across several statistical interpretations of reproducibility,
i.e., reproducibility metrics. *ReproMeta* allows meta-researchers to
compute a whole set of reproducibility metrics systematically across the
same set of original- replication study pairs within a single and
consistent workflow, without having to rely on multiple packages and
metric-specific implementations.

> Why does a unified interface matter? Metrics disagree in practice due
> to their underlying definitions, assumptions and limitations (see DOI:
> <https://doi.org/10.1098/rsos.242076>). In Mita et al. (2016 – to
> appear), we systematically investigated the mathematical deifnition,
> assumptions and limitation of known reproducibility metrics, and
> established a unified mathematical notation to allow comparison.
> *ReproMeta* provides the implementation of the metrics in R. This
> directly helps applied researchers and meta-researchers to compute
> several of metrics side by side on the same study pair.

## Installation

ReproMeta is hosted on the UZH GitLab instance. To install the package,
copy the following into your console:

``` r
# install.packages("remotes")
remotes::install_gitlab("crsuzh/reprometa", host = "gitlab.uzh.ch")
```

## Usage – quick start

``` r
# After installation
library(ReproMeta)

# Prepare data
replication_data <- prepare_effect_sizes(r_o = 0.35, n_o = 100,
                                         r_r = 0.21, n_r = 120, 
                                         effectsize = "COR")

replication_results <- assess_replication_outcome(replication_data,
  metric = c("significance_r", "pScept", "BFr")
)

resqual <- get_replication_results(replication_results, type = "qualitative")
resquan <- get_replication_results(replication_results, type = "quantitative")
```

The qualitative results can be read of in the table below
<table>

<thead>

<tr>

<th style="text-align:right;">

Study
</th>

<th style="text-align:left;">

Metric
</th>

<th style="text-align:left;">

Outcome
</th>

<th style="text-align:left;">

Details
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1
</td>

<td style="text-align:left;">

Significance criterion
</td>

<td style="text-align:left;">

Replication success
</td>

<td style="text-align:left;">

replication effect is significant and same direction of original
</td>

</tr>

<tr>

<td style="text-align:right;">

1
</td>

<td style="text-align:left;">

Sceptical p-value (controlled)
</td>

<td style="text-align:left;">

Replication success
</td>

<td style="text-align:left;">

Evidence for replication success
</td>

</tr>

<tr>

<td style="text-align:right;">

1
</td>

<td style="text-align:left;">

Replication Bayes factor
</td>

<td style="text-align:left;">

Replication success
</td>

<td style="text-align:left;">

Evidence for replication success
</td>

</tr>

</tbody>

</table>

## Contributing and support

For bug reports and feature requests, open an issue at
<https://gitlab.uzh.ch/crsuzh/reprometa/-/issues> or send an email to
XXXX.

## Related work

- [ReplicationSuccess](https://CRAN.R-project.org/package=ReplicationSuccess)
  – design and analysis of replication studies
- [FReD](https://forrt.org/fred/) – provides various interfaces to the
  FORRT Replication Database
