

## mpd data
mpb_non_thin <- read.csv("~/Downloads/mpb_non_thinned.csv")
mpb_pct <- read.csv("~/Downloads/mpb_precommercial_thin.csv")

## fir engraver
fe_non_thin <- read.csv("~/Downloads/firengraver_non_thinned.csv")
fe_ct <- read.csv("~/Downloads/firengraver_commercial_thin.csv")
fe_salvage <- read.csv("~/Downloads/firengraver_salvage_thin.csv")

mpb_non_thin_n <- nrow(mpb_non_thin)
mpb_pct_n <- nrow(mpb_pct)

mpb_non_thin_mean <- mean(mpb_non_thin$y)
mpd_pct_mean <- mean(mpb_pct$y)

mpb_non_thin_sd <- sd(mpb_non_thin$y)
mpb_pct_sd <- sd(mpb_pct$y)

mpb_non_thin_se <- mpb_non_thin_sd / sqrt(mpb_non_thin_n)
mpb_pct_se <- mpb_pct_sd / sqrt(mpb_pct_n)


fe_non_thin_n <- nrow(fe_non_thin)
fe_ct_n <- nrow(fe_ct)
fe_salvage_n <- nrow(fe_salvage)

fe_non_thin_mean <- mean(fe_non_thin$y)
mpd_ct_mean <- mean(fe_ct$y)
mpd_salvage_mean <- mean(fe_salvage$y)

fe_non_thin_sd <- sd(fe_non_thin$y)
fe_ct_sd <- sd(fe_ct$y)
fe_salvage_sd <- sd(fe_salvage$y)

fe_non_thin_se <- fe_non_thin_sd / sqrt(fe_non_thin_n)
fe_ct_se <- fe_ct_sd / sqrt(fe_ct_n)
fe_salvage_se <- fe_salvage_sd / sqrt(fe_salvage_n)
