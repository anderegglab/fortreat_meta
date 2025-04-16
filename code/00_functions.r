
## define functions
overall_mean <- function(mean, n) {
  return(sum(mean * n) / sum(n))
}

overall_sd <- function(mean, n, sd) {
  return({
    sqrt((1 / (sum(n) - 1)) *
         (sum((n - 1) * sd^2) + sum(n * (mean - overall_mean(mean, n)))))
    })
}

combine_se  <- function(mean, n, se, N) {
  return({
    sqrt(((N[1] * se[1]^2 + N[2] * se[2]^2) / (N[1] + N[2])) +
         ((n[1] * n[1] * (mean[1] - mean[2])^2) / ((n[1] + n[2]) * (N[1] + N[2]))))
  })
}

overall_se <- function(mean, n, se) {
  while(TRUE) {
    N <- n^2 - n
    ov_se <- combine_se(mean[1:2], n[1:2], se[1:2], N[1:2])
    if (length(mean) == 2) {
      break
    }
    mean <- c(mean(mean[1:2]), mean[3:length(mean)])
    n <- c(sum(n[1:2]), n[3:length(n)])
    se <- c(ov_se, se[3:length(se)])
  }
  return(ov_se)
}

group_data <- function(data, remove_zeros = FALSE, add_constant = FALSE) {

  if (remove_zeros) {
    data <- data[data$treatment_se != 0 & data$control_se != 0,]
  } else if(add_constant) {
    data[data$se_treatment == 0 & !is.na(data$se_treatment), c("se_treatment", "sd_treatment")] <- 1e-10
    data[data$se_control == 0 & !is.na(data$se_control), c("se_control", "sd_control")] <- 1e-10
  }

  data_grouped <- data[0,]

  ## loop through studies
  for (id in unique(data$studyID)) {

    sdata <- data[data$studyID == id,] ## pull all data from ID

    if (!is.na(sdata[1,"grouping_flags"])) {

      for (grp in unique(sdata$grouping_flags)) {

        gdata <- sdata[sdata$grouping_flags == grp,] ## pull data with same group
        ndata <- sdata[1,]

        ndata[,c("mean_treatment", "n_treatment", "sd_treatment", "se_treatment")] <-
         c(overall_mean(gdata$mean_treatment, gdata$n_treatment),
           sum(gdata$n_treatment),
           overall_sd(gdata$mean_treatment, gdata$n_treatment, gdata$sd_treatment),
           overall_se(gdata$mean_treatment, gdata$n_treatment, gdata$se_treatment))

        ndata[,c("mean_control", "n_control", "sd_control", "se_control")] <-
          c(overall_mean(gdata$mean_control, gdata$n_control),
            sum(gdata$n_control),
            overall_sd(gdata$mean_control, gdata$n_control, gdata$sd_control),
            overall_se(gdata$mean_control, gdata$n_control, gdata$se_control))

        data_grouped <- rbind(data_grouped, ndata)
      }

    } else {

      data_grouped <- rbind(data_grouped, sdata)

    }
  }
  return(data_grouped)
}


## fn to calculate log response ratio
lrr <- function(mean_1, mean_2) {
  return(log(mean_1 / mean_2))
}

## fn to calculate standard error of log response ratio
lrr_se <- function(mean_1, mean_2, se_1, se_2) {
  return({
    sqrt((se_1^2 / mean_1^2) + (se_2^2  / mean_2^2))
  })
}

table_gen <- function(pool, outfile = "default.csv") {

  summary <- data.frame(pool)
  output <- data.frame(variable = as.character(summary$term),
                       estimate = summary$estimate,
                       se = summary$std.error,
                       ci_low = summary$estimate - 1.97 * summary$std.error,
                       ci_high = summary$estimate + 1.97 * summary$std.error,
                       p_value = summary$p.value)
  for (i in 1:nrow(output)) {
    output[i, "variable"] <- gsub("disturbance_type", "", output[i, "variable"])
    output[i, "variable"] <- gsub("_bin", "", output[i, "variable"])
  }
  write.csv(output, paste0("tables/", outfile))

}
