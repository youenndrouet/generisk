#' summary statistics table of the number of at-risk and the number of cancers by age and genotype from a generisk object
#'
#' @param x an object returned by the generisk function
#'
#' @return to be completed
#' @export
#'
#' @examples
#' #to be completed
#'
#' @import dplyr
#'
statdesc_table <- function(x){

  ndis  <- length(x$params$Ft.pop)
  nloci <- length(x$params$fA)
  cc_event <- seq((6+nloci+2)-1, (6+nloci+2*ndis)-1, 2) # columns coding for diseases events 0/1
  cc_ages  <- seq((6+nloci+2), (6+nloci+2*ndis), 2)  # columns coding for diseases ages at diagnosis or last news

  OUT <- NULL

  outcomes <- names(x$params$FIT.pars)
  outcomes.Ftpop <- names(x$params$Ft.pop)
  age_var_names <- colnames(x$DATA_generisk)[cc_ages]
  event_var_names <- colnames(x$DATA_generisk)[cc_event]

  for (dis in seq_along(1:ndis)){

    age_var <- age_var_names[dis]
    event_var <- event_var_names[dis]
    outcome <- outcomes[dis]
    outcome.Ftpop <- outcomes.Ftpop[dis]

    out <- statdesc(x, age_var = age_var, event_var = event_var, outcome = outcome.Ftpop) %>%
      mutate("OUTCOME" = outcome)

    OUT <- bind_rows(OUT, out)

  }

  OUT <- OUT %>% relocate(.data$OUTCOME)

  return(OUT)

}


statdesc <- function(x, age_var, event_var, outcome){

  out <- x$DATA_generisk %>%
    mutate(

      "SEX" = x$DATA_generisk[,3],
      "PROBAND_FLAG" = x$DATA_generisk[,6],
      "GENO" = x$DATA_generisk[,7],

      "SEX" = if_else(.data$SEX == 1, "Male", "Female"),
      "AGEVAR" = !!as.name(age_var),
      "EVENTVAR" = !!as.name(event_var),
      "TYP_IND" = case_when(
        PROBAND_FLAG == 1 ~ "Proband (Mt)",
        PROBAND_FLAG == 0 & GENO == 1 ~ "Relative (Mt)",
        PROBAND_FLAG == 0 & GENO == 0 ~ "Relative (wt)",
        PROBAND_FLAG == 0 & GENO == 4 ~ "Untested"
      )) %>%
    group_by(.data$TYP_IND, .data$SEX) %>%
    summarize(ukn_age = sum(.data$AGEVAR == 0),
              "at_risk_1y"  = sum(.data$AGEVAR >= 1),
              "at_risk_10y" = sum(.data$AGEVAR >= 10),
              "at_risk_20y" = sum(.data$AGEVAR >= 20),
              "at_risk_30y" = sum(.data$AGEVAR >= 30),
              "at_risk_40y" = sum(.data$AGEVAR >= 40),
              "at_risk_50y" = sum(.data$AGEVAR >= 50),
              "at_risk_60y" = sum(.data$AGEVAR >= 60),
              "at_risk_70y" = sum(.data$AGEVAR >= 70),
              "at_risk_80y" = sum(.data$AGEVAR >= 80),
              "at_risk_90y" = sum(.data$AGEVAR >= 90),
              #n_evt_ukn_age = sum(EVENTVAR == 1 & AGEVAR ==0),
              "n_evt_1_10y" = sum(.data$EVENTVAR == 1 & .data$AGEVAR > 0 & .data$AGEVAR <= 10),
              "n_evt_11_20y" = sum(.data$EVENTVAR == 1 & .data$AGEVAR > 10 & .data$AGEVAR <= 20),
              "n_evt_21_30y" = sum(.data$EVENTVAR == 1 & .data$AGEVAR > 20 & .data$AGEVAR <= 30),
              "n_evt_31_40y" = sum(.data$EVENTVAR == 1 & .data$AGEVAR > 30 & .data$AGEVAR <= 40),
              "n_evt_41_50y" = sum(.data$EVENTVAR == 1 & .data$AGEVAR > 40 & .data$AGEVAR <= 50),
              "n_evt_51_60y" = sum(.data$EVENTVAR == 1 & .data$AGEVAR > 50 & .data$AGEVAR <= 60),
              "n_evt_61_70y" = sum(.data$EVENTVAR == 1 & .data$AGEVAR > 60 & .data$AGEVAR <= 70),
              "n_evt_71_80y" = sum(.data$EVENTVAR == 1 & .data$AGEVAR > 70 & .data$AGEVAR <= 80),
              "n_evt_81_90y" = sum(.data$EVENTVAR == 1 & .data$AGEVAR > 80 & .data$AGEVAR <= 90),
              "n_evt_91_100y" = sum(.data$EVENTVAR == 1 & .data$AGEVAR > 90 & .data$AGEVAR <= 100),
              '(0-10]'  = ifelse(.data$n_evt_1_10y == 0, paste0(.data$at_risk_1y), paste0(.data$at_risk_1y, " (", .data$n_evt_1_10y, ")")),
              '(10-20]' = ifelse(.data$n_evt_11_20y == 0, paste0(.data$at_risk_10y), paste0(.data$at_risk_10y, " (", .data$n_evt_11_20y, ")")),
              '(20-30]' = ifelse(.data$n_evt_21_30y == 0, paste0(.data$at_risk_20y), paste0(.data$at_risk_20y, " (", .data$n_evt_21_30y, ")")),
              '(30-40]' = ifelse(.data$n_evt_31_40y == 0, paste0(.data$at_risk_30y), paste0(.data$at_risk_30y, " (", .data$n_evt_31_40y, ")")),
              '(40-50]' = ifelse(.data$n_evt_41_50y == 0, paste0(.data$at_risk_40y), paste0(.data$at_risk_40y, " (", .data$n_evt_41_50y, ")")),
              '(50-60]' = ifelse(.data$n_evt_51_60y == 0, paste0(.data$at_risk_50y), paste0(.data$at_risk_50y, " (", .data$n_evt_51_60y, ")")),
              '(60-70]' = ifelse(.data$n_evt_61_70y == 0, paste0(.data$at_risk_60y), paste0(.data$at_risk_60y, " (", .data$n_evt_61_70y, ")")),
              '(70-80]' = ifelse(.data$n_evt_71_80y == 0, paste0(.data$at_risk_70y), paste0(.data$at_risk_70y, " (", .data$n_evt_71_80y, ")")),
              '(80-90]' = ifelse(.data$n_evt_81_90y == 0, paste0(.data$at_risk_80y), paste0(.data$at_risk_80y, " (", .data$n_evt_81_90y, ")")),
              '(90-100]'= ifelse(.data$n_evt_91_100y == 0, paste0(.data$at_risk_90y), paste0(.data$at_risk_90y, " (", .data$n_evt_91_100y, ")")),
              .groups = "keep"
    )

  if(x$params$Ft.pop[[outcome]][120,"m"] == 0){
    out <- out %>% filter(.data$SEX == "Female")
  }else{
    if(x$params$Ft.pop[[outcome]][120,"f"] == 0){
      out <- out %>% filter(.data$SEX == "Male")
    }
  }

  return(out)

}
