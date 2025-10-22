library(excode)
library(testthat)
library(dplyr)
library(lubridate)
suppressWarnings(library(surveillance))

data(mort_df_germany)
data(shadar_df)


set.seed(123456)

# Poisson x Harmonic x Linear x 2 states
res_har_pois <- run_excode(surv_ts = shadar_df,
            timepoints = 295,
            distribution = "Poisson",
            states = 2,
            periodic_model = "Harmonic",
            time_trend = "Linear",
            set_baseline_state = TRUE)
sum_har_pois <- summary(res_har_pois)
expected_result <- data.frame(posterior0=0.86692061802787512281298631933168508112430572509765625,
                              posterior1=0.133079381972124988209316143183968961238861083984375,
                              pval=0.188346850907013363407571659990935586392879486083984375,
                              zscore=1.1411988070776402093287060779402963817119598388671875,
                              date=as_date(13388),
                              timepoint=295,
                              observed=6,
                              mu0=3.82631505801585003467835122137330472469329833984375,
                              mu1=12.4394715096714509883213395369239151477813720703125,
                              BIC=1149.24366997181959959561936557292938232421875,
                              AIC=1066.033945491270060301758348941802978515625)
expect_equal(sum_har_pois, expected_result)

# Poisson x Custom x 2 states
cov_data <- res_har_pois@emission@excode_formula@data
res_custom_pois <- run_excode(surv_ts = shadar_df[35:295,1:2],
                           timepoints = 261,
                           covariate_df = cov_data,
                           distribution = "Poisson",
                           states = 2,
                           periodic_model = "Custom",
                           time_trend = "None",
                           set_baseline_state = TRUE)
sum_custom_pois <- summary(res_custom_pois)
sum_custom_pois$timepoint <- 295
expect_equal(sum_custom_pois, sum_har_pois)

# NegBinom x FarringtonNoufaily x None x 2 states
res_fn_nb <- run_excode(surv_ts = shadar_df,
                           timepoints = 295,
                           distribution = "NegBinom",
                           states = 2,
                           periodic_model = "FarringtonNoufaily",
                           time_trend = "None",
                           set_baseline_state = TRUE) 
sum_fn_nb <- summary(res_fn_nb)
expected_result <- data.frame(posterior0=0.94066750386438069408967521667364053428173065185546875,
                              posterior1=0.059332496135619326727006495048044598661363124847412109375,
                              pval=0.22816206357078261390824991394765675067901611328125,
                              zscore=0.82602018894115214475704078722628764808177947998046875,
                              date=as_date(13388),
                              timepoint=295,
                              observed=6,
                              mu0=4.08391097281145132313895373954437673091888427734375,
                              mu1=13.809012170386353091089404188096523284912109375,
                              size=6770.80316358002528431825339794158935546875,
                              BIC=1203.57776920385322227957658469676971435546875,
                              AIC=1074.72236216208739278954453766345977783203125)
expect_equal(sum_fn_nb, expected_result)

# NegBinom x Harmonic x Spline2 x 3 states
res_har_nb_2 <- run_excode(surv_ts = mort_df_germany,
                         timepoints = 325,
                         distribution = "NegBinom",
                         states = 3,
                         periodic_model = "Harmonic",
                         time_trend = "Spline2") 
sum_har_nb_2 <- summary(res_har_nb_2)
expected_result <- data.frame(posterior0=0.8841155288357229746765142408548854291439056396484375,
                              posterior1=0.11588446366601724835110331923715420998632907867431640625,
                              posterior2=0.0000000074982598782866864382812277511902721016667783260345458984375,
                              pval=0.05204923224775277745823842678873916156589984893798828125,
                              zscore=1.9959303424259022818887387984432280063629150390625,
                              date=as_date(19071),
                              timepoint=325,
                              observed=21357,
                              mu0=20342.0921900660541723482310771942138671875,
                              mu1=22228.27810775602483772672712802886962890625,
                              mu2=25566.80308253416660591028630733489990234375,
                              size=1144.232172671881244241376407444477081298828125,
                              BIC=4251.2924495814104375313036143779754638671875,
                              AIC=4145.3891638788927593850530683994293212890625)
expect_equal(sum_har_nb_2, expected_result)


if (FALSE) {
  result_to_code <- function(result) {
    paste0("expected_result <- data.frame(", 
           paste(lapply(names(result), 
                        function(x) paste0(x, "=", ifelse(x=="date", "as_date(", ""), 
                                           ifelse(is.numeric(result[[x]]), 
                                                  gsub("\\.*0+$", "", formatC(result[[x]], digits = 300, format='f')),
                                                  result[[x]]),
                                           ifelse(x=="date", ")", ""))),
                 collapse=",\n"), ")") %>% cat()
  }
  result_to_code(sum_har_pois)
  result_to_code(sum_fn_nb)
  result_to_code(sum_har_nb_2)
}


