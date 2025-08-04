# Predefined list of available metrics
.available_nca_metrics <- c(
  "N_Samples", "Rsq", "Rsq_adjusted", "Corr_XY", "No_points_lambda_z", "Lambda_z", "Lambda_z_intercept", "Lambda_z_lower", "Lambda_z_upper", "Halflife_Lambda_z", "Span", "Tmax", "Cmax", "Cmax_D", "Tlast", "Clast", "Clast_pred", "AUClast", "AUCall", "AUCINF_obs", "AUC_Extrap_obs", "AUCINF_pred", "AUC_Extrap_pred", "AUMClast", "AUMCINF_obs", "AUMC_Extrap_obs", "AUMCINF_pred", "AUMC_Extrap_pred", "MRTlast", "MRTINF_obs", "MRTINF_pred", "AUC7day", "AUC28day", "AUC3month", "Vz_F_obs", "Cl_F_obs", "Vz_F_pred", "Cl_F_pred"
)

#' Validate requested metrics
#' @keywords internal
validate_metrics <- function(requested_metrics, available_metrics) {
  if ("all" %in% requested_metrics) {
    return(available_metrics)
  }

  invalid <- setdiff(requested_metrics, available_metrics)
  if (length(invalid) > 0) {
    stop("Invalid metric(s): ", paste(invalid, collapse = ", "))
  }

  return(requested_metrics)
}

#' Calculate NCA metrics for T-cell persistence data
#'
#' @param time Numeric vector of time points
#' @param conc Numeric vector of concentrations (e.g., pg/mL, VCN)
#' @param dose Dosage
#' @param metrics Character vector of metrics to calculate, or "all"
#' @return Named list of NCA results
#' @export
calculate_nca_tcell <- function(time, conc, dose, metrics = "all") {
  metrics_to_calc <- validate_metrics(metrics, .available_nca_metrics)

  # data preprocessing
  tmp <- dataQc(time = time, conc = conc, sampletype = "tcell")
  tmp.lambdaz <- getLambdaZ(data = tmp)

  # result
  result <- list()

  if ("N_Samples" %in% metrics_to_calc) {
    result$N_Samples <- N_Samples = getSampleNumber(data = tmp)
  }

  if ("Rsq" %in% metrics_to_calc) {
    result$Rsq <- tmp.lambdaz$rsq
  }

  if ("Rsq_adjusted" %in% metrics_to_calc) {
    result$Rsq_adjusted <- tmp.lambdaz$rsqAdj
  }

  if ("Corr_XY" %in% metrics_to_calc) {
    result$Corr_XY <- tmp.lambdaz$corrXY  # from lambdaZ.R
  }

  if ("No_points_lambda_z" %in% metrics_to_calc) {
    result$No_points_lambda_z <- tmp.lambdaz$nSample
  }

  if ("Lambda_z" %in% metrics_to_calc) {
    result$Lambda_z <- tmp.lambdaz$lambdaZ
  }

  if ("Lambda_z_intercept" %in% metrics_to_calc) {
    result$Lambda_z_intercept <- tmp.lambdaz$Intercept  
  }

  if ("Lambda_z_lower" %in% metrics_to_calc) {
    result$Lambda_z_lower <- tmp.lambdaz$lambdazLower
  }

  if ("Lambda_z_upper" %in% metrics_to_calc) {
    result$Lambda_z_upper <- tmp.lambdaz$lambdazUpper
  }

  if ("Halflife_Lambda_z" %in% metrics_to_calc) {
    result$Halflife_Lambda_z <- tmp.lambdaz$halfLife
  }

  if ("Span" %in% metrics_to_calc) {
    result$Span <- tmp.lambdaz$lambdazSpan
  }

  if ("Tmax" %in% metrics_to_calc) {
    result$Tmax <- getTmax(data = tmp)
  }

  if ("Cmax" %in% metrics_to_calc) {
    result$Cmax <- getCmax(data = tmp)
  }

  if ("Cmax_D" %in% metrics_to_calc) {
    result$Cmax_D <- getCmaxD(data = tmp, Dose = dose)
  }

  if ("Tlast" %in% metrics_to_calc) {
    result$Tlast <- getTlast(data = tmp)
  }

  if ("Clast" %in% metrics_to_calc) {
    result$Clast <- getClast(data = tmp)
  }

  if ("Clast_pred" %in% metrics_to_calc) {
    result$Clast_pred <- tmp.lambdaz$clastPred
  }

  if ("AUClast" %in% metrics_to_calc) {
    result$AUClast <- getAucLast(data = tmp, method = "LinearUpLogDown")
  }

  if ("AUCall" %in% metrics_to_calc) {
    result$AUCall <- getAucAll(data = tmp, method = "LinearUpLogDown")
  }

  if ("AUCINF_obs" %in% metrics_to_calc) {
    result$AUCINF_obs <- getAucInf(data = tmp, type = "obs", method = "LinearUpLogDown")
  }

  if ("AUC_Extrap_obs" %in% metrics_to_calc) {
    result$AUC_Extrap_obs <- getAucPrecentExtrap(data = tmp, type = "obs", method = "LinearUpLogDown")
  }

  if ("AUCINF_pred" %in% metrics_to_calc) {
    result$AUCINF_pred <- getAucInf(data = tmp, type = "pred", method = "LinearUpLogDown")
  }

  if ("AUC_Extrap_pred" %in% metrics_to_calc) {
    result$AUC_Extrap_pred <- getAucPrecentExtrap(data = tmp, type = "pred", method = "LinearUpLogDown")
  }

  if ("AUMClast" %in% metrics_to_calc) {
    result$AUMClast <- getAumcLast(data = tmp, method = "LinearUpLogDown")
  }

  if ("AUMCINF_obs" %in% metrics_to_calc) {
    result$AUMCINF_obs <- getAumcInf(data = tmp, type = "obs", method = "LinearUpLogDown")
  }  

  if ("AUMC_Extrap_obs" %in% metrics_to_calc) {
    result$AUMC_Extrap_obs <- getAumcPrecentExtrap(data = tmp, type = "obs", method = "LinearUpLogDown")
  }

  if ("AUMCINF_pred" %in% metrics_to_calc) {
    result$AUMCINF_pred <- getAumcInf(data = tmp, type = "pred", method = "LinearUpLogDown")
  }

  if ("AUMC_Extrap_pred" %in% metrics_to_calc) {
    result$AUMC_Extrap_pred <- getAumcPrecentExtrap(data = tmp, type = "pred", method = "LinearUpLogDown")
  }

  if ("MRTlast" %in% metrics_to_calc) {
    result$MRTlast <- getMrtLast(data = tmp, method = "LinearUpLogDown")
  }

  if ("MRTINF_obs" %in% metrics_to_calc) {
    result$MRTINF_obs <- getMrtInf(data = tmp, type = "obs", method = "LinearUpLogDown")
  }

  if ("MRTINF_pred" %in% metrics_to_calc) {
    result$MRTINF_pred <- getMrtInf(data = tmp, type = "pred", method = "LinearUpLogDown")
  }

  if ("AUC7day" %in% metrics_to_calc) {
    result$AUC7day <- getAuc7D(data = tmp, method = "LinearUpLogDown")
  }

  if ("AUC28day" %in% metrics_to_calc) {
    result$AUC28day <- getAuc28D(data = tmp, method = "LinearUpLogDown")
  }

  if ("AUC3month" %in% metrics_to_calc) {
    result$AUC3month <- getAuc3M(data = tmp, method = "LinearUpLogDown")
  }

  if ("AUC6month" %in% metrics_to_calc) {
    result$AUC6month <- getAuc6M(data = tmp, method = "LinearUpLogDown")
  }

  if ("Vz_F_obs" %in% metrics_to_calc) {
    result$Vz_F_obs <- getVzF(data = tmp, type = "obs", Dose = tmp.dose, method = "LinearUpLogDown")
  }

  if ("Cl_F_obs" %in% metrics_to_calc) {
    result$Cl_F_obs <- getClF(data = tmp, type = "obs", Dose = tmp.dose, method = "LinearUpLogDown")
  }

  if ("Vz_F_pred" %in% metrics_to_calc) {
    result$Vz_F_pred <- getVzF(data = tmp, type = "pred", Dose = tmp.dose, method = "LinearUpLogDown")
  }

  if ("Cl_F_pred" %in% metrics_to_calc) {
    result$Cl_F_pred <- getClF(data = tmp, type = "pred", Dose = tmp.dose, method = "LinearUpLogDown")
  }  

  return(result)
}

calculate_nca_cytokine <- function(time, conc, dose, metrics = "all") {
  metrics_to_calc <- validate_metrics(metrics, .available_nca_metrics)

  # data preprocessing
  tmp <- dataQc(time = time, conc = conc, sampletype = "cytokine")
  tmp.lambdaz <- getLambdaZ(data = tmp)

  # result
  result <- list()

  if ("N_Samples" %in% metrics_to_calc) {
    result$N_Samples <- N_Samples = getSampleNumber(data = tmp)
  }

  if ("Rsq" %in% metrics_to_calc) {
    result$Rsq <- tmp.lambdaz$rsq
  }

  if ("Rsq_adjusted" %in% metrics_to_calc) {
    result$Rsq_adjusted <- tmp.lambdaz$rsqAdj
  }

  if ("Corr_XY" %in% metrics_to_calc) {
    result$Corr_XY <- tmp.lambdaz$corrXY  # from lambdaZ.R
  }

  if ("No_points_lambda_z" %in% metrics_to_calc) {
    result$No_points_lambda_z <- tmp.lambdaz$nSample
  }

  if ("Lambda_z" %in% metrics_to_calc) {
    result$Lambda_z <- tmp.lambdaz$lambdaZ
  }

  if ("Lambda_z_intercept" %in% metrics_to_calc) {
    result$Lambda_z_intercept <- tmp.lambdaz$Intercept  
  }

  if ("Lambda_z_lower" %in% metrics_to_calc) {
    result$Lambda_z_lower <- tmp.lambdaz$lambdazLower
  }

  if ("Lambda_z_upper" %in% metrics_to_calc) {
    result$Lambda_z_upper <- tmp.lambdaz$lambdazUpper
  }

  if ("Halflife_Lambda_z" %in% metrics_to_calc) {
    result$Halflife_Lambda_z <- tmp.lambdaz$halfLife
  }

  if ("Span" %in% metrics_to_calc) {
    result$Span <- tmp.lambdaz$lambdazSpan
  }

  if ("Tmax" %in% metrics_to_calc) {
    result$Tmax <- getTmax(data = tmp)
  }

  if ("Cmax" %in% metrics_to_calc) {
    result$Cmax <- getCmax(data = tmp)
  }

  if ("Cmax_D" %in% metrics_to_calc) {
    result$Cmax_D <- getCmaxD(data = tmp, Dose = dose)
  }

  if ("Tlast" %in% metrics_to_calc) {
    result$Tlast <- getTlast(data = tmp)
  }

  if ("Clast" %in% metrics_to_calc) {
    result$Clast <- getClast(data = tmp)
  }

  if ("Clast_pred" %in% metrics_to_calc) {
    result$Clast_pred <- tmp.lambdaz$clastPred
  }

  if ("AUClast" %in% metrics_to_calc) {
    result$AUClast <- getAucLast(data = tmp, method = "LinearUpLogDown")
  }

  if ("AUCall" %in% metrics_to_calc) {
    result$AUCall <- getAucAll(data = tmp, method = "LinearUpLogDown")
  }

  if ("AUCINF_obs" %in% metrics_to_calc) {
    result$AUCINF_obs <- getAucInf(data = tmp, type = "obs", method = "LinearUpLogDown")
  }

  if ("AUC_Extrap_obs" %in% metrics_to_calc) {
    result$AUC_Extrap_obs <- getAucPrecentExtrap(data = tmp, type = "obs", method = "LinearUpLogDown")
  }

  if ("AUCINF_pred" %in% metrics_to_calc) {
    result$AUCINF_pred <- getAucInf(data = tmp, type = "pred", method = "LinearUpLogDown")
  }

  if ("AUC_Extrap_pred" %in% metrics_to_calc) {
    result$AUC_Extrap_pred <- getAucPrecentExtrap(data = tmp, type = "pred", method = "LinearUpLogDown")
  }

  if ("AUMClast" %in% metrics_to_calc) {
    result$AUMClast <- getAumcLast(data = tmp, method = "LinearUpLogDown")
  }

  if ("AUMCINF_obs" %in% metrics_to_calc) {
    result$AUMCINF_obs <- getAumcInf(data = tmp, type = "obs", method = "LinearUpLogDown")
  }  

  if ("AUMC_Extrap_obs" %in% metrics_to_calc) {
    result$AUMC_Extrap_obs <- getAumcPrecentExtrap(data = tmp, type = "obs", method = "LinearUpLogDown")
  }

  if ("AUMCINF_pred" %in% metrics_to_calc) {
    result$AUMCINF_pred <- getAumcInf(data = tmp, type = "pred", method = "LinearUpLogDown")
  }

  if ("AUMC_Extrap_pred" %in% metrics_to_calc) {
    result$AUMC_Extrap_pred <- getAumcPrecentExtrap(data = tmp, type = "pred", method = "LinearUpLogDown")
  }

  if ("MRTlast" %in% metrics_to_calc) {
    result$MRTlast <- getMrtLast(data = tmp, method = "LinearUpLogDown")
  }

  if ("MRTINF_obs" %in% metrics_to_calc) {
    result$MRTINF_obs <- getMrtInf(data = tmp, type = "obs", method = "LinearUpLogDown")
  }

  if ("MRTINF_pred" %in% metrics_to_calc) {
    result$MRTINF_pred <- getMrtInf(data = tmp, type = "pred", method = "LinearUpLogDown")
  }

  if ("AUC7day" %in% metrics_to_calc) {
    result$AUC7day <- getAuc7D(data = tmp, method = "LinearUpLogDown")
  }

  if ("AUC28day" %in% metrics_to_calc) {
    result$AUC28day <- getAuc28D(data = tmp, method = "LinearUpLogDown")
  }

  if ("AUC3month" %in% metrics_to_calc) {
    result$AUC3month <- getAuc3M(data = tmp, method = "LinearUpLogDown")
  }

  if ("AUC6month" %in% metrics_to_calc) {
    result$AUC6month <- getAuc6M(data = tmp, method = "LinearUpLogDown")
  }

  if ("Vz_F_obs" %in% metrics_to_calc) {
    result$Vz_F_obs <- getVzF(data = tmp, type = "obs", Dose = tmp.dose, method = "LinearUpLogDown")
  }

  if ("Cl_F_obs" %in% metrics_to_calc) {
    result$Cl_F_obs <- getClF(data = tmp, type = "obs", Dose = tmp.dose, method = "LinearUpLogDown")
  }

  if ("Vz_F_pred" %in% metrics_to_calc) {
    result$Vz_F_pred <- getVzF(data = tmp, type = "pred", Dose = tmp.dose, method = "LinearUpLogDown")
  }

  if ("Cl_F_pred" %in% metrics_to_calc) {
    result$Cl_F_pred <- getClF(data = tmp, type = "pred", Dose = tmp.dose, method = "LinearUpLogDown")
  }  
  
  return(result)
}
