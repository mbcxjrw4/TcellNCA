# Predefined list of available metrics
.available_nca_metrics <- c(
    "NoSamples", "Rsq", "RsqAdjusted", "CorrXY", "NoPointsLambdaz", "Lambdaz", "LambdazIntercept", "LambdazLower", "LambdazUpper", "HalfLife", "Span", "Tmax", "Cmax", "CmaxD", "Tlast", "Clast", "ClastPred", "AucLast", "AucAll", "AucInfObs", "AucExtrapObs", "AucInfPred", "AucExtrapPred", "AumcLast", "AumcInfObs", "AumcExtrapObs", "AumcInfPred", "AumcExtrapPred", "MrtLast", "MrtInfObs", "MrtInfPred", "AucD7", "AucD28", "AucM3", "AucM6", "VolumeOfDistributionObs", "ClearanceObs", "VolumeOfDistributionPred", "ClearancePred"
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
#' @param conc Numeric vector of concentrations (e.g., VCN, absolute T-cell numbers)
#' @param dose Dosage
#' @param metrics Character vector of metrics to calculate, or "all"
#' @return Named list of NCA results
#' @import data.table
#' @export
calculate_nca_tcell <- function(time, conc, dose, metrics = "all") {
    metrics_to_calc <- validate_metrics(metrics, .available_nca_metrics)

    # data preprocessing
    tmp <- dataQc(time = time, conc = conc, sampletype = "tcell")
    tmp.lambdaz <- getLambdaZ(data = tmp)

    # result
    result <- list()

    if ("NoSamples" %in% metrics_to_calc) {
        result$NoSamples <- getSampleNumber(data = tmp)
    }

    if ("Rsq" %in% metrics_to_calc) {
        result$Rsq <- tmp.lambdaz$rsq
    }

    if ("RsqAdjusted" %in% metrics_to_calc) {
        result$RsqAdjusted <- tmp.lambdaz$rsqAdj
    }

    if ("CorrXY" %in% metrics_to_calc) {
        result$CorrXY <- tmp.lambdaz$corrXY  # from lambdaZ.R
    }

    if ("NoPointsLambdaz" %in% metrics_to_calc) {
        result$NoPointsLambdaz <- tmp.lambdaz$nSample
    }

    if ("Lambdaz" %in% metrics_to_calc) {
        result$Lambdaz <- tmp.lambdaz$lambdaZ
    }

    if ("LambdazIntercept" %in% metrics_to_calc) {
        result$LambdazIntercept <- tmp.lambdaz$Intercept
    }

    if ("LambdazLower" %in% metrics_to_calc) {
        result$LambdazLower <- tmp.lambdaz$lambdazLower
    }

    if ("LambdazUpper" %in% metrics_to_calc) {
        result$LambdazUpper <- tmp.lambdaz$lambdazUpper
    }

    if ("HalfLife" %in% metrics_to_calc) {
        result$HalfLife <- tmp.lambdaz$halfLife
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

    if ("CmaxD" %in% metrics_to_calc) {
        result$CmaxD <- getCmaxD(data = tmp, Dose = dose)
    }

    if ("Tlast" %in% metrics_to_calc) {
        result$Tlast <- getTlast(data = tmp)
    }

    if ("Clast" %in% metrics_to_calc) {
        result$Clast <- getClast(data = tmp)
    }

    if ("ClastPred" %in% metrics_to_calc) {
        result$ClastPred <- tmp.lambdaz$clastPred
    }

    if ("AucLast" %in% metrics_to_calc) {
        result$AucLast <- getAucLast(data = tmp, method = "LinearUpLogDown")
    }

    if ("AucAll" %in% metrics_to_calc) {
        result$AucAll <- getAucAll(data = tmp, method = "LinearUpLogDown")
    }

    if ("AucInfObs" %in% metrics_to_calc) {
        result$AucInfObs <- getAucInf(data = tmp, type = "obs", method = "LinearUpLogDown")
    }

    if ("AucExtrapObs" %in% metrics_to_calc) {
        result$AucExtrapObs <- getAucPrecentExtrap(data = tmp, type = "obs", method = "LinearUpLogDown")
    }

    if ("AucInfPred" %in% metrics_to_calc) {
        result$AucInfPred <- getAucInf(data = tmp, type = "pred", method = "LinearUpLogDown")
    }

    if ("AucExtrapPred" %in% metrics_to_calc) {
        result$AucExtrapPred <- getAucPrecentExtrap(data = tmp, type = "pred", method = "LinearUpLogDown")
    }

    if ("AumcLast" %in% metrics_to_calc) {
        result$AumcLast <- getAumcLast(data = tmp, method = "LinearUpLogDown")
    }

    if ("AumcInfObs" %in% metrics_to_calc) {
        result$AumcInfObs <- getAumcInf(data = tmp, type = "obs", method = "LinearUpLogDown")
    }

    if ("AumcExtrapObs" %in% metrics_to_calc) {
        result$AumcExtrapObs <- getAumcPrecentExtrap(data = tmp, type = "obs", method = "LinearUpLogDown")
    }

    if ("AumcInfPred" %in% metrics_to_calc) {
        result$AumcInfPred <- getAumcInf(data = tmp, type = "pred", method = "LinearUpLogDown")
    }

    if ("AumcExtrapPred" %in% metrics_to_calc) {
        result$AumcExtrapPred <- getAumcPrecentExtrap(data = tmp, type = "pred", method = "LinearUpLogDown")
    }

    if ("MrtLast" %in% metrics_to_calc) {
        result$MrtLast <- getMrtLast(data = tmp, method = "LinearUpLogDown")
    }

    if ("MrtInfObs" %in% metrics_to_calc) {
        result$MrtInfObs <- getMrtInf(data = tmp, type = "obs", method = "LinearUpLogDown")
    }

    if ("MrtInfPred" %in% metrics_to_calc) {
        result$MrtInfPred <- getMrtInf(data = tmp, type = "pred", method = "LinearUpLogDown")
    }

    if ("AucD7" %in% metrics_to_calc) {
        result$AucD7 <- getAuc7D(data = tmp, method = "LinearUpLogDown")
    }

    if ("AucD28" %in% metrics_to_calc) {
        result$AucD28 <- getAuc28D(data = tmp, method = "LinearUpLogDown")
    }

    if ("AucM3" %in% metrics_to_calc) {
        result$AucM3 <- getAuc3M(data = tmp, method = "LinearUpLogDown")
    }

    if ("AucM6" %in% metrics_to_calc) {
        result$AucM6 <- getAuc6M(data = tmp, method = "LinearUpLogDown")
    }

    if ("VolumeOfDistributionObs" %in% metrics_to_calc) {
        result$VolumeOfDistributionObs <- getVzF(data = tmp, type = "obs", Dose = dose, method = "LinearUpLogDown")
    }

    if ("ClearanceObs" %in% metrics_to_calc) {
        result$ClearanceObs <- getClF(data = tmp, type = "obs", Dose = dose, method = "LinearUpLogDown")
    }

    if ("VolumeOfDistributionPred" %in% metrics_to_calc) {
        result$VolumeOfDistributionPred <- getVzF(data = tmp, type = "pred", Dose = dose, method = "LinearUpLogDown")
    }

    if ("ClearancePred" %in% metrics_to_calc) {
        result$ClearancePred <- getClF(data = tmp, type = "pred", Dose = dose, method = "LinearUpLogDown")
    }

    return(result)
}


#' Calculate NCA metrics for serum cytokine data
#'
#' @param time Numeric vector of time points
#' @param conc Numeric vector of concentrations (e.g., pg/mL)
#' @param dose Dosage
#' @param metrics Character vector of metrics to calculate, or "all"
#' @return Named list of NCA results
#' @import data.table
#' @export
calculate_nca_cytokine <- function(time, conc, dose, metrics = "all") {
    metrics_to_calc <- validate_metrics(metrics, .available_nca_metrics)

    # data preprocessing
    tmp <- dataQc(time = time, conc = conc, sampletype = "cytokine")
    tmp.lambdaz <- getLambdaZ(data = tmp)

    # result
    result <- list()

    if ("NoSamples" %in% metrics_to_calc) {
        result$NoSamples <- getSampleNumber(data = tmp)
    }

    if ("Rsq" %in% metrics_to_calc) {
        result$Rsq <- tmp.lambdaz$rsq
    }

    if ("RsqAdjusted" %in% metrics_to_calc) {
        result$RsqAdjusted <- tmp.lambdaz$rsqAdj
    }

    if ("CorrXY" %in% metrics_to_calc) {
        result$CorrXY <- tmp.lambdaz$corrXY  # from lambdaZ.R
    }

    if ("NoPointsLambdaz" %in% metrics_to_calc) {
        result$NoPointsLambdaz <- tmp.lambdaz$nSample
    }

    if ("Lambdaz" %in% metrics_to_calc) {
        result$Lambdaz <- tmp.lambdaz$lambdaZ
    }

    if ("LambdazIntercept" %in% metrics_to_calc) {
        result$LambdazIntercept <- tmp.lambdaz$Intercept
    }

    if ("LambdazLower" %in% metrics_to_calc) {
        result$LambdazLower <- tmp.lambdaz$lambdazLower
    }

    if ("LambdazUpper" %in% metrics_to_calc) {
        result$LambdazUpper <- tmp.lambdaz$lambdazUpper
    }

    if ("HalfLife" %in% metrics_to_calc) {
        result$HalfLife <- tmp.lambdaz$halfLife
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

    if ("CmaxD" %in% metrics_to_calc) {
        result$CmaxD <- getCmaxD(data = tmp, Dose = dose)
    }

    if ("Tlast" %in% metrics_to_calc) {
        result$Tlast <- getTlast(data = tmp)
    }

    if ("Clast" %in% metrics_to_calc) {
        result$Clast <- getClast(data = tmp)
    }

    if ("ClastPred" %in% metrics_to_calc) {
        result$ClastPred <- tmp.lambdaz$clastPred
    }

    if ("AucLast" %in% metrics_to_calc) {
        result$AucLast <- getAucLast(data = tmp, method = "Linear")
    }

    if ("AucAll" %in% metrics_to_calc) {
        result$AucAll <- getAucAll(data = tmp, method = "Linear")
    }

    if ("AucInfObs" %in% metrics_to_calc) {
        result$AucInfObs <- getAucInf(data = tmp, type = "obs", method = "Linear")
    }

    if ("AucExtrapObs" %in% metrics_to_calc) {
        result$AucExtrapObs <- getAucPrecentExtrap(data = tmp, type = "obs", method = "Linear")
    }

    if ("AucInfPred" %in% metrics_to_calc) {
        result$AucInfPred <- getAucInf(data = tmp, type = "pred", method = "Linear")
    }

    if ("AucExtrapPred" %in% metrics_to_calc) {
        result$AucExtrapPred <- getAucPrecentExtrap(data = tmp, type = "pred", method = "Linear")
    }

    if ("AumcLast" %in% metrics_to_calc) {
        result$AumcLast <- getAumcLast(data = tmp, method = "Linear")
    }

    if ("AumcInfObs" %in% metrics_to_calc) {
        result$AumcInfObs <- getAumcInf(data = tmp, type = "obs", method = "Linear")
    }

    if ("AumcExtrapObs" %in% metrics_to_calc) {
        result$AumcExtrapObs <- getAumcPrecentExtrap(data = tmp, type = "obs", method = "Linear")
    }

    if ("AumcInfPred" %in% metrics_to_calc) {
        result$AumcInfPred <- getAumcInf(data = tmp, type = "pred", method = "Linear")
    }

    if ("AumcExtrapPred" %in% metrics_to_calc) {
        result$AumcExtrapPred <- getAumcPrecentExtrap(data = tmp, type = "pred", method = "Linear")
    }

    if ("MrtLast" %in% metrics_to_calc) {
        result$MrtLast <- getMrtLast(data = tmp, method = "Linear")
    }

    if ("MrtInfObs" %in% metrics_to_calc) {
        result$MrtInfObs <- getMrtInf(data = tmp, type = "obs", method = "Linear")
    }

    if ("MrtInfPred" %in% metrics_to_calc) {
        result$MrtInfPred <- getMrtInf(data = tmp, type = "pred", method = "Linear")
    }

    if ("AucD7" %in% metrics_to_calc) {
        result$AucD7 <- getAuc7D(data = tmp, method = "Linear")
    }

    if ("AucD28" %in% metrics_to_calc) {
        result$AucD28 <- getAuc28D(data = tmp, method = "Linear")
    }

    if ("AucM3" %in% metrics_to_calc) {
        result$AucM3 <- getAuc3M(data = tmp, method = "Linear")
    }

    if ("AucM6" %in% metrics_to_calc) {
        result$AucM6 <- getAuc6M(data = tmp, method = "Linear")
    }

    if ("VolumeOfDistributionObs" %in% metrics_to_calc) {
        result$VolumeOfDistributionObs <- getVzF(data = tmp, type = "obs", Dose = dose, method = "Linear")
    }

    if ("ClearanceObs" %in% metrics_to_calc) {
        result$ClearanceObs <- getClF(data = tmp, type = "obs", Dose = dose, method = "Linear")
    }

    if ("VolumeOfDistributionPred" %in% metrics_to_calc) {
        result$VolumeOfDistributionPred <- getVzF(data = tmp, type = "pred", Dose = dose, method = "Linear")
    }

    if ("ClearancePred" %in% metrics_to_calc) {
        result$ClearancePred <- getClF(data = tmp, type = "pred", Dose = dose, method = "Linear")
    }

    return(result)
}
