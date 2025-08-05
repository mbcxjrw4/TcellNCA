# R/generalParameters.R

#' Data checking and pre-treatment
#' Sorting: prior to analysis, data within each profile is sorted in ascending time order
#' Inserting initial time points
#' Data exclusions: 1. Missing values, 2. Data points preceding the dose time
#'
#' @param time Numeric vector of time points
#' @param conc Numeric vector of concentrations
#' @param sampletype sample type: "tcell" or "cytokine"
#' @return QCed data.table
#' @import data.table
#' @keywords internal
dataQc <- function(time, conc, sampletype){
    # construct a data.table
    data <- data.table(time = time, conc = conc)

    # data exclusion
    data <- data[!is.na(conc)]

    # inserting initial time points
    # for T-cell persistence data, the concentration is 1 VCN/ATC.
    if(sampletype == "tcell"){
        data <- data[time >= 0]

        if(!(0 %in% data$time)){
            initial_tp <- list(0, 1)
            data <- rbind(data, initial_tp)
        }
    }

    # for serum cytokine data, the concentration is pre-infusion collection concentration
    if(sampletype == "cytokine"){
        data[time == (-0.5), time := 0]

        if(!(0 %in% data$time)){
            initial_tp <- list(0, data$conc[data$time<0][1])
        }
    }

    data <- data[time >= 0]

    # sorting
    data <- data[order(time)]

    return(data)
}

# parameters that do not require Lambda Z estimation

#' Extract the parameter reports the number of non-missing observations used in the analysis of the profile
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @return A single numeric value representing the number of observations
#' @keywords internal
getSampleNumber <- function(data){
    return(nrow(data))
}

#' Extract the parameter reports the maximum observed concentration
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @return A single numeric value representing the maximum observed concentration
#' @keywords internal
getCmax <- function(data){
    max_val <- max(data$conc, na.rm = T)
    peak <- NA_real_
    if(!is.na(max_val)){
        peak <- max_val
    }
    return(peak)
}

#' Extract the parameter reports the time of maximum observed concentration
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @return A single numeric value representing the time of maximum observed concentration
#' @keywords internal
getTmax <- function(data){
    max_val <- max(data$conc, na.rm = T)
    peak_time <- NA_real_
    if(!is.na(max_val)){
        peak_time <- min(data$time[data$conc==max_val], na.rm = T)
    }
    return(peak_time)
}

#' Calculate the parameter reports the ratio of Cmax and dosage
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param Dose Dosage
#' @return A single numeric value representing Cmax/Dose
#' @keywords internal
getCmaxD <- function(data, Dose){
    return(getCmax(data)/Dose)
}

#' Extract the parameter reports the time of last measurable (positive) observed concentration
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @return A single numeric value representing the time of last measurable (positive) observed concentration
#' @keywords internal
getTlast <- function(data){
    return(max(data$time))
}

#' Extract the parameter reports the observed concentration corresponding to Tlast
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @return A single numeric value representing the observed concentration corresponding to Tlast
#' @keywords internal
getClast <- function(data){
    return(tail(data$conc, n=1))
}
