# auc.R

#' Log-linear concentration interpolation (internal)
#'
#' Estimates a concentration at a given timepoint using log-linear interpolation
#' between two known concentrations and their corresponding timepoints.
#'
#' @param c1 Numeric. Concentration at time `t1`.
#' @param c2 Numeric. Concentration at time `t2`.
#' @param t1 Numeric. First timepoint.
#' @param t2 Numeric. Second timepoint.
#' @param timepoint Numeric. The timepoint at which to estimate the concentration.
#'
#' @return Estimated concentration at the specified `timepoint`.
#' @keywords internal
logestimate <- function(c1, c2, t1, t2, timepoint){
    return(c1 * exp(log(c2/c1) * (timepoint - t1) / (t2 - t1)))
}

#' Calculate AUC using Linear Trapezoidal Method
#'
#' @param c1 Numeric. Concentration at time `t1`.
#' @param c2 Numeric. Concentration at time `t2`.
#' @param t1 Numeric. First timepoint.
#' @param t2 Numeric. Second timepoint.
#'
#' @return A single numeric value representing calculated AUC
#' @keywords internal
linearauc <- function(c1, c2, t1, t2){
    return(0.5 * (c1 + c2) * (t2 - t1))
}

#' Calculate AUMC (Area under the first moment curve) using Linear Trapezoidal Method
#'
#' @param c1 Numeric. Concentration at time `t1`.
#' @param c2 Numeric. Concentration at time `t2`.
#' @param t1 Numeric. First timepoint.
#' @param t2 Numeric. Second timepoint.
#'
#' @return A single numeric value representing calculated AUMC
#' @keywords internal
linearaumc <- function(c1, c2, t1, t2){
    return(0.5 * (c1 * t1 + c2 * t2) * (t2 - t1))
}

#' Calculate AUC using Logarithmic Trapezoidal Method
#'
#' @param c1 Numeric. Concentration at time `t1`.
#' @param c2 Numeric. Concentration at time `t2`.
#' @param t1 Numeric. First timepoint.
#' @param t2 Numeric. Second timepoint.
#'
#' @return A single numeric value representing calculated AUC
#' @keywords internal
logauc <- function(c1, c2, t1, t2){
    return((c1 - c2) * (t2 - t1) / (log(c1) - log(c2)))
}

#' Calculate AUMC using Logarithmic Trapezoidal Method
#'
#' @param c1 Numeric. Concentration at time `t1`.
#' @param c2 Numeric. Concentration at time `t2`.
#' @param t1 Numeric. First timepoint.
#' @param t2 Numeric. Second timepoint.
#'
#' @return A single numeric value representing calculated AUMC
#' @keywords internal
logaumc <- function(c1, c2, t1, t2){
    return( (t2 * c2 - t1 * c1) * (t2 - t1) / log(c2/c1) - (t2 - t1)^2 * (c2 - c1) / (log(c2 / c1)^2) )
}

#' Calculate AUC using Linear-Log Trapezoidal Method
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param auc.type `auc` or `aumc`
#' @param startpoint Numeric. First timepoint.
#' @param endpoint Numeric. Second timepoint.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#'
#' @return A single numeric value representing calculated AUC/AUMC
#' @import data.table
#' @keywords internal
auc_cal <- function(data, startpoint, endpoint, auc.type, method){
    x <- data$time
    y <- data$conc
    len <- length(x)
    if(len < 2){
        return(NaN)
    }else{
        if(endpoint < min(x)){
            return(NaN)
        }else{
            # order x and y by ascending x
            y <- y[order(x)]
            x <- x[order(x)]

            # exclude data prior to start point
            y <- y[x >= startpoint]
            x <- x[x >= startpoint]

            if(endpoint > max(x)){
                # if the end time occurs after the last numeric observation and Lambda Z is estimable, Lambda Z is used to estimate the corresponding Y
                tmp.lambdaz <- getLambdaZ(data)
                y_end <- exp(tmp.lambdaz$Intercept - tmp.lambdaz$lambdaZ * endpoint)
                rm(tmp.lambdaz)
            }
            else if(!(endpoint %in% x)){
                # if end time point is not in x, then the value will be estimated by linear interpolation (while increasing) or a logarithmic decline
                # get two adjacent points
                x_adj <- (x - endpoint)
                left_point <- max(x_adj[x_adj < 0])
                left_index <- which(x_adj == left_point)
                x_left <- x[left_index]
                y_left <- y[left_index]
                x_right <- x[left_index + 1]
                y_right <- y[left_index + 1]

                if(method == "LinearUpLogDown"){
                    if(y_left <= y_right){
                        y_end <- stats::approx(x=c(x_left, x_right), y=c(y_left, y_right), xout=endpoint)$y
                    }else{
                        y_end <- logestimate(c1 = y_left, t1 = x_left, c2 = y_right, t2 = x_right, timepoint = endpoint)
                    }
                }

                if(method == "Linear"){
                    y_end <- stats::approx(x=c(x_left, x_right), y=c(y_left, y_right), xout=endpoint)$y
                }

                if(method == "Log"){
                    if(y_left == y_right){
                        y_end <- y_left
                    }else{
                        y_end <- logestimate(c1 = y_left, t1 = x_left, c2 = y_right, t2 = x_right, timepoint = endpoint)
                    }
                }

                rm(x_adj, left_point, left_index, x_left, x_right, y_left, y_right)
            }else{
                y_end <- y[x == endpoint]
            }

            # update x and y
            y <- c(y[x<endpoint], y_end)
            x <- c(x[x<endpoint], endpoint)
            len <- length(x)

            tmp <- data.frame(C1 = y[1:(len-1)], C2 = y[2:len], T1 = x[1:(len-1)], T2 = x[2:len])

            tmp <- data.table::as.data.table(tmp)
            tmp[, Diff := C2 - C1]

            if(auc.type == "auc"){
                if(method == "LinearUpLogDown"){
                    tmp[Diff >= 0, subAUC := linearauc(c1 = C1, c2 = C2, t1 = T1, t2 = T2)]
                    tmp[Diff <  0, subAUC := logauc(c1 = C1, c2 = C2, t1 = T1, t2 = T2)]
                }

                if(method == "Linear"){
                    tmp[, subAUC := linearauc(c1 = C1, c2 = C2, t1 = T1, t2 = T2)]
                }

                if(method == "Log"){
                    tmp[Diff == 0, subAUC := linearauc(c1 = C1, c2 = C2, t1 = T1, t2 = T2)]
                    tmp[Diff !=  0, subAUC := logauc(c1 = C1, c2 = C2, t1 = T1, t2 = T2)]
                }
            }

            if(auc.type == "aumc"){
                if(method == "LinearUpLogDown"){
                    tmp[Diff >= 0, subAUC := linearaumc(c1 = C1, c2 = C2, t1 = T1, t2 = T2)]
                    tmp[Diff <  0, subAUC := logaumc(c1 = C1, c2 = C2, t1 = T1, t2 = T2)]
                }

                if(method == "Linear"){
                    tmp[, subAUC := linearaumc(c1 = C1, c2 = C2, t1 = T1, t2 = T2)]
                }

                if(method == "Log"){
                    tmp[Diff == 0, subAUC := linearaumc(c1 = C1, c2 = C2, t1 = T1, t2 = T2)]
                    tmp[Diff !=  0, subAUC := logaumc(c1 = C1, c2 = C2, t1 = T1, t2 = T2)]
                }
            }

            return(sum(tmp$subAUC))
        }
    }
}

#' Calculate AUC from day0 to day7
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#'
#' @return A single numeric value representing calculated AUC
#' @keywords internal
getAuc7D <- function(data, method){
    return(auc_cal(data = data, startpoint = 0, endpoint = 7, auc.type = "auc", method = method))
}

#' Calculate AUC from day0 to day28
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#'
#' @return A single numeric value representing calculated AUC
#' @keywords internal
getAuc28D <- function(data, method){
    return(auc_cal(data = data, startpoint = 0, endpoint = 28, auc.type = "auc", method = method))
}

#' Calculate AUC from day0 to month3
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#'
#' @return A single numeric value representing calculated AUC
#' @keywords internal
getAuc3M <- function(data, method){
    return(auc_cal(data = data, startpoint = 0, endpoint = 84, auc.type = "auc", method = method))
}

#' Calculate AUC from day0 to month6
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#'
#' @return A single numeric value representing calculated AUC
#' @keywords internal
getAuc6M <- function(data, method){
    return(auc_cal(data = data, startpoint = 0, endpoint = 168, auc.type = "auc", method = method))
}

#' Calculate AUC from day0 to month12
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#'
#' @return A single numeric value representing calculated AUC
#' @keywords internal
getAuc12M <- function(data, method){
    return(auc_cal(data = data, startpoint = 0, endpoint = 366, auc.type = "auc", method = method))
}

#' Calculate AUC from the time of dosing to the time of the last measurable (positive) concentration (Tlast)
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#'
#' @return A single numeric value representing calculated AUC
#' @keywords internal
getAucLast <- function(data, method){
    return(auc_cal(data = data, startpoint = 0, endpoint = getTlast(data), auc.type = "auc", method = method))
}


# AUClast_D: AUClast/Dose
#' Calculate AUClast/Dose
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#' @param Dose Dosage
#'
#' @return A single numeric value representing calculated AUC
#' @keywords internal
getAucLastD <- function(data, Dose, method){
    return(getAucLast(data, method)/Dose)
}

#' Calculate AUC from the time of dosing to the time of the last observation. If the last concentration is positive, AUClast=AUCall.
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#'
#' @return A single numeric value representing calculated AUC
#' @keywords internal
getAucAll <- function(data, method){
    return(auc_cal(data = data, startpoint = 0, endpoint = getTlast(data), auc.type = "auc", method = method))
}

#' Calculate area under the moment curve from the time of dosing to the last measurable (positive) concentration
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#'
#' @return A single numeric value representing calculated AUMC
#' @keywords internal
getAumcLast <- function(data, method){
    return(auc_cal(data = data, startpoint = 0, endpoint = getTlast(data), auc.type = "aumc", method = method))
}

#' Calculate AUC from time of dosing extrapolated to infinity, based on the last observed concentration (_obs) or last predicted concentration (_pred). AUClast + (Clast/Lambda_z)
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#' @param type `pred` or `obs`
#'
#' @return A single numeric value representing calculated AUC
#' @keywords internal
getAucInf <- function(data, type, method){
    tmp <- getLambdaZ(data)
    res <- auc_cal(data = data, startpoint = 0, endpoint = getTlast(data), auc.type = "auc", method = method)

    if(type == "pred"){
        res <- res + (tmp$clastPred / tmp$lambdaZ)
    }

    if(type == "obs"){
        res <- res + (getClast(data) / tmp$lambdaZ)
    }
    return(res)
}

#' Calculate AUCINF/Dose
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#' @param type `pred` or `obs`
#' @param Dose dosage
#'
#' @return A single numeric value representing calculated AUCINF/Dose
#' @keywords internal
getAucInfD <- function(data, type, Dose, method){
    return( getAucInf(data, type, method) / Dose)
}

#' Calculate percentage of AUCINF(_obs, _pred) due to extrapolation from Tlast to infinity, 100[(AUCINF– AUClast)/AUCINF]
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#' @param type `pred` or `obs`
#'
#' @return A single numeric value representing calculated percentage of AUCINF(_obs, _pred) due to extrapolation from Tlast to infinity
#' @keywords internal
getAucPrecentExtrap <- function(data, type, method){
    aucInf <- getAucInf(data, type = type, method = method)
    aucLast <- getAucLast(data, method = method)

    return(100 * (aucInf - aucLast) / aucInf )
}

#' Calculate Area under the first moment curve (AUMC) extrapolated to infinity, based on the last observed concentration (obs) or the last predicted concentration (pred)
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#' @param type `pred` or `obs`
#'
#' @return A single numeric value representing calculated AUMC
#' @keywords internal
getAumcInf <- function(data, type, method){
    tmp <- getLambdaZ(data)
    res <- auc_cal(data = data, startpoint = 0, endpoint = getTlast(data), auc.type = "aumc", method = method)

    if(type == "pred"){
        res <- res + (getTlast(data) * tmp$clastPred / tmp$lambdaZ) + (tmp$clastPred / (tmp$lambdaZ^2))
    }

    if(type == "obs"){
        res <- res + (getTlast(data) * getClast(data) / tmp$lambdaZ) + getClast(data) / ((tmp$lambdaZ^2))
    }
    return(res)
}

#' Calculate Percent of AUMCINF(_obs, _pred) that is extrapolated
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#' @param type `pred` or `obs`
#'
#' @return A single numeric value representing Percent of AUMCINF(_obs, _pred) that is extrapolated
#' @keywords internal
getAumcPrecentExtrap <- function(data, type, method){
    aumcInf <- getAumcInf(data, type = type, method = method)
    aumcLast <- getAumcLast(data, method = method)

    return(100 * (aumcInf - aumcLast) / aumcInf )
}

#' Calculate mean residence time from the time of dosing to the time of the last measurable concentration. AUC_0_last / AUC_infinity
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#'
#' @return A single numeric value representing calculated MRT
#' @keywords internal
getMrtLast <- function(data, method){
    aumcLast <- getAumcLast(data, method = method)
    aucLast <- getAucLast(data, method = method)
    return(aumcLast/aucLast)
}

#' Calculate mean residence time (MRT) extrapolated to infinity based on AUCINF(_obs, _pred)
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#' @param type `pred` or `obs`
#'
#' @return A single numeric value representing calcuated MRT
#' @keywords internal
getMrtInf <- function(data, type, method){
    aumcInf <- getAumcInf(data, type, method = method)
    aucInf <- getAucInf(data, type, method = method)
    return(aumcInf/aucInf)
}

#' Calculate volume of distribution based on the terminal phase, Dose / (AucInf * LambdaZ)
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#' @param type `pred` or `obs`
#' @param Dose dosage
#'
#' @return A single numeric value representing calcuated volume of distribution
#' @keywords internal
getVzF <- function(data, type, Dose, method){
    tmp <- getLambdaZ(data)
    lambdaz <- tmp$lambdaZ
    return( Dose / (getAucInf(data, type, method = method) * lambdaz))
}

#' Calculate total body clearance for extravascular administration, Dose/AUCINF
#'
#' @param data A `data.table` with columns `time` and `conc`.
#'             Time must be numeric and strictly increasing.
#' @param method `LinearUpLogDown`, `Linear`, `Log`
#' @param type `pred` or `obs`
#' @param Dose dosage
#'
#' @return A single numeric value representing calcuated clearance
#' @keywords internal
getClF <- function(data, type, Dose, method){
    return(Dose / getAucInf(data, type, method = method))
}
