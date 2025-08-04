# R/generalParameters.R

# Data checking and pre-treatment
# Sorting: prior to analysis, data within each profile is sorted in ascending time order
# Inserting initial time points
# Data exclusions: 1. Missing values, 2. Data points preceding the dose time
dataQc <- function(time, conc, sampletype){
	# construct a data.table
	data <- data.table::data.table(time = time, conc = conc)

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
# N_Samples: This parameter reports the number of non-missing observations used in the analysis of the profile
getSampleNumber <- function(data){
    return(nrow(data))
}

# Tlag: Time of observation prior to the first observation with a measurable (non-zero) concentration

# Cmax: Maximum observed concentration
getCmax <- function(data){
    max_val <- max(data$conc, na.rm = T)
    peak <- NA_real_
    if(!is.na(max_val)){
        peak <- max_val
    }
    return(peak)
}

# Tmax: Time of maximum observed concentration
getTmax <- function(data){
    max_val <- max(data$conc, na.rm = T)
    peak_time <- NA_real_
    if(!is.na(max_val)){
        peak_time <- min(data$time[data$conc==max_val], na.rm = T)
    }
    return(peak_time)
}

# Cmax_D: Cmax/Dose
getCmaxD <- function(data, Dose){
    return(getCmax(data)/Dose)
}

# Tlast: Time of last measurable (positive) observed concentration	
getTlast <- function(data){
    return(max(data$time))
}

# Clast: Observed concentration corresponding to Tlast
getClast <- function(data){
    return(tail(data$conc, n=1))
}
