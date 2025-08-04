# Rsq: Goodness of fit statistic for the terminal elimination phase
calculateRsq <- function(X, Y){
    # fit the linear model
    model <- lm(Y~X)
    
    # get the summary of the model
    r.sq <- summary(model)$r.squared
    
    return(r.sq)
}

# Rsq_adjusted: Goodness of fit statistic for the terminal elimination phase, adjusted for the number of points used in the estimation of Lambda Z
calculateRsqAdj <- function(X, Y){
    # get the R-squared value
    r.sq <- calculateRsq(X, Y)
    
    # calculate adjusted r-squared value
    N <- length(Y)
    r.sq.adj <- 1 - ((1 - r.sq) * (N - 1) / (N - 2))
    
    return(r.sq.adj)
}

# parameters that are estimated when Lambda Z is estimated
getLambdaZ <- function(data){
    # points satisfying the condition:Time ≥ Tmax
    t.max <- getTmax(data)
    df <- data[time > t.max]
    row.num <- nrow(df)
    
    if((row.num < 3) | (sd(df$conc)) == 0){
        res <- list(nSample = NaN, rsq = NaN, rsqAdj = NaN, corrXY = NaN, lambdaZ = NaN, Intercept = NaN, lambdazLower = NaN, lambdazUpper = NaN, halfLife = NaN, lambdazSpan = NaN, clastPred = NaN)
    }else{
        # transform concentration for lambdaZ calculation
        df[, lnValue := log(conc)]
    
        # data used to calculate Lambda Z
        for(i in seq(3, row.num, by = 1)){
            tmp <- df[(row.num - i + 1) : row.num]
            
            if(sd(tmp$conc) != 0){
                # check lambda Z
                # fit the linear model
                tmp.model <- lm(tmp$lnValue ~ tmp$time)
            
                # Extract the rate constant
                tmp.lambdaz <- (-summary(tmp.model)$coefficients[2])
            
                if(tmp.lambdaz > 0){
                    tmp.r.sq.adj <- calculateRsqAdj(X = tmp$time, Y = tmp$lnValue)
        
                    if(!exists("max.r.sq.adj")){
                        max.r.sq.adj <- tmp.r.sq.adj
                        max.i <- i
                    }else{
                        if(tmp.r.sq.adj > max.r.sq.adj){
                            max.r.sq.adj <- tmp.r.sq.adj
                            max.i <- i
                        }
                    }
                }
            }
        }
        
        if(exists("max.i")){
            tmp <- df[(row.num - max.i + 1) : row.num]
    
            # parameters related to Lambda Z
            n.samples.lambdaz <- getSampleNumber(tmp) # No_points_lambda_z: Number of points used in computing Lambda Z 
    
            rsq.lambdaz <- calculateRsq(X = tmp$time, Y = tmp$lnValue) # Rsq
    
            rsq.adj.lambdaz <- calculateRsqAdj(X = tmp$time, Y = tmp$lnValue) # Rsq_adjusted
    
            # Corr_XY: Correlation between time (X) and log concentration (Y) for the points used in the estimation of Lambda Z.
            corr_xy <- cor(tmp$time, tmp$lnValue, method = "pearson") 
        
            # Lambda_z: First-order rate constant associated with the terminal (log-linear) portion of the curve. Estimated by linear regression of time vs log concentration.
            # fit the linear model
            model <- lm(tmp$lnValue ~ tmp$time)
    
            # Extract the rate constant
            lambdaz <- (-summary(model)$coefficients[2])
        
            # Lambda_z_intercept: Intercept on log scale estimated via linear regression of time vs. log concentration.
            intercept <- summary(model)$coefficients[1]
    
            # Lambda_z_lower: Lower limit on time for values to be included in the calculation of Lambda Z
            lambdaz.lower <- min(tmp$time)
    
            # Lambda_z_upper:	Upper limit on time for values to be included in the calculation of Lambda Z
            lambdaz.upper <- max(tmp$time)
    
            # HL_Lambda_z: Terminal half-life
            halflife <- log(2)/lambdaz
    
            # Span: (Lambda_z_upper – Lambda_z_lower)/HL_Lambda_z 
            lambdaz.span <- (lambdaz.upper - lambdaz.lower) / halflife
    
            # Clast_pred: Predicted concentration at Tlast
            tlast <- getTlast(tmp)
            clast.pred <- exp(intercept - lambdaz*tlast)
    
            res <- list(nSample = n.samples.lambdaz, rsq = rsq.lambdaz, rsqAdj = rsq.adj.lambdaz, corrXY = corr_xy, lambdaZ = lambdaz, Intercept = intercept, lambdazLower = lambdaz.lower, lambdazUpper = lambdaz.upper, halfLife = halflife, lambdazSpan = lambdaz.span, clastPred = clast.pred)
        }else{
            res <- list(nSample = NaN, rsq = NaN, rsqAdj = NaN, corrXY = NaN, lambdaZ = NaN, Intercept = NaN, lambdazLower = NaN, lambdazUpper = NaN, halfLife = NaN, lambdazSpan = NaN, clastPred = NaN)
        }
    }
    return(res)
}
