

DR.est <- function(df,  treatment, t1, response, X, support = 0.99, alpha = 0.05, r1 = F){
    full.n   <- nrow(df)
    miss     <- as.logical(rowSums(is.na(df[,c(treatment, response, X)])))
    if(any(miss)) df <- df[!miss,]
    
    if(length(table(df[,treatment])) != 2) stop('Treatment variable must be binary.')
    df[, treatment] <- ifelse(df[,treatment] == t1, 1, 0)
    
    modelZ <- glm(reformulate(X, treatment), data = df, family = binomial)
    e      <- fitted(modelZ, df)
    trim   <- e < (1 - support)/2 | e > (1 + support)/2
    if(any(trim)){
        df <- df[!trim,]
        e  <- e[!trim]
    }
    
    n.ylev <- length(table(df[,response]))
    if(n.ylev == 2){
        cat('Response assumed binary.','\n')
        if(r1 == F & !is.numeric(df[, response])){
            r1 <- names(table(df[, response]))[1]
            cat(paste0('Using indicator of response = ', r1,'.'),'\n')
        }
        df[, response] <- ifelse(df[, response] == r1, 1, 0)
        fam <- 'binomial'
    } else if(!is.numeric(df[, treatment])){
        stop('Response is character/factor with !=2 levels. Response must be binary or continuous.')
    } else {
        cat('Response assumed continuous.','\n')
        fam <- 'gaussian'
    }
    trt <- df[, treatment] == 1
    model0 <- glm(reformulate(X, response), data = df[!trt, ], family = get(fam))
    model1 <- glm(reformulate(X, response), data = df[trt, ],  family = get(fam))
    
    m0 <- predict(model0, df) 
    m1 <- predict(model1, df)  
    e0 <- e[!trt]
    e1 <- e[trt]
    Z  <- df[, treatment]
    Y  <- df[, response]
    if(fam == 'binomial'){
        m0 <- plogis(m0)
        m1 <- plogis(m1)
    }
    
    d0    <- ((1 - Z)*Y + (Z - e)*m0)/(1 - e)
    d1    <- (Z*Y - (Z - e)*m1)/e 
    mu0   <- mean(d0)
    mu1   <- mean(d1)
    delta <- mu1 - mu0
    SE    <- sqrt(sum(((d1 - d0) - delta)^2))/nrow(df)
    p     <- 2*min(pnorm(c(-1,1)*delta/SE))
    conf  <- delta + c(-1,1)*qnorm(1 - alpha/2)*SE
    names(conf) <- c('lower','upper')
    
    est <- round(c(mu0 = mu0, mu1 = mu1, delta = delta, conf, SE = SE, p = p, alpha = alpha), 4)
    print(est); cat('\n')
    rem <- c(missing = sum(miss), extreme.ps = sum(trim))
    rem <- rbind(n = rem, '%' = round(100*rem/full.n, 1))
    cat('Removed Observations:','\n')
    print(rem)
    
    return(list(estimates = est, removed = rem, model0 = model0, model1 = model1, 
                modelZ = modelZ, e0 = e0, e1 = e1))
}


