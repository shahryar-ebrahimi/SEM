
library(lavaan)
library(Matrix) # Contains the bdiag function
mydata <- read.csv(file.choose(new = FALSE))

# Part 1: Combining them all together ------------------------------------

start <- rep(1.5, 16)
start[14:16] <- 0 #If covariances are equal to variances, matrix will not be invertible

semCustom <- function(start, dat, wish = FALSE, estimator = "ML", parameterNumb){
  
  # functions -----------------------------------
  #  
  
  # calculating cost function
  setFunction <- function(start, dat, wishart = wish, fit = estimator){
    
    #Λ is our loadings matrix. This is iteratively changed to minimize our fit function
    #We break lambda out for exogenous (lambda_x) and endogenous (lambda_y) variables
    lambda_y <- bdiag(start[1:4])
    lambda_x <- bdiag(start[5:6])
    
    # The Beta matrix contains the regression coefficients between endogenous variables.This is iteratively changed to minimize our fit function
    beta <- t(t(0))
    
    # The Gamma matrix contains the regression coefficients from exogenous to endogenous variables. This is iteratively changed to minimize our fit function
    gamma <- t(t(start[7])) #exogenous regression coefficients
    
    theta <- diag(start[8:13],ncol=6,nrow=6)
    theta[3,2] <- start[14]
    theta[5,4] <- start[15]
    theta[6,4] <- start[16]
    theta[upper.tri(theta)] <- t(theta)[upper.tri(theta)]
    
    psi <- diag(1)
    phi <- diag(1)
    
    # Get the total number of factors
    nfactors <- ncol(lambda_x)+ncol(lambda_y)
    # P is the total number of observed variables
    p <- ncol(dat)
    # N is the sample size
    N <- nrow(dat)
    
    theta_epsilon <- theta[1:nrow(lambda_y),1:nrow(lambda_y)]
    theta_delta <- theta[(p-nrow(lambda_x)+1):p,(p-nrow(lambda_x)+1):p]
    
    
    #Sigma is the model-implied covariance matrix; it's a supermatrix of 4 smaller matrices
    sigma_YY <- lambda_y %*% solve(diag(1,nrow=nrow(beta),ncol=ncol(beta))-beta) %*%((gamma%*%phi%*%t(gamma)+psi)%*% t(solve(diag(1,nrow=nrow(beta),ncol=ncol(beta))-beta))%*%t(lambda_y))+theta_epsilon
    sigma_XY <- lambda_x %*% phi %*% t(gamma) %*% t(solve(diag(1,nrow=nrow(beta),ncol=ncol(beta))-beta)) %*% t(lambda_y)
    sigma_XX <- lambda_x %*% phi %*% t(lambda_x) + theta_delta
    sigma <- cbind(rbind(sigma_YY,sigma_XY),rbind(t(sigma_XY), sigma_XX))
    
    if (wishart==TRUE){
      S <- cov(dat)
    } else{
      S <- cov(dat)*(nrow(dat)-1)/(nrow(dat))
    }
    
    
    #(Negative) Log likelihood
    tmp <- matrix(array(S%*%solve(sigma)), p ,p)
    LL0 <- (-N*p/2)*log(2*pi) - (N/2)*log(det(sigma)) - (N/2)*(sum(diag(tmp)))
    LL1 <- (-N*p/2)*log(2*pi) - (N/2)*log(det(S)) - (N/2)*p
    #nlminb seems to give an error when we use the negative log liklihood, but works with the positie log likelihood. Weird.
    
    LL0 <- -LL0
    LL1 <- -LL1
    
    
    if (fit == "ML"){
      #Maximum Likelihood Fit Function
      fit <- log(det(sigma)) - log(det(S)) + sum(diag(solve(sigma)%*%S)) - p
    } else if (fit == "ULS") {
      #Unweighted Least Squares Fit Function
      tmp <- S-sigma
      fit <- 0.5*sum(tmp^2)
    } else if (fit == "GLS") {
      #Generalized Least Squares Fit Function
      tmp <- (S-sigma)%*%solve(sigma)
      fit <- 0.5*sum(tmp^2)
    }else if (fit == "LL0") {
      fit <- LL0
    }
    
    
    #Standard errors
    #se <-sqrt(diag(solve(lavInspect(fit2, "information.expected")*75)))
    
    
    return(fit)
    
  }
  
  # calculating model dfs
  DF <- function(start, dat){    
    
    # P is the total number of observed variables
    p <- ncol(dat)
    
    Basedf <- p*(p+1)/2
    Nulldf <- Basedf-p
    Modeldf <- Basedf-length(start)
    
    return(c(Basedf, Nulldf, Modeldf))
    
  }
  
  # calculating chi-square tests
  ChiTest <- function(start, dat, fit, wishart = wish){
    
    if (wishart==TRUE){
      S <- cov(dat)
      chisqmodel <- fit*(nrow(dat)-1)
      chisqnull <- (log(det(diag(diag(S))))-log(det(S))+sum(ncol(dat))-ncol(dat))*(nrow(dat)-1)
    } else{
      S <- cov(dat)*(nrow(dat)-1)/(nrow(dat))
      chisqmodel <- fit*nrow(dat)
      chisqnull <- (log(det(diag(diag(S))))-log(det(S))+sum(ncol(dat))-ncol(dat))*nrow(dat)
    }
    
    return(c(chisqmodel, chisqnull))
    
  }
  
  # calculating fit indices
  Findices <- function(start, dat, fit, wishart = wish){
    
    chiout      <- ChiTest(start, dat, fit, wishart = wishart)
    dfout       <- DF(start, dat)
    
    chisqmodel <- chiout[1]
    chisqnull  <- chiout[2]
    
    Basedf     <- dfout[1]
    Nulldf     <- dfout[2]
    Modeldf    <- dfout[3]
    
    N          <- nrow(dat)
    
    LL0        <- setFunction(start, dat, wishart = wishart, fit = "LL0")
    
    #Fit indices
    CFI <- 1- (chisqmodel-Modeldf)/(chisqnull-Nulldf)
    # Tucker-Lewis Index/Non-Normed Fit Index
    TLI <- ((chisqnull/Nulldf)-(chisqmodel/Modeldf))/((chisqnull/Nulldf)-1)
    # Bentler-Bonnet Normed Fit Index
    BBNFI <- 1 - (chisqmodel/chisqnull)
    # Parsimony Normed Fit Index
    PNFI <- (Modeldf/Nulldf)*BBNFI
    # Bollen's Relative Fit Index
    BRFI <- ((chisqnull/Nulldf)-(chisqmodel/Modeldf))/(chisqnull/Nulldf)
    # Bollen's Incremental Fit Index
    BIFI <- (chisqnull-chisqmodel)/(chisqnull-Modeldf)
    #Relative Noncentrality Index
    RNI <- ((chisqnull-Nulldf) - (chisqmodel-Modeldf))/(chisqnull-Nulldf)
    
    #RMSEA
    RMSEA <- sqrt((chisqmodel-Modeldf)/(Modeldf*N))
    
    #AIC
    AIC <- (-2*(-LL0))+2*Nulldf
    #BIC
    BIC <- (-2*(-LL0))+Nulldf*log(N)
    
    #SSABIC
    SSABIC <- (-2*(-LL0))+Nulldf*log((N+2)/24)
    #ECVI
    ECVI <- fit +(2*(Basedf-Modeldf))/N
    
    #Hoelter's Critical N, alpha=.05
    HCN05 <- qchisq(.05, Modeldf, lower.tail=FALSE)/fit+1
    
    #Hoelter's Critical N, alpha=.01
    HCN01 <- qchisq(.01, Modeldf, lower.tail=FALSE)/fit+1
    
    #McDonald's Fit Index
    MFI <- exp(-.5*((chisqmodel-Modeldf)/N))
    
    return(c(SSABIC, BIC, AIC, RMSEA, RNI, BIFI, BRFI, PNFI, BBNFI, TLI, CFI, ECVI, HCN05, HCN01, MFI))
    
  }
  
  # printing result
  printResult <- function(output.estimate, output.df, output.chi, output.indices, parameterNumb){
    
    parameterName = c("loading", "Gamma", "Variance", "Covariance")
    IndicesName   = c("SSABIC", "BIC", "AIC", "RMSEA", "RNI", "BIFI", "BRFI", "PNFI", "BBNFI", "TLI", "CFI", "ECVI", "HCN05", "HCN01", "MFI")
    
    writeLines(sprintf(""))
    writeLines(sprintf("--------------------------------"))
    
    for(rep in 1:length(parameterNumb)){
      
      if (rep==1){
        rs <- 1
        re <- parameterNumb[1]
      }else{
        rs <- re+1
        re <- re+parameterNumb[rep]
      }
      
      writeLines(sprintf("%s %s = %f", parameterName[rep], 1:parameterNumb[rep], output.estimate$par[rs:re]))
      writeLines(sprintf("--------------------------------"))
      
    }  
    
    writeLines(sprintf("Base  DF = %f", output.df[1]))
    writeLines(sprintf("Null  DF = %f", output.df[2]))
    writeLines(sprintf("Model DF = %f", output.df[3]))
    
    writeLines(sprintf("--------------------------------"))
    
    writeLines(sprintf("Model Chi-Sq = %f", output.chi[1]))
    writeLines(sprintf("Null  Chi-Sq = %f", output.chi[2]))
    
    writeLines(sprintf("--------------------------------"))
    
    for (rep in 1:length(IndicesName)){
      writeLines(sprintf("%s  = %f", IndicesName[rep], output.indices[rep]))
    }
    
    writeLines(sprintf("--------------------------------"))
    writeLines(sprintf(""))
    
  } 
  
  # Outputs -------------------------------------
  # 
  
  # saving outputs
  output.estimate <- nlminb(start, setFunction, dat=dat)
  output.df       <- DF(start, dat)
  output.chi      <- ChiTest(output.estimate$par, dat, fit = output.estimate$objective, wishart = wish)
  output.indices  <- Findices(output.estimate$par, dat, fit = output.estimate$objective, wishart = wish)
  
  # forcing first loading to 1
  
  #output.estimate$par[1:4] <- output.estimate$par[1:4]/output.estimate$par[1]
  #output.estimate$par[5:6] <- output.estimate$par[5:6]/output.estimate$par[5]
  
  
  printResult(output.estimate, output.df, output.chi, output.indices, parameterNumb)
  
  return(list(output.estimate, output.df, output.chi, output.indices))
  
}

res <- semCustom(start, dat = mydata, wish = FALSE, estimator = "ULS", parameterNumb = c(6, 1, 6, 3))

