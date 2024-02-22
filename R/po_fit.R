po_fit <- function(L,R,truncation = F,C,x,order = 3,equal_space = T,nknot, myknots,conv_cri = 1e-7){
  if(truncation == T){
    library(nleqslv)
    library(Matrix)
    library(ggplot2)
    if(equal_space == T){
      #the number of interior knots
      ik <- nknot-2
      #equal space knots
      #knots at the right end + 0.01 because all the basis equals to 0 at the end points if not
      mx <-  max(setdiff(c(R,L),Inf)) + 0.01
      knots <- seq(0,mx,length.out = nknot)
    }else{
      knots <- myknots
      ik <- length(knots) - 2
    }

    #Degree for Ispline
    dgr <- order
    #number of parameters
    K <- dgr + ik
    #iv is the X's
    iv <- as.matrix(x)
    P <- ncol(iv)
    n <- nrow(iv)
    #delta
    d_1 <- rep(0,n)
    d_2 <- rep(0,n)
    d_3 <- rep(0,n)
    d_0 <- as.numeric(L == R)
    for (i in which(d_0 == 0)) {
      if (L[i] == C[i]) {
        d_1[i] <- 1
      }
      else if (R[i] == Inf) {
        d_3[i] <- 1
      }
      else {
        d_2[i] <- 1
      }
    }

    Ml <- Mspline(L,dgr,knots)
    Mr <- Mspline(R,dgr,knots)
    Il <- Ispline(L,dgr,knots)
    Ir <- Ispline(R,dgr,knots)
    #rate for interval censored data
    I <- matrix(0,K,n)
    I[,d_1 == 1 | d_2 == 1] <- Ir[,d_1 + d_2 == 1] - Il[,d_1 == 1 | d_2 == 1]
    bt <- rep(1, P)
    gama <- rep(1, K)
    df <- rep(1, P)
    ite <- 0

    #left truncation time spline
    Ic <- Ispline(C,dgr,knots)

    while( t(df)%*%df > conv_cri & ite < 20000){
      #exp(x_i*beta) is n*1 vector
      exb <- exp(iv%*%bt)
      #Lambda_0(R) and Lambda_0(L), n*1 vector
      Bsr <- t(Ir)%*%gama
      Bsl <- t(Il)%*%gama

      #Ephi and Epsi are n*1 vectors
      Ephi <- as.vector(1/(Bsl*exb + 1))
      Epsi <- as.vector(1/(Bsr*exb + 1))
      Epsi[d_3 == 1] = rep(0,sum(d_3))

      #EU is a n*K matrix ( can use the ^del to get rid of the 0 denomator in right censored data)
      EU <- t(I*gama)/as.vector(t(I)%*%gama)
      EU[d_0 == 1,] <- t(Mr[,d_0 == 1]*gama)/as.vector(t(Mr[,d_0 == 1])%*%gama)
      EU[d_3 == 1,] <- rep(0,K)

      #EU is a n*K matrix ( can use the ^del to get rid of the 0 denomator in right censored data)
      EV_1 <- 1/(1 + as.vector(t(Ic)%*%gama)*exb)
      EV_2 <- (t(Ic*gama)*as.vector(exb))/as.vector(1 + as.vector(t(Ic)%*%gama)*exb)
      EV <- cbind(EV_1,EV_2)

      #The equations needed to be solved, x in the function to be solved is the beta
      #A is a K*1 vector
      A <- colSums(EV_2) + colSums(EU)
      #B is a n*K matrix
      B <- t(Ir)*Epsi*(d_0 + d_1 + d_2) + t(Il)*Ephi

      btf <- function(x){
        y <- numeric(length(bt))
        for(h in 1:length(bt)){
          y[h] <- (d_0 + d_1 + d_2 + rowSums(EV_2))%*%iv[,h] - t(B%*%(A/as.vector(t(B)%*%exp(iv%*%x)))) %*% (iv[,h]*exp(iv%*%x))
        }
        y
      }
      btstart <- rep(0,length(bt))
      #solve and get the updated bt
      sol <- nleqslv(btstart,btf,method="Newton")
      btnew <- sol$x

      gamanew <-  A/( as.vector(t(B)%*%exp(iv%*%btnew)) )

      df <- btnew - bt
      bt <- as.vector(btnew)
      gama <- as.vector(gamanew)
      ite <- ite + 1

    }

    ################calculate covariance matrix##################
    #First calculate all the expectationes and variances given the theta_hat
    ######Expectations######
    #x_i*beta is n*1 vector
    xb <- iv%*%bt
    #exp(x_i*beta) is n*1 vector
    exb <- exp(xb)
    #Lambda_0(R) and Lambda_0(L), n*1 vector
    Bsr <- t(Ir)%*%gama
    Bsl <- t(Il)%*%gama

    #Ephi and Epsi are n*1 vectors
    Ephi <- as.vector(1/(Bsl*exb + 1))
    Epsi <- as.vector(1/(Bsr*exb + 1))
    Epsi[d_3 == 1] = rep(0,sum(d_3))

    #EU is a n*K matrix ( can use the ^del to get rid of the 0 denomator in right censored data)
    EU <- t(I*gama)/as.vector(t(I)%*%gama)
    EU[d_0 == 1,] <- t(Mr[,d_0 == 1]*gama)/as.vector(t(Mr[,d_0 == 1])%*%gama)
    EU[d_3 == 1,] <- rep(0,K)

    #EU is a n*K matrix ( can use the ^del to get rid of the 0 denomator in right censored data)
    EV_1 <- as.vector( 1/(1 + as.vector(t(Ic)%*%gama)*exb) )
    EV_2 <- (t(Ic*gama)*as.vector(exb))/as.vector(1 + as.vector(t(Ic)%*%gama)*exb)


    #########variances#############
    #variance phi is n*1 vector
    Vphi <- as.vector(1/(Bsl*exb + 1))^2
    Vpsi <- as.vector(1/(Bsr*exb + 1))^2
    Vpsi[d_3 == 1] = rep(0,sum(d_3))

    #VU variance of U_il is a n*K matrix, the il_th entry is var(u_il)
    VU <- EU*(1- EU)
    VV <- EV_2*(1- EV_2)

    ####part 1 Q########################################
    Q <- matrix(0,(K+P),(K+P))
    #same A and B as before (in the EM part). A is a K*1 vector, B is a n*K matrix
    A <- colSums(EV_2) + colSums(EU)
    B <- t(Ir)*Epsi*(d_0 + d_1 + d_2) + t(Il)*Ephi

    #elements in xx^t as a vector, four entries as a group
    e <- as.vector(t(cbind(iv,iv))) * rep(as.vector(t(iv)),each = P)
    #coefficients for each matrix, repeat each elememt p times so it match the vector above
    co <-  rep(B%*%gama*exb, each = P^2)
    #fill in the corresponding place in Q
    Q[seq(1,P),seq(1,P)] <- - matrix(apply(matrix(e*co,P^2,n),1,sum),P,P)

    Q[seq(P+1,K+P),seq(1,P)] <- - t(B)%*%(iv*as.vector(exb))

    Q[seq(1,P),seq(P+1,P+K)] <- t( Q[seq(P+1,K+P),seq(1,P)])

    diag(Q[seq(1+P,P+K),seq(1+P,P+K)]) <- - gama^(-2)*A

    ####part 2 VC ###########################################
    vc <- matrix(0,(K+P),(K+P))

    #cov(l_c/beta,l_c/beta )
    e_covc <- Vpsi*(Bsr*exb)^2 + Vphi*(Bsl*exb)^2 + EV_1*(1 - EV_1)
    covc <-  rep(e_covc, each = P^2)
    vc[seq(1,P),seq(1,P)] <- matrix(apply(matrix(e*covc,P^2,n),1,sum),P,P)

    #coefficients for cov(l_c/beta,l_c/gama )
    co_t <- t(Ir)*as.vector(Bsr*exb^2)*Vpsi + t(Il)*as.vector(Bsl*exb^2)*Vphi + t(t(EV_1*EV_2)/gama)
    vc[seq(1,P),seq(P+1,P+K)] <- t(iv)%*%co_t
    vc[seq(P+1,K+P),seq(1,P)] <- t(vc[seq(1,P),seq(P+1,P+K)])

    #coefficients for cov(l_c/gama, l_c/gama) on diagnal
    dg1 <- t(VU + VV)/gama^2
    dg2 <- ( t(Ir^2)*Vpsi + Vphi*t(Il^2) )*as.vector(exb^2)

    diag(vc[seq(P+1,K+P),seq(P+1,K+P)]) <- apply(t(dg1) + dg2,2,sum)
    #coefficients for cov(l_c/gama_l, l_c/gama_k) off diagnal
    for(l in 1:K){
      for(m in 1:K){
        if(l != m){
          part_1 <- -as.vector(EV_2[,l]*EV_2[,m] + EU[,l]*EU[,m])/as.vector(gama[l]*gama[m])
          part_2 <- (Ir[l,]*Ir[m,]*Vpsi + Il[l,]*Il[m,]*Vphi)*exb^2
          vc[P+l,P+m] <- sum( part_1 + part_2 )
        }
      }
    }

    v <- -(Q + vc)

    ####part3
    tol <- 1e-7
    if(sum(gama < tol) == 0){
      vv <- v
    } else {
      rg=1:K
      index=rg[gama < tol]+P
      vv=v[-index,-index]
    }

    if( rcond(vv) > .Machine$double.eps ){
      se_theta = sqrt(diag(solve(vv))[1:P])
    } else {
      se_theta = sqrt(diag(solve(vv + diag(1e-4,nrow(vv),ncol(vv))))[1:P])
    }


    ###############CI###########################
    CI_lower <- bt - 1.96*se_theta
    CI_upper <- bt + 1.96*se_theta
    ci <- cbind(CI_lower,CI_upper)

    #define function log??likelihood AND calculate AIC and BIC
    llhd <- sum(log((1+t(Ic)%*%gama*exp(iv%*%bt))*((t(Ml)%*%gama) * exp(iv%*%bt)/( t(Il)%*%gama * exp(iv%*%bt) + 1 )^2 )^d_0 * ( 1/( t(Il)%*%gama * exp(iv%*%bt) + 1) - 1/( t(Ir)%*%gama * exp(iv%*%bt) + 1))^(d_2+d_1) * ( 1/( t(Il)%*%gama * exp(iv%*%bt) + 1) )^(d_3)))

    AIC <- 2*(P+K) - 2*llhd
    BIC <- (P+K)*log(n) - 2*llhd

    #########################################################
    # plots:odds, survival,hazard
    ########################################################
    #grids of time points
    tgrids <- seq(0,mx,0.1)
    #calculate baseline odds
    b <- Ispline(tgrids,dgr,knots)
    m <- Mspline(tgrids,dgr,knots)
    odds <- t(as.matrix(gama)) %*% b
    ff <- t(as.matrix(gama))%*%m

    #check baseline hazard rate
    hzd <- ff/(1+odds)
    bsl_hz = ggplot() + geom_line(data = as.data.frame( cbind(tgrids,t(hzd) ) ),aes(tgrids,t(hzd) )) + labs(x="t",y="h(t)")
    #check baseline odds
    bsl_odds = ggplot() + geom_line(data = as.data.frame( cbind(tgrids,t(odds) ) ),aes(tgrids,t(odds) )) + labs(x="t",y="Baseline odds")
    #check baseline survival
    sur <- 1/(1+odds)
    bsl_surv = ggplot() + geom_line(data = as.data.frame( cbind(tgrids,t(sur) ) ),aes(tgrids,t(sur) )) + labs(x="t",y="S(t)")

    #Return the result
    list(coefficient_est = as.data.frame(cbind(bt,se_theta,ci),row.names = colnames(x)), spline_coef = gama, knots = knots, AIC = AIC, BIC = BIC, Baseline_Surv = bsl_surv, Baseline_hazard = bsl_hz, Baseline_odds = bsl_odds )
  }else{
    library(nleqslv)
    library(Matrix)
    library(ggplot2)
    if(equal_space == T){
      #the number of interior knots
      ik <- nknot-2
      #equal space knots
      #knots at the right end + 0.01 because all the basis equals to 0 at the end points if not
      mx <-  max(setdiff(c(R,L),Inf)) + 0.01
      knots <- seq(0,mx,length.out = nknot)
    }else{
      knots <- myknots
      ik <- length(knots) - 2
    }

    #Degree for Ispline
    dgr <- order
    #number of parameters
    K <- dgr + ik
    #iv is the X's
    iv <- as.matrix(x)
    P <- ncol(iv)
    n <- nrow(iv)
    #delta
    d_1 <- rep(0,n)
    d_2 <- rep(0,n)
    d_3 <- rep(0,n)
    d_0 <- as.numeric(L == R)
    for (i in which(d_0 == 0)) {
      if (L[i] == 0) {
        d_1[i] <- 1
      }
      else if (R[i] == Inf) {
        d_3[i] <- 1
      }
      else {
        d_2[i] <- 1
      }
    }

    Ml <- Mspline(L,dgr,knots)
    Mr <- Mspline(R,dgr,knots)
    Il <- Ispline(L,dgr,knots)
    Ir <- Ispline(R,dgr,knots)
    #rate for interval censored data
    I <- matrix(0,K,n)
    I[,d_1 == 1 | d_2 == 1] <- Ir[,d_1 + d_2 == 1] - Il[,d_1 == 1 | d_2 == 1]

    #set an initial value
    bt <- rep(1,P)
    gama <- rep(1,K)

    ################################
    df <- rep(1,P)
    ite <- 0

    while( t(df)%*%df > conv_cri & ite < 20000){
      #exp(x_i*beta) is n*1 vector
      exb <- exp(iv%*%bt)
      #Lambda_0(R) and Lambda_0(L), n*1 vector
      Bsr <- t(Ir)%*%gama
      Bsl <- t(Il)%*%gama

      #Ephi and Epsi are n*1 vectors
      Ephi <- as.vector(1/(Bsl*exb + 1))
      Epsi <- as.vector(1/(Bsr*exb + 1))
      Epsi[d_3 == 1] = rep(0,sum(d_3))

      #EU is a n*K matrix ( can use the ^del to get rid of the 0 denomator in right censored data)
      EU <- t(I*gama)/as.vector(t(I)%*%gama)
      EU[d_0 == 1,] <- t(Mr[,d_0 == 1]*gama)/as.vector(t(Mr[,d_0 == 1])%*%gama)
      EU[d_3 == 1,] <- rep(0,K)

      #The equations needed to be solved, x in the function to be solved is the beta
      #A is a K*1 vector
      A <- colSums(EU)
      #B is a n*K matrix
      B <- t(Ir)*Epsi*(d_0 + d_1 + d_2) + t(Il)*Ephi

      btf <- function(x){
        y <- numeric(length(bt))
        for(h in 1:length(bt)){
          y[h] <- (d_0 + d_1 + d_2)%*%iv[,h] - t(B%*%(A/as.vector(t(B)%*%exp(iv%*%x)))) %*% (iv[,h]*exp(iv%*%x))
        }
        y
      }
      btstart <- rep(0,length(bt))
      #solve and get the updated bt
      sol <- nleqslv(btstart,btf,method="Newton")
      btnew <- sol$x

      gamanew <-  A/( as.vector(t(B)%*%exp(iv%*%btnew)) )

      df <- btnew - bt
      bt <- as.vector(btnew)
      gama <- as.vector(gamanew)
      ite <- ite + 1
    }

    ################calculate covariance matrix##################
    #First calculate all the expectationes and variances given the theta_hat
    ######Expectations######
    #x_i*beta is n*1 vector
    xb <- iv%*%bt
    #exp(x_i*beta) is n*1 vector
    exb <- exp(xb)
    #Lambda_0(R) and Lambda_0(L), n*1 vector
    Bsr <- t(Ir)%*%gama
    Bsl <- t(Il)%*%gama

    #Ephi and Epsi are n*1 vectors
    Ephi <- as.vector(1/(Bsl*exb + 1))
    Epsi <- as.vector(1/(Bsr*exb + 1))
    Epsi[d_3 == 1] = rep(0,sum(d_3))

    #EU is a n*K matrix ( can use the ^del to get rid of the 0 denomator in right censored data)
    EU <- t(I*gama)/as.vector(t(I)%*%gama)
    EU[d_0 == 1,] <- t(Mr[,d_0 == 1]*gama)/as.vector(t(Mr[,d_0 == 1])%*%gama)
    EU[d_3 == 1,] <- rep(0,K)

    #########variances#############
    #variance phi is n*1 vector
    Vphi <- as.vector(1/(Bsl*exb + 1))^2
    Vpsi <- as.vector(1/(Bsr*exb + 1))^2
    Vpsi[d_3 == 1] = rep(0,sum(d_3))

    #VU variance of U_il is a n*K matrix, the il_th entry is var(u_il)
    VU <- EU*(1- EU)

    ####part 1 Q########################################
    Q <- matrix(0,(K+P),(K+P))
    #same A and B as before (in the EM part). A is a K*1 vector, B is a n*K matrix
    A <- colSums(EU)
    B <- t(Ir)*Epsi*(d_0 + d_1 + d_2) + t(Il)*Ephi

    #elements in xx^t as a vector, four entries as a group
    e <- as.vector(t(do.call(cbind,replicate(P,iv,simplify = F)))) * rep(as.vector(t(iv)),each = P)
    #coefficients for each matrix, repeat each elememt p times so it match the vector above
    co <-  rep(B%*%gama*exb, each = P^2)
    #fill in the corresponding place in Q
    Q[seq(1,P),seq(1,P)] <- - matrix(apply(matrix(e*co,P^2,n),1,sum),P,P)

    Q[seq(P+1,K+P),seq(1,P)] <- - t(B)%*%(iv*as.vector(exb))

    Q[seq(1,P),seq(P+1,P+K)] <- t( Q[seq(P+1,K+P),seq(1,P)])

    diag(Q[seq(1+P,P+K),seq(1+P,P+K)]) <- - gama^(-2)*A

    ####part 2 VC ###########################################
    vc <- matrix(0,(K+P),(K+P))

    #cov(l_c/beta,l_c/beta )
    e_covc <- Vpsi*(Bsr*exb)^2 + Vphi*(Bsl*exb)^2
    covc <-  rep(e_covc, each = P^2)
    vc[seq(1,P),seq(1,P)] <- matrix(apply(matrix(e*covc,P^2,n),1,sum),P,P)

    #coefficients for cov(l_c/beta,l_c/gama )
    co_t <- t(Ir)*as.vector(Bsr*exb^2)*Vpsi + t(Il)*as.vector(Bsl*exb^2)*Vphi
    vc[seq(1,P),seq(P+1,P+K)] <- t(iv)%*%co_t
    vc[seq(P+1,K+P),seq(1,P)] <- t(vc[seq(1,P),seq(P+1,P+K)])

    #coefficients for cov(l_c/gama, l_c/gama) on diagnal
    dg1 <- t(VU)/gama^2
    dg2 <- ( t(Ir^2)*Vpsi + Vphi*t(Il^2) )*as.vector(exb^2)

    diag(vc[seq(P+1,K+P),seq(P+1,K+P)]) <- apply(t(dg1) + dg2,2,sum)
    #coefficients for cov(l_c/gama_l, l_c/gama_k) off diagnal
    for(l in 1:K){
      for(m in 1:K){
        if(l != m){
          part_1 <- -as.vector(EU[,l]*EU[,m])/as.vector(gama[l]*gama[m])
          part_2 <- (Ir[l,]*Ir[m,]*Vpsi + Il[l,]*Il[m,]*Vphi)*exb^2
          vc[P+l,P+m] <- sum( part_1 + part_2 )
        }
      }
    }

    v <- -(Q + vc)

    ####part3
    tol <- 1e-7
    if(sum(gama < tol) == 0){
      vv <- v
    } else {
      rg=1:K
      index=rg[gama < tol]+P
      vv=v[-index,-index]
    }

    if( rcond(vv) > .Machine$double.eps ){
      se_theta = sqrt(diag(solve(vv))[1:P])
    } else {
      se_theta = sqrt(diag(solve(vv + diag(1e-4,nrow(vv),ncol(vv))))[1:P])
    }

    ###############CI###########################
    CI_lower <- bt - 1.96*se_theta
    CI_upper <- bt + 1.96*se_theta
    ci <- cbind(CI_lower,CI_upper)

    #define function log??likelihood AND calculate AIC and BIC
    llhd <- sum(log(((t(Ml) %*% gama) * exp(iv %*% bt)/(t(Il) %*%
                                                          gama * exp(iv %*% bt) + 1)^2)^d_0 * (1 - 1/(t(Ir) %*%
                                                                                                        gama * exp(iv %*% bt) + 1))^(d_1) * (1/(t(Il) %*% gama *
                                                                                                                                                  exp(iv %*% bt) + 1) - 1/(t(Ir) %*% gama * exp(iv %*%
                                                                                                                                                                                                  bt) + 1))^(d_2) * (1/(t(Il) %*% gama * exp(iv %*% bt) +
                                                                                                                                                                                                                          1))^(d_3)))
    AIC <- 2*(P+K) - 2*llhd
    BIC <- (P+K)*log(n) - 2*llhd

    #########################################################
    # plots:odds, survival,hazard
    ########################################################
    #grids of time points
    tgrids <- seq(0,mx,0.1)
    #calculate baseline odds
    b <- Ispline(tgrids,dgr,knots)
    m <- Mspline(tgrids,dgr,knots)
    odds <- t(as.matrix(gama)) %*% b
    ff <- t(as.matrix(gama))%*%m

    #check baseline hazard rate
    hzd <- ff/(1+odds)
    bsl_hz = ggplot() + geom_line(data = as.data.frame( cbind(tgrids,t(hzd) ) ),aes(tgrids,t(hzd) )) + labs(x="t",y="h(t)")
    #check baseline odds
    bsl_odds = ggplot() + geom_line(data = as.data.frame( cbind(tgrids,t(odds) ) ),aes(tgrids,t(odds) )) + labs(x="t",y="Baseline odds")
    #check baseline survival
    sur <- 1/(1+odds)
    bsl_surv = ggplot() + geom_line(data = as.data.frame( cbind(tgrids,t(sur) ) ),aes(tgrids,t(sur) )) + labs(x="t",y="S(t)")

    #Return the result
    list(coefficient_est = as.data.frame(cbind(bt,se_theta,ci),row.names = colnames(x)), spline_coef = gama, knots = knots, AIC = AIC, BIC = BIC, Baseline_Surv = bsl_surv, Baseline_hazard = bsl_hz, Baseline_odds = bsl_odds )
  }

}
