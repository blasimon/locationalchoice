
#======================  Personal Space Project =============================================

### MODEL ESTIMATION

# Heterogeneous estimation
# Allows for:
# - different number of choices per respondent
# - continous X and/or dummy coded X
# - different model asumptions (see "flag" values)
# - covariates in PS and L components (same covariates)

library(bayesm)
library(Rcpp)

#dyn.load("C:/Users/Tanya/Desktop/Personal space/Cpp optimization/Release/PS_locational.dll")  # home
#dyn.load("C:/Users/Tanya/Desktop/Personal space/Cpp optimization/Debug/PS_locational.dll")  # home
#dyn.load("C:/Users/Tanya/Desktop/Personal space/Cpp optimization/Release/PS_locational_withBeta.dll")  # home
#dyn.load("C:/Users/user/Desktop/PS/DLLs/PS_locational_universal.dll")   # office with Cov
#dyn.load("C:/Users/user/Desktop/PS/DLLs/PS_locational_311.dll")   # office
#dyn.load("C:/Users/user/Desktop/PS/DLLs/PS_locational_new.dll")   # office
#dyn.unload("C:/Users/user/Desktop/PS/DLLs/PS_locational - withBeta.dll")   # office

# Renamed to reflect the regression nature, and not mixture

FUN_PS_LPSreg_het_Cpp_withCov_diffNumCh = function (Data, Prior, Mcmc)
{
  #  flags from C++
  #const int Beta3 = 1; //locational curved 3-parameter
	#const int Alpha = 2; //personal space
	#const int Beta2 = 4; //locational curved 2-parameter (Kumaraswammy)
	#const int Beta1 = 8; //locational linear

  #  flag values in R
  # 1 = locational quadratic
  # 2 = PS only
  # 3 = locational quadratic + PS
  # 4 = locational 2 par only (Kumaraswammy)
  # 6 = locational 2 par (Kumaraswammy) + PS
  # 8 = locational linear only
  # 10 = locational linear + PS

  # Load data
    lgtdata = Data$lgtdata
    nlgt = length(lgtdata)
    #Z = Data$Z
    Z = matrix(rep(1, nlgt), ncol = 1)
    nz = ncol(Z)

  # Priors 
    nu_b = Prior$nu_b
    V_b = Prior$V_b
    ADelta = Prior$ADelta
    Deltabar = Prior$Deltabar
   
  # MCCM settings
    sbeta_vec = Mcmc$sbeta_vec
    keep = Mcmc$keep
    R = Mcmc$R
    set.seed(Mcmc$seed)
    switchR = ifelse(is.null(Mcmc$switchR)==T,1000, Mcmc$switchR)
    flag = as.integer(Mcmc$flag)
    #PSflag = as.integer(Mcmc$PSflag) 
    PS_var = Mcmc$PS_var
    L_var = Mcmc$L_var
    Cond_var = Mcmc$Cond_var
    
    #nvar_cond = length(lgtdata[[1]][[1]]$X_cond[Cond_var])      # XcondL and XcondP is the same for now
    # use if have covariates
    nvar_cond = 1     # XcondL and XcondP is the same for now

    # this number of parameters that is needed for C, not the actual number of parameters
    nvar_PS =  length(lgtdata[[1]][[1]]$X_PS2[1,PS_var])
    nvar_L =   ifelse(flag==8|flag==10,length(lgtdata[[1]][[1]]$X_L[1,L_var]),     # linear
               #ifelse(flag==8|flag==10,ncol(lgtdata[[1]][[1]]$X_L-1),             # hard coded to remove col=18 (var 29) for now
               ifelse(flag==4|flag==6,ncol(lgtdata[[1]][[1]]$X_L)+2,               # curved 2 par original line (Kuma)
               #ifelse(flag==4|flag==6,ncol(lgtdata[[1]][[1]]$X_L)+4,              # curved 2 par new with no b3=b1b2 and b0x 
               ifelse(flag==1|flag==3,length(lgtdata[[1]][[1]]$X_L[1,L_var])+3,    # quadratic locational
               ifelse(flag==2,length(lgtdata[[1]][[1]]$X_PS2[1,PS_var]),            # PS
                 1))))
    npar = ifelse(flag==3|flag==6|flag==10,nvar_L+nvar_PS,   # combined L+PS
           ifelse(flag==2,nvar_PS,                           # PS only
                          nvar_L))                           # L only
    
    betadraw = array(double(floor(R/keep) * nlgt * npar * nvar_cond ), dim = c(nlgt,npar*nvar_cond, floor(R/keep)))
    Vbetadraw = matrix(double(floor(R/keep) * npar * npar * nvar_cond* nvar_cond ), ncol = npar * npar * nvar_cond* nvar_cond )
    Deltadraw = matrix(double(floor(R/keep) * npar * nvar_cond  * nz), ncol = npar * nz * nvar_cond )
    
    oldbetas = matrix(double(nlgt * npar  * nvar_cond ), ncol = npar  * nvar_cond  )
    #oldbetas = matrix((runif(nlgt * npar,-0.1,0.1)), ncol = npar)
    newbetas = matrix(double(nlgt * npar  * nvar_cond ), ncol = npar * nvar_cond  )
    oldVbeta = diag(npar * nvar_cond  )
    oldVbetai = diag(npar  * nvar_cond )
    oldDelta = matrix(double(npar * nz  * nvar_cond ), ncol = npar  * nvar_cond )

    #oldbetas = betastart
    #oldVbeta = 5 * diag((npar* nvar_cond))
    #oldVbetai = chol2inv(chol(oldVbeta))
    #oldDelta = deltastart
    

    reject = matrix(0,ncol=1,nrow=floor(R/keep))
    llike = array(0, dim = c(R/keep))
    logl = matrix(0,ncol=1,nrow=nlgt)
    dummy = double(1)

    CData = integer(2)
    CData = .C("Initialize", CData=CData, as.integer(nlgt), as.integer(nvar_L), as.integer(ncol(lgtdata[[1]][[1]]$X_L[,L_var])), as.integer(nvar_PS))$CData  
    #CData = .C("Initialize", CData=CData, as.integer(nlgt), as.integer(nvar_L+nvar_PS), as.integer(ncol(lgtdata[[1]][[1]]$X_L+ncol(lgtdata[[1]][[1]]$X_PS))), as.integer(nvar_PS), as.integer(2))$CData  

    # wC are parameters that go into C after accounting for the conditions/covariates
    wC_all_old = list()  #  code for different # of choices for all respondents
    X_cond = list()
    nch_hh = c(rep(0,nlgt))
    wC_old_table = NULL
    for(hh in 1:nlgt)                                                                                                   
    {   nch_hh[hh] = length(lgtdata[[hh]])
        XL_hh = NULL
        XP_hh = NULL
        y_hh = c(rep(0,nch_hh[hh]))
        Xsize_hh = c(rep(0,nch_hh[hh]))
        X_cond[[hh]] = matrix(0,nrow=nch_hh[hh],ncol=nvar_cond)
        for(j in 1:nch_hh[hh])
         {  XtempL=lgtdata[[hh]][[j]]$X_L[,L_var]  # need to point to variables in PS in vector L_var
            XL_hh = rbind(XL_hh,XtempL)  # original code
            XP_hh = rbind(XP_hh,matrix(lgtdata[[hh]][[j]]$X_PS2[,PS_var],ncol=nvar_PS) ) # need to point to variables in PS in vector PS_var
            Xsize_hh[j] = as.integer(nrow(XtempL))
            y_hh[j] = as.integer(lgtdata[[hh]][[j]]$y)
            if(is.null(filetouse[[hh]][[1]]$X_cond))
            {  X_cond[[hh]][j,] = c(rep(1,nvar_cond))
              }
            else{  X_cond[[hh]][j,] = lgtdata[[hh]][[j]]$X_cond[Cond_var]  # need to point to variables in X_cond in vector Cond_var
            }
         }
        XP_hh = matrix(as.double(XP_hh),ncol=nvar_PS) 
        .C("StoreX", CData, t(matrix(as.numeric(XL_hh),ncol=length(L_var))), t(XP_hh), as.integer(Xsize_hh), as.integer(y_hh), as.integer(hh), as.integer(nch_hh[hh]))    
        #.C("StoreX", CData, t(matrix(as.numeric(XL_hh),ncol=nvar_L)), t(XP_hh), as.integer(Xsize_hh), as.integer(y_hh), as.integer(hh), as.integer(nch_hh[hh]))  
        #.C("StoreX", CData, t(matrix(as.numeric(cbind(XL_hh,XP_hh)),ncol=nvar_L+nvar_PS)), t(XP_hh), as.integer(Xsize_hh), as.integer(y_hh), as.integer(hh), as.integer(nch_hh[hh])) 
        #.C("StoreX", CData, t(matrix(as.numeric(XL_hh),ncol=4)), t(XP_hh), as.integer(Xsize_hh), as.integer(y_hh), as.integer(hh), as.integer(nch_hh[hh]))    
        wC_all_old[[hh]] = X_cond[[hh]] %*%  t( matrix(oldbetas[hh,],ncol=nvar_cond, nrow=npar, byrow=T))
        wC_old_table = rbind(wC_old_table,wC_all_old[[hh]])
    }
    
  #  flag values in R
  # 1 = locational quadratic
  # 2 = PS only
  # 3 = locational quadratic + PS
  # 4 = locational 2 par only (Kumaraswammy)
  # 6 = locational 2 par (Kumaraswammy) + PS
  # 8 = locational linear only
  # 10 = locational linear + PS

    if(flag==1|flag==8)      {  betaC =  wC_old_table[,1:(nvar_L)]         # L quadratic (or linear)
                                alphaC = dummy
    }    
    if(flag==3|flag==10)     {  betaC =  wC_old_table[,1:(nvar_L)]             # combination of L quadratic (or linear) and PS 
                                alphaC = wC_old_table[,-1:-(nvar_L)]
    }
    if(flag==4)     {  betaC = exp( wC_old_table[,1:(nvar_L)])   # L kuma  
                       alphaC = dummy
    }    
    if(flag==6)     {  betaC = exp(wC_old_table[,1:(nvar_L)])    # combination of L Kuma + PS  
                       alphaC = wC_old_table[,-1:-(nvar_L)]
    }
    if(flag==2 )    {  betaC = dummy                       # PS only
                       alphaC = wC_old_table
    } 
    
    ll_out_cpp = .C("LogLikelihood", CData, t(betaC), t(alphaC), ll=logl, as.integer(flag), PredictiveFlag=0)
    logl = ll_out_cpp$ll

    itime = proc.time()[3]
    cat("MCMC Iteration (est time to end - min)", fill = TRUE)
    #fsh()

    wC_new_table = wC_old_table
    
    for (r in 1:R) {
       
       # locational and PS betas
       rej_b = 0
      if(1) {
       sbeta=sbeta_vec[1]
       if(r>switchR) {sbeta=sbeta_vec[2]}
       sV = sbeta * oldVbeta
       root = t(chol(sV))
       #wC_new_table = NULL
       for (i in 1:nlgt) {
           newbetas[i,] = oldbetas[i, ] + root %*% rnorm(npar*nvar_cond)
           wC_new_table[(sum(nch_hh[1:i])-sum(nch_hh[i])+1):sum(nch_hh[1:i]),] = X_cond[[i]] %*%  t( matrix(newbetas[i,],ncol=nvar_cond, nrow=npar, byrow=T))    # 5 times faster
           #wC_new_table= rbind(wC_new_table, X_cond[[i]] %*%  t( matrix(newbetas[i,],ncol=nvar_cond, nrow=npar, byrow=T)))  
       }
       if(flag==1|flag==8)      {  betaC = wC_new_table[,1:nvar_L]         # L quadratic (or linera)
                                   alphaC = dummy
       }    
       if(flag==3|flag==10)     {  betaC = (wC_new_table[,1:nvar_L])       # combination of L quadratic (or linear) and PS 
                                   alphaC = wC_new_table[,-1:-nvar_L]
       }
       if(flag==4)     {  betaC = exp( wC_new_table[,1:nvar_L])   # L kuma 
                          alphaC = dummy
       }    
       if(flag==6)     {  betaC = exp(wC_new_table[,1:nvar_L])    # combination of L Kuma + PS
                          alphaC = wC_new_table[,-1:-nvar_L]
       }
       if(flag==2 )    {  betaC = dummy                       # PS only
                          alphaC = wC_new_table
       } 
       
       ll_out_cpp = .C("LogLikelihood", CData, t(betaC), t(alphaC), ll=logl, as.integer(flag), PredictiveFlag=0)
       #ll_out_cpp$ll
       for (i in 1:nlgt) {
          lognew = ll_out_cpp$ll[i]
          logold = logl[i]
          logknew = -0.5 * (t(newbetas[i,]) - Z[i, ] %*% oldDelta) %*% oldVbetai %*% (newbetas[i,] - t(Z[i, ] %*% oldDelta))
          logkold = -0.5 * (t(oldbetas[i,]) - Z[i, ] %*% oldDelta) %*% oldVbetai %*% (oldbetas[i,] - t(Z[i, ] %*% oldDelta))
            alpha_r = exp(lognew + logknew - logold - logkold)
            if (alpha_r == "NaN")
                alpha_r = -1
            u = runif(n = 1, min = 0, max = 1)
            if (u < alpha_r) {  oldbetas[i, ] = newbetas[i,]
                                logl[i] =  lognew }
            else             {  rej_b = rej_b + 1      }
        }
       }
       out = rmultireg(oldbetas, Z, Deltabar, ADelta, nu_b, V_b)
       oldDelta = out$B
       oldVbeta = out$Sigma
       oldVbetai = chol2inv(chol(oldVbeta))
        
       current_rej = c(rej_b/nlgt)
       current_ll = sum(logl)

       mkeep = r/keep
       if (mkeep * keep == (floor(mkeep) * keep)) {
          Deltadraw[mkeep, ] = as.vector(oldDelta)
          Vbetadraw[mkeep, ] = as.vector(oldVbeta)
          betadraw[, ,mkeep ] = oldbetas
          llike[mkeep]= current_ll
          reject[mkeep] = current_rej
        }
    if (r%%keep == 0) {
    #if ((r%%(keep/10)) == 0) {
            ctime = proc.time()[3]
            timetoend = ((ctime - itime)/r) * (R - r)
            cat(" ", r, " (", round(timetoend/60, 1), ")", "LL= ",current_ll,
                "RR= ", format(reject[mkeep],digits=3),"D= ",format(oldDelta,digits=1), fill = TRUE)
            #fsh()
        }    
    }
    ctime = proc.time()[3]
    cat(" Total Time Elapsed:", round((ctime - itime)/60, 2),fill = TRUE)
    .C("Cleanup", CData)
    return(list(betadraw = betadraw, Vbetadraw = Vbetadraw, Deltadraw = Deltadraw,
                llike = llike, reject = reject,flag=flag,switchR=switchR,
                seed=Mcmc$seed,ADelta = ADelta,nu_b=nu_b, V_b=V_b,Deltabar=Deltabar,
                keep=keep,R=R,sbeta_vec=sbeta_vec))
}




