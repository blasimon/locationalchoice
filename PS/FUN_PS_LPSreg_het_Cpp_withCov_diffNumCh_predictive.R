#======================  Personal Space Project =============================================

### MODEL

# Can be used for both continous X and dummy coded X
# Created:        2016-10-07


library(bayesm)
library(Rcpp)

#dyn.unload("C:/Users/Tanya/Desktop/Personal space/Cpp optimization/Debug/PS_locational.dll")
#dyn.load("C:/Users/Tanya/Desktop/Personal space/Cpp optimization/Release/PS_locational.dll")    # home
#dyn.load("C:/Users/user/Desktop/PS/DLLs/PS_locational.dll")

# Renamed to reflect the regression nature, and not mixture

FUN_PS_LPSreg_het_Cpp_withCov_diffNumCh_predictive = function (Data, Mcmc)
{
  #  flags from C++
  #const int Beta3 = 1; //locational curved 3-parameter
  #const int Alpha = 2; //personal space
  #const int Beta2 = 4; //locational curved 2-parameter
  #const int Beta1 = 8; //locational linear

  #  flag values in R
  # 1 = locational 3 par only
  # 2 = PS only
  # 3 = locational 3 par + PS
  # 4 = locational 2 par only
  # 6 = locational 2 par + PS
  # 8 = locational linear only
  # 10 = locational linear + PS

    lgtdata = Data$lgtdata
    nlgt = length(lgtdata)
    #Z = Data$Z
    Z = matrix(rep(1, nlgt), ncol = 1)
    nz = ncol(Z)
    flag = as.integer(Mcmc$flag)
    PredictiveFlag = as.integer(Mcmc$PredictiveFlag)
    oldbetas = Mcmc$startBeta                      # means of posterior draws
    PSflag = as.integer(Mcmc$PSflag) 
    dummy = double(1)
    PS_var = Mcmc$PS_var
    L_var = Mcmc$L_var
    Cond_var = Mcmc$Cond_var
    
    if(is.null(lgtdata[[1]][[1]]$X_cond[Cond_var]))
      { nvar_cond = 1 }    # XcondL and XcondP is the same for now
    else{ nvar_cond = length(lgtdata[[1]][[1]]$X_cond[Cond_var])      # XcondL and XcondP is the same for now
    }
    
    # this number of parameters that is needed for C, not the actual number of parameters
    nvar_PS =  length(lgtdata[[1]][[1]]$X_PS2[1,PS_var])
    nvar_L =   ifelse(flag==8|flag==10,length(lgtdata[[1]][[1]]$X_L[1,L_var]),     # linear
               #ifelse(flag==8|flag==10,ncol(lgtdata[[1]][[1]]$X_L-1),             # hard coded to remove col=18 (var 29) for now
               ifelse(flag==4|flag==6,ncol(lgtdata[[1]][[1]]$X_L)+2,               # curved 2 par original line (Kuma)
               #ifelse(flag==4|flag==6,ncol(lgtdata[[1]][[1]]$X_L)+4,              # curved 2 par new with no b3=b1b2 and b0x 
               ifelse(flag==1|flag==3,length(lgtdata[[1]][[1]]$X_L[1,L_var])+3,    # quadratic locational
               ifelse(flag==2,length(lgtdata[[1]][[1]]$X_PS[1,PS_var]),            # PS
                 1))))
    npar = ifelse(flag==3|flag==6|flag==10,nvar_L+nvar_PS,   # combined L+PS
           ifelse(flag==2,nvar_PS,                           # PS only
                          nvar_L))                           # L only
    
    CData = integer(2)
    #CData = .C("Initialize", CData=CData, as.integer(nlgt), as.integer(nvar_L), as.integer(ncol(lgtdata[[1]][[1]]$X_L)), as.integer(nvar_PS))$CData
    CData = .C("Initialize", CData=CData, as.integer(nlgt), as.integer(nvar_L), as.integer(ncol(lgtdata[[1]][[1]]$X_L[,L_var])), as.integer(nvar_PS))$CData  

     wC_all_old = list()  #  code for different # of choices for all respondents
     X_cond = list()
     nch_hh = c(rep(0,nlgt))
     totalXsize = 0
     Xsize = list()
     y_all = list()
     wC_old_table = NULL
     for(hh in 1:nlgt)                                                                                                   
     {  nch_hh[hh] = length(lgtdata[[hh]])
        XL_hh = NULL
        XP_hh = NULL
        y_hh = c(rep(0,nch_hh[hh]))
        Xsize_hh = c(rep(0,nch_hh[hh]))   # vector of the number of alternatives in each choice
        X_cond[[hh]] = matrix(0,nrow=nch_hh[hh],ncol=nvar_cond)
        for(j in 1:nch_hh[hh])
         {  XtempL=lgtdata[[hh]][[j]]$X_L[,L_var]
            XL_hh = rbind(XL_hh,XtempL)  # original code
            XP_hh = rbind(XP_hh,matrix(lgtdata[[hh]][[j]]$X_PS2[,PS_var],ncol=nvar_PS)) # removed distinction between PS and PS2, 
                                                                                        # just need to point to variables in PS, which included all PS variables
            Xsize_hh[j] = as.integer(nrow(XtempL))
            y_hh[j] = as.integer(lgtdata[[hh]][[j]]$y)
            if(is.null(filetouse[[hh]][[1]]$X_cond))
            {  X_cond[[hh]][j,] = c(rep(1,nvar_cond))
            }
            else{ X_cond[[hh]][j,] = lgtdata[[hh]][[j]]$X_cond[Cond_var]
            }
            #wC_all_old[[hh]][[j]] = as.vector(exp( matrix(oldlambdas[hh,],nrow=nlambda,byrow=T) %*% ( lgtdata[[hh]][[j]]$X_cond[1:ncov] ) ))
            #wC_all_old[[hh]][[j]]  = rep(1,nlambda)
         }
        XP_hh = matrix(as.double(XP_hh),ncol=nvar_PS) 
        .C("StoreX", CData, t(matrix(as.numeric(XL_hh),ncol=length(L_var))), t(XP_hh), as.integer(Xsize_hh), as.integer(y_hh), as.integer(hh), as.integer(nch_hh[hh]))    
        Xsize[[hh]] = Xsize_hh
        y_all[[hh]] = y_hh
        totalXsize = totalXsize+ sum(Xsize_hh)
        wC_all_old[[hh]] = X_cond[[hh]] %*%  t( matrix(oldbetas[hh,],ncol=nvar_cond, nrow=npar, byrow=T))
        wC_old_table = rbind(wC_old_table,wC_all_old[[hh]])
      }
      probvec = matrix(0,ncol=1,nrow=totalXsize)   

  #  flag values in R
  # 1 = locational quadratic
  # 2 = PS only
  # 3 = locational quadratic + PS
  # 4 = locational 2 par only (Kuma)
  # 6 = locational 2 par (Kuma) + PS
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

      ll_out_cpp = .C("LogLikelihood", CData, t(betaC), t(alphaC), ll=probvec, as.integer(flag), PredictiveFlag=PredictiveFlag)
      all_prob = ll_out_cpp$ll
          
      result = list()
      sumXsizes=0
      Hit_flag_count = 0
      Prob_sel_toMean = NULL
      for(hh in 1:nlgt)
      {
        Prob_hh = all_prob[ (1+sumXsizes):(sumXsizes+sum(Xsize[[hh]])) ]  # get probabilities for all choices and all alternatives from the long vector from C++
        sumXsizes = sumXsizes+ sum(Xsize[[hh]])   # total number of Xsizes for one respondent across all choices
        Prob_Sel = c(rep(0,nch_hh[hh]))  # length is the number of choices
        Hit_flag = c(rep(0,nch_hh[hh]))  # length is the number of choices
        Prob_ch=list()
        sumXch = 0
        for(ch_num in 1:nch_hh[hh])
        {  Prob_ch[[hh]]=list()                          # list of length nlgt of probabilities for alternatives for one choice for that respondents
           Prob_ch[[hh]][[ch_num]] = Prob_hh[ (1+sumXch):(sumXch+Xsize[[hh]][ch_num]) ]  # vector of all alternatives' probabilities for one choice for that respondent
           sumXch = sumXch + Xsize[[hh]][ch_num]
           Prob_Sel[ch_num] = Prob_ch[[hh]][[ch_num]][ y_all[[hh]][ch_num] ]        # probability of the selected alternative
           Hit_flag[ch_num] = ifelse(which.max(Prob_ch[[hh]][[ch_num]])==y_all[[hh]][ch_num] ,1,0)
        }
        Hit_flag_count = Hit_flag_count + sum(Hit_flag)
        Prob_sel_toMean = c(Prob_sel_toMean, Prob_Sel)
        result[[hh]] = list(Prob_Sel=Prob_Sel,Hit_flag=Hit_flag,Prob_hh=Prob_hh, Xsize_hh=Xsize[[hh]],y_hh=y_all[[hh]]  )
      }
    .C("Cleanup", CData)
    return(list(startBeta = oldbetas, result = result, Hit_flag_count=Hit_flag_count, MeanProb_sel = mean(Prob_sel_toMean,na.rm = TRUE),totalChoices=sum(nch_hh)))
}




