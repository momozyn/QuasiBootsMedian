median_boots <- function(bx,bxse,by,byse,cormat,same_population, IV_list, boots_iteration = 10000, seed=314159265){
  set.seed(seed)

  library("Matrix")
  df = data.frame(IV_list, bx, bxse, by, byse)
  #df=df[order(df$IV_list),]
  # step 1 replace NA with 0 if exists
  
  cormat<-cor.smooth(ifelse(is.na(cormat), 0,cormat))


  # step 2 calculate causal estimate for each IV
  df$bxy = df$by/df$bx
  
  if(same_population == TRUE){
    df$sexy = sqrt((df$byse^2)/(df$bx^2)+(df$by^2*(df$bxse)^2)/(df$bx^4)-2*(df$by*cov(df$by/df$byse, df$bx/df$bxse))/(df$bx^3))
    
    if(any(is.na(df$sexy)) == TRUE){
      #change to 1st order for those with NA sexy
      idxs_se = which(is.na(df$sexy))
        for(idx_se in idxs_se){
          df$sexy[idx_se] = sqrt(((df$byse[idx_se])^2)/(df$bx[idx_se]^2))
        }

      
      
    }
  }else{
    df$sexy = sqrt((df$byse^2)/(df$bx^2)+(df$by^2*(df$bxse)^2)/(df$bx^4))
  }
  
  # step 3 update cormat
  if(nrow(cormat)>1){
    cormat <- cormat[order(rownames(cormat)), order(colnames(cormat))]
    cormat <- cormat[df$IV_list,df$IV_list]
  }else{
    rownames(cormat)<-df$IV_list
    colnames(cormat)<-df$IV_list
  }
  
  # Step 4 New method: bootstrap to get new SE
  C = df$sexy%*%t(df$sexy)*cormat
  median_estimate = median(df$bxy)
  
  # center bxy
  mean_bxy = mean(df$bxy)
  df$Zn = df$bxy - mean_bxy
  
  # chol for C
  
  if (class(try(chol(C))) != "try-error"){
    L = chol(C)
  }else{
    # deal with the non-PD case
    # assume C is PSD
    # V = eigen(C)$values # this produced negative values, not PSD
    
    #     
    # convert C to the nearst PD matrix
    C_nearest = nearPD(C)

    L = chol(C_nearest$mat)
    #Q = eigen(C)$vectors
    #L = Q%*%(diag(sqrt(V)))
    
    
  }
  Un = df$Zn %*% solve(L)
  
  # 
  med = NULL
  for (boot_iter in 1:boots_iteration){
    idx = sample(length(df$Zn), replace = T)
    Un_boot = Un[idx]
    Zn_boot = Un_boot %*% L
    bxy_boot = Zn_boot + mean_bxy
    
    med[boot_iter] = median(bxy_boot)
    
  }
  
  se_med = sd(med) 
  
  # assume normal
  # later: arguments assuming t-distribution 
  boot_lowerCI = median_estimate - qnorm(0.975)*se_med
  boot_upperCI = median_estimate + qnorm(0.975)*se_med
    
  pvalue = 2*pnorm(-abs(median_estimate/se_med))
  
  method = "median Bootstrap"
  results = data.frame(method, median_estimate, se_med, boot_lowerCI, boot_upperCI, pvalue)

  return(results)
}