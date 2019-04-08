########################################################
########################################################
########################################################

noninf_ANOVA <- function(y_vec, x_vec, epsilon=NA, delta=NA, alpha=0.05){

complete_data<-na.omit(data.frame(y_vec, x_vec))
y_vec <-complete_data$y_vec
x_vec <-complete_data$x_vec

    group <- factor(x_vec)
    size <- table(x_vec)
    ng <- length(size)
 
 if(is.na(epsilon)){ epsilon <- sqrt(length(size)*delta/(1-delta))}
 if(is.na(delta)){ delta <- (epsilon^2)/(length(size)+ (epsilon^2))}
  
    owt <- oneway.test(y_vec ~ x_vec, var.equal = TRUE)
    f <- owt$statistic
    N <- length(y_vec)   
	ybar_vec<-apply(cbind(1:ng) , 1, function(j){ mean(y_vec[x_vec==unique(x_vec)[j]])})
    SSw<-sum(unlist(apply(cbind(1:ng), 1,function(j) {(y_vec[x_vec==unique(x_vec)[j]] - ybar_vec[j])^2})))


#   Point estimate for cribbie2:    
 	cribbie2_hat <-   ((size/mean(size))%*%((ybar_vec-(size%*%ybar_vec/(sum(size))))^2)) /   
						( ((length(y_vec)-ng)^(-1))*(SSw))

#   Point estimate for eta2, using notation as in Okada (2013): 
	SSb <- sum((size%*%((ybar_vec-(size%*%ybar_vec/(sum(size))))^2)))
	SSt <- SSb + SSw
	eta2_hat <-  SSb/SSt

    
    pval_cribbie <- pf(f, df1 = ng - 1, df2 = owt$parameter[2], 
        ncp = mean(size) * (epsilon^2) )
	names(pval_cribbie)<-NULL
	
	
	pval_eta2 <-( pf(	f, ng - 1, df2 = owt$parameter[2],
					ncp = (delta * N) / (1 - delta)) )
	names(pval_eta2)<-NULL


	return(list(epsilon= epsilon , delta = delta , pval_cribbie = pval_cribbie,  pval_eta2= pval_eta2, cribbie2_hat= cribbie2_hat, cribbie_hat= sqrt(cribbie2_hat), eta2_hat= eta2_hat))
}

########################################################

noninf_wellekANOVA <- function(y_vec, x_vec, epsilon=NA, delta=NA, alpha=0.05){

complete_data<-na.omit(data.frame(y_vec, x_vec))
y_vec <-complete_data$y_vec
x_vec <-complete_data$x_vec

    group <- factor(x_vec)
    size <- table(x_vec)
    ng <- length(size)
    N <- length(y_vec) 
 
 if(is.na(epsilon)){ epsilon <- sqrt(length(size)*delta/(1-delta))}
 if(is.na(delta)){ delta <- (epsilon^2)/(length(size)+ (epsilon^2))}
  
    ybar_vec<-apply(cbind(1:ng) , 1, function(j){ mean(y_vec[x_vec==unique(x_vec)[j]])})    
	s2_vec<-apply(cbind(1:ng) , 1, function(j){ var(y_vec[x_vec==unique(x_vec)[j]])})        
    w_vec <- size/s2_vec
    ybarprime<-sum(w_vec%*%ybar_vec)/sum(w_vec)  
    
    Fprime_top <- w_vec%*%((ybar_vec-ybarprime)^2)/(ng-1)
    Fprime_bottom <- 1 + (2*(ng-2)/(ng^2-1)) * sum((1/(size-1))%*%(1- w_vec/sum(w_vec) )^2)   
    Fprime<-Fprime_top/Fprime_bottom
    
    dfprimetop <-		ng^2-1 
    dfprimebottom <- 	3*sum((1/(size-1))%*%(1- w_vec/sum(w_vec) )^2)
    dfprime <- dfprimetop/dfprimebottom
      
	ybar_vec<-apply(cbind(1:ng) , 1, function(j){ mean(y_vec[x_vec==unique(x_vec)[j]])})
    SSw<-sum(unlist(apply(cbind(1:ng), 1,function(j) {(y_vec[x_vec==unique(x_vec)[j]] - ybar_vec[j])^2})))
    
	cribbie2_hat <-   Fprime*((ng-1)/mean(size)) 
    eta2_hat <- cribbie2_hat/(length(size)+ cribbie2_hat)

#### p-value calculation:        
#    owt <- oneway.test(y_vec ~ x_vec, var.equal = FALSE)        
#    Fprime<-owt$statistic
#    dfprime<-  owt$parameter[2]  

    pval_cribbie <-( pf(	Fprime, ng - 1, df2 = dfprime,
					ncp = mean(size) * (epsilon^2) ))    
    names(pval_cribbie)<-NULL


    pval_eta2 <-( pf(	Fprime, ng - 1, df2 = dfprime,
					ncp = (delta * N) / (1 - delta) ))    
    names(pval_cribbie)<-NULL

	return(list(epsilon= epsilon , delta = delta , pval_cribbie = pval_cribbie,  pval_eta2= pval_eta2, cribbie2_hat= cribbie2_hat, cribbie_hat= sqrt(cribbie2_hat), eta2_hat= eta2_hat))
}


########################################################
########################################################