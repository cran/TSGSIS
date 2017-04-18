#' @title Two Stage-Grouped Sure Independence Screening
#' @description The package is a beta version that provides a high-dimensional grouped variable selection approach for detection of whole-genome SNP effects and SNP-SNP interactions, as described in Fang et al. (2017, under review). The proposed TSGSIS is developed to study interactions that may not have marginal effects.
#' @param XA The \eqn{N \times P} matrix of XA. There are \eqn{N} individuals and \eqn{P} variables in matrix, with one individual in each row and one genotype in each column.
#' @param Y The \eqn{N \times 1} matrix of Y. It can be real number or binary outcome.
#' @param Gene_list The \eqn{a \times d} matrix of the Gene_list. a is the maximal number of gene size in the Gene_list which other values are denoted by 0. d is the number of genes.
#' @param ntest The ntest (\eqn{<N}) is the number of testing data for evaluation of MSE.
#' @param lambda The lambda is the parameter of Lasso regression.
#' @param Method "Reg" for quantitative trait modeling, "LR" for disease trait modeling.
#' 
#' @return Returns a result of screening
#' @return \item{result}{First element of the result is the MSE of testing data, the rest elements are the important SNP effects and SNP-SNP interactions after TSGSIS modeling.}
#' @references Yao-Hwei Fang, Jie-Huei Wang and Chao A. Hsiung (2017). TSGSIS: A High-dimensional Grouped Variable Selection Approach for Detection of Whole-genome SNP-SNP Interactions. (Under revision in Bioinformatics)
#' @note The missing value (NA) in the XA and Y is not allowed in this beta version.
#' @import glmnet MASS stats
#' @export
#' @examples
#' #We investigate the performance of TS-GSIS under model 1 with intra-gene correlation rho = 0.2, 
#' #trait dispersion sigma^2 = 1, effect size k = 3 and homogeneous MAF. 
#' #Given 100 SNPs with gene size d = 10, 500 unrelated individuals are simulated.
#' #(Please refer to the Figure 3 of the reference)
#'
#' library(glmnet)
#' library(MASS)
#' 
#' set.seed(1)# Set seed
#' #Parameter setting
#' ntotal = 500
#' p = 100
#' n.pred = 10 #Gene sizes
#' rho = 0.2 #Intra-gene correlation in block
#' k = 3 #Effect size
#' vari = 1 #Sigma2
#' lambda = 0.5 #For lasso parameter
#' ntest = 150 #For evaluation
#' Method="Reg"#For quantitative trait

#' #Heterogeneous MAF: randomly set to 0.35, 0.2 or 0.1 with equal likelihood.
#' MAF = matrix(0,2,3) 
#' MAF[,1] = c(0.1225,0.5775)
#' MAF[,2] = c(0.04,0.36)
#' MAF[,3] = c(0.01,0.19)
#' #Trait Y
#' modelY = "k*XA[,1] - k*(sqrt(rho))*XA[,5] + k*XA[,31]*XA[,5] + rnorm(ntotal,0,vari)"
#' 
#' PAS1 = function(z){ g = paste("A",z,sep = "") 
#' 	return(g)
#' }#Define colname fun.
#' norm = function(a) (a-mean(a))/sd(a) #Define standardization fun.
#'
#' #The codes of simulated data for quantitative trait are listed in the following. We use mvrnorm
#' #function to simulate the genotype data. Y is continuous with normal distribution, all errors are
#' #assumed to be normally distributed with a mean of zero and a variance of one (vari = 1).
#' out = array(0, dim=c(n.pred)) #For LOOCV
#' corrmat = diag(rep(1-rho, n.pred)) + matrix(rho, n.pred, n.pred) #Create covariance matrix with rho
#' corrmat[,5] = sqrt(rho)
#' corrmat[5,] = sqrt(rho)
#' corrmat[5,5] = 1
#' L = array(0, dim=c(n.pred, n.pred, (p/n.pred)))
#' L[,,1] = corrmat
#' for(i in 2:(p/n.pred)){ 
#' 	L[,,i] = diag(rep(1-rho, n.pred)) + matrix(rho, n.pred, n.pred)
#' }
#' temp = "bdiag(L[,,1]"
#' for (i in 2:(p/n.pred)){
#' 	temp = paste(temp,",","L[,,",i,"]", sep="")
#' }
#' temp = paste(temp,")", sep="")
#' corrmat2 = eval(parse(text=temp))
#' 
#' beta0 = matrix(0,p,1) #Simulate genotype
#' X = matrix(0,ntotal,p)
#' X = mvrnorm(ntotal, beta0, corrmat2 , tol=1e-8, empirical=FALSE)
#' XA = data.frame(X); colnames(XA) <- c(sapply(1:p,PAS1))
#' C1 = matrix(0,1,p)
#' C2 = matrix(0,1,p)
#' tempMAF = sample(3,1)
#' for (i in 1:p){
#' 	C2[1,i] = quantile(X[,i], MAF[1,tempMAF])
#' 	C1[1,i] = quantile(X[,i], MAF[2,tempMAF])
#' 	XA[X[,i] > C1[,i],i] = 1
#' 	XA[X[,i] <= C1[,i] & X[,i] >= C2[,i],i] = 0
#' 	XA[X[,i] < C2[,i],i] = -1
#' }
#' XA = apply(XA, 2, norm) #Standardization
#' 
#' Y = eval(parse(text=modelY)) #Simulate gaussian response
#' 
#' temp = 1:p
#' Gene_list = matrix(temp,nrow=n.pred) #Create Gene-based SNP set
#' #Run TSGSIS fun. with XA, Y, Gene_list, ntest (for predicted model), lambda of lasso regression,
#' #Method types: "Reg" for quantitative trait; "LR" for disease trait.
#' Screen_result = TSGSIS(XA, Y, Gene_list, ntest, lambda, Method) 

###################################################################
#              CREATE TSGSIS FUNCTION for R PACKAGE               #
###################################################################

TSGSIS <- function(XA,Y,Gene_list,ntest,lambda,Method){
  ## colnames fun.
  PAS <- function(z){ g <- paste("X",z,sep="")
  return(g)}
  PAS3 <- function(z){g <- paste("G",z,sep="")
  return(g)}
  p=dim(XA)[2]
  ntotal=dim(XA)[1]
  Block = dim(Gene_list)[2]
  n.pred <- p/Block
  Ori_XA = XA; colnames(Ori_XA) = c(sapply(1:p,PAS))#copy
  ##testdata
  n = ntotal-ntest
  XAtest = data.frame(Ori_XA)
  XAtest = XAtest[-(1:n),];rownames(XAtest)=1:ntest
  XA = XA[-((n+1):ntotal),]
  Ytest = Y[(n+1):ntotal]
  Y = Y[-((n+1):ntotal)]
  Ori_XA = Ori_XA[-((n+1):ntotal),]
  #Step2 GSIS for RSS
  RY=Y[sample(n,n)]
  #sfExport('Ori_XA','Y','XA','Gene_list','RY','Block','n.pred','n',"Method")
  ROW <- sapply(1:Block , function(i) {#Gene-based collapsing
	#stepwise regression
	from=paste(colnames(XA[,Gene_list[,i]]),collapse="+")
	from=paste("~", from , sep="")
	if (Method=="Reg"){
		model=step(glm(Y ~ 1, data=data.frame(cbind(XA[,Gene_list[,i]],Y))), direction='forward', scope=eval(parse(text=from)),steps=2,trace=0)
	} else if (Method=="LR"){
		model=step(glm(Y ~ 1, data=data.frame(cbind(XA[,Gene_list[,i]],Y))), direction='forward', scope=eval(parse(text=from)),steps=2,trace=0)
	}
	return(sum(model$residuals^2)/model$df.residual)#RSS/df
  })
  names(ROW)=c(sapply(1:Block,PAS3))
  #Permutation
  ROW2 <- sapply(1:Block , function(i) {
    from=paste(colnames(XA[,Gene_list[,i]]),collapse="+")
	from=paste("~", from , sep="")
	if (Method=="Reg"){
		model=step(glm(RY ~ 1, data=data.frame(cbind(XA[,Gene_list[,i]],RY))), direction='forward', scope=eval(parse(text=from)),steps=2,trace=0)
	} else if (Method=="LR"){
		model=step(glm(RY ~ 1, data=data.frame(cbind(XA[,Gene_list[,i]],RY))),direction='forward', scope=eval(parse(text=from)),steps=2,trace=0)
	}
	return(sum(model$residuals^2)/model$df.residual)
  })
  model1 = names(ROW)[ROW < min(ROW2)]
  #step3 pairwise interaction term & permutation data
  if (length(model1) != 0){ 
    Int_DATA<-data.frame(Ori_XA[,1])
	Int_DATA_Perm<-data.frame(Ori_XA[,1])
	RSS<-data.frame(1)
	RSS_Perm<-data.frame(1)
	dim.init1=0
	for (i in as.numeric(substr(model1,2,100000)) ) {#main effects
		for (j in setdiff(1:Block,i)){#orders
			dim.init2=0
			for (o in 1:n.pred){
				Int_DATA[(dim.init2 + 1):(dim.init2 + n.pred)] <- XA[,Gene_list[,i][o]] * XA[,Gene_list[,j]]
				Int_DATA_Perm[(dim.init2 + 1):(dim.init2 + n.pred)] <- XA[,Gene_list[,i][o]] * XA[sample(n,n),Gene_list[,j]]
				colnames(Int_DATA)[(dim.init2 + 1):(dim.init2 + n.pred)]         <- paste(colnames(XA)[Gene_list[,i][o]],":",colnames(XA)[Gene_list[Gene_list[,j]!=0,j]],sep="")
				colnames(Int_DATA_Perm)[(dim.init2 + 1):(dim.init2 + n.pred)]    <- paste(colnames(XA)[Gene_list[,i][o]],":",colnames(XA)[Gene_list[Gene_list[,j]!=0,j]],sep="")
				dim.init2=dim.init2+n.pred
			}
			from=paste(colnames(data.frame(Int_DATA)),collapse="+")#create formula
			from=paste("~", from , sep="")
			if (Method=="Reg"){
				model=step(glm(Y ~ 1, data=cbind(data.frame(Int_DATA,Y))),direction='forward', scope=eval(parse(text=from)),steps=2,trace=0)
			} else if (Method=="LR"){
				model=step(glm(Y ~ 1, data=cbind(data.frame(Int_DATA,Y))),direction='forward', scope=eval(parse(text=from)),steps=2,trace=0)
			}
			RSS[dim.init1 + 1]<-sum(model$residuals^2)/model$df.residual;colnames(RSS)[dim.init1 + 1]<- paste(names(ROW)[c(i,j)],collapse=":")
			#conditional permutation
			if (Method=="Reg"){
				model2=step(glm(Y ~ 1, data=cbind(data.frame(Int_DATA_Perm,Y))), direction='forward', scope=eval(parse(text=from)),steps=2,trace=0)
			} else if (Method=="LR"){
				model2=step(glm(Y ~ 1, data=cbind(data.frame(Int_DATA_Perm,Y))), direction='forward', scope=eval(parse(text=from)),steps=2,trace=0)
			}
			RSS_Perm[dim.init1 + 1]<-sum(model2$residuals^2)/model2$df.residual;colnames(RSS_Perm)[dim.init1 + 1]<- paste(names(ROW)[c(i,j)],collapse=":")
			dim.init1=dim.init1+1
		}
	}
	modelI=colnames(RSS)[RSS<min(RSS_Perm)]
	if (length(modelI)!=0){
		modelII=modelI#delete replicate
		replicat <- sapply(1:length(modelII) , function(i) {
			return(c(strsplit(modelII[i],split=":",fixed=T))[[1]])
		})
		replicat1=replicat
		aa=replicat[1,]
		replicat1[1,]=replicat[2,]
		replicat1[2,]=aa
		detl=""
		for (i in 1:length(modelI)){
			for (j in setdiff(i:length(modelI),i)){
				if (all(replicat[,i]%in%replicat1[,j])){
					detl=c(detl,j)
				}
			}
		}
		detl=detl[-1]
		if (length(detl)!=0){
			modelI=modelII[-(as.numeric(detl))]
		} else{modelI=modelII}
		#print(modelI)
		GGG=strsplit(modelI,split=":",fixed=T)	
		#Interaction -> SNPs		
		SIS_G <- lapply(1:length(GGG) , function(i) {
			glmmod<-substr(GGG[[i]],2,100000)#
		})
		Gene.loop <- function(o) {
			return(Ori_XA[,c(Gene_list[,as.numeric(SIS_G[[o]])])])
		}
		Da <- lapply(1:length(SIS_G),Gene.loop)
		#test data
		Gene.loop2 <- function(o) {
			return(XAtest[,c(Gene_list[,as.numeric(SIS_G[[o]])])])
		}
		Datest <- lapply(1:length(SIS_G),Gene.loop2)
		#create Step5's data
		temp1<-data.frame(Ori_XA[,Gene_list[,as.numeric(substr(model1,2,100000))[1]]])
		temp1_test<-data.frame(XAtest[,Gene_list[,as.numeric(substr(model1,2,100000))[1]]])
		from=paste(colnames(temp1), collapse="+")
		from=paste("~(", from ,")^2", sep="")
		a=model.matrix(as.formula(from), data = data.frame(temp1))
		temp1[1:(ncol(a)-1)] <-a[,2:ncol(a)]#delete first column
		colnames(temp1)[1:(ncol(a)-1)]         <- colnames(a)[2:ncol(a)]
		#test
		a=model.matrix(as.formula(from), data = data.frame(temp1_test))
		temp1_test[1:(ncol(a)-1)] <-a[,2:ncol(a)]#delete first column
		colnames(temp1_test)[1:(ncol(a)-1)]         <- colnames(a)[2:ncol(a)]
		dim.init1=ncol(temp1)
		if (length(model1)>1){
			for (i in 2:length(model1)){
				aa<-data.frame(Ori_XA[,Gene_list[,as.numeric(substr(model1,2,100000))[i]]])
				from=paste(colnames(aa), collapse="+")
				from=paste("~(", from ,")^2", sep="")
				a=model.matrix(as.formula(from), data = data.frame(aa))
				temp1[(ncol(temp1)+1):(ncol(temp1)+ncol(a)-1)] <-a[,2:ncol(a)]#delete first column
				colnames(temp1)[(dim.init1+1):(dim.init1+ncol(a)-1)]         <- colnames(a)[2:ncol(a)]
				#test
				aa<-data.frame(XAtest[,Gene_list[,as.numeric(substr(model1,2,100000))[i]]])
				from=paste(colnames(aa), collapse="+")
				from=paste("~(", from ,")^2", sep="")
				a=model.matrix(as.formula(from), data = data.frame(aa))
				temp1_test[(ncol(temp1_test)+1):(ncol(temp1_test)+ncol(a)-1)] <-a[,2:ncol(a)]#
				colnames(temp1_test)[(dim.init1+1):(dim.init1+ncol(a)-1)]         <- colnames(a)[2:ncol(a)]			
				dim.init1=ncol(temp1)
			}
		}
		dim.init2=ncol(temp1)
		for (i in 1:length(SIS_G) ) {
			for (j in setdiff(Gene_list[,as.numeric(SIS_G[[i]])[1]],0)){
				temp1[(dim.init2 + 1):(dim.init2 + n.pred)] <- Ori_XA[,j] * Ori_XA[,Gene_list[,as.numeric(SIS_G[[i]])[2]]]
				colnames(temp1)[(dim.init2 + 1):(dim.init2 + n.pred)] <- paste(colnames(Ori_XA)[j],":",colnames(Ori_XA)[Gene_list[,as.numeric(SIS_G[[i]])[2]]],sep="")
				#test                                                         
				temp1_test[(dim.init2 + 1):(dim.init2 + n.pred)] <- XAtest[,j] * XAtest[,Gene_list[,as.numeric(SIS_G[[i]])[2]]]
				colnames(temp1_test)[(dim.init2 + 1):(dim.init2 + n.pred)]         <- paste(colnames(XAtest)[j],":",colnames(XAtest)[Gene_list[,as.numeric(SIS_G[[i]])[2]]],sep="")
				dim.init2=dim.init2+n.pred
			}
		}	
	} else {#no interaction
		SIS_G <- lapply(1:length(model1) , function(i) {
			glmmod<-substr(model1[i],2,100000)#create main effects
		})
		Gene.loop <- function(o) {
			return(Ori_XA[,c(Gene_list[,as.numeric(SIS_G[[o]])])])
		}
		Da <- lapply(1:length(SIS_G),Gene.loop)
		temp1=do.call(cbind,Da)
		Gene.loop2 <- function(o) {
			return(XAtest[,c(Gene_list[,as.numeric(SIS_G[[o]])])])#test data
		}
		Datest <- lapply(1:length(SIS_G),Gene.loop2)
		temp1_test=do.call(cbind,Datest)
	}
    #Step5:Lasso
    MM = "(Intercept)"
	if (Method=="Reg"){
		glmmod<-glmnet(data.matrix(temp1),Y,alpha=1,lambda=lambda[1])
	}else if (Method=="LR"){
		glmmod=glmnet(data.matrix(temp1),Y,alpha=1,,lambda=lambda[1],family='binomial')
	}
	MM=union(MM,rownames(coef(glmmod))[grep("TRUE",coef(glmmod)!=0)])
  } else { #no model1
    MM = "no genetic effect"
  }
  #prediction
  if (Method=="Reg"){
	if (length(MM)!=1){
		MM[1] = paste("MSE:",as.character(mean((Ytest - coef(glmmod)[1] -   as.matrix(temp1_test)   %*%  coef(glmmod)[2:length(coef(glmmod))]   )^2)), sep="")
	} else{
		MM[1] = paste("MSE:",as.character(mean((Ytest - mean(Ytest))^2)), sep="")
	}
  } else if (Method=="LR"){
	if (length(MM)!=1){
		predicted_y <- predict(glmmod,as.matrix(temp1_test), type='response')
		predicted_y[predicted_y>0.5]  <- 1
		predicted_y[predicted_y<=0.5] <- 0
		confusion_matrix<- ftable(Ytest, predicted_y)
		MM[1] = paste("MSE:",as.character(sum(diag(confusion_matrix)) / length(predicted_y)), sep="")
	} else{
		MM[1] = paste("MSE:",as.character(length(grep("TRUE",Ytest==rep(0.5*(sign(mean(Ytest)-0.5)+1),ntest)))/ntest), sep="")
	}
  }
  return(MM)
}
