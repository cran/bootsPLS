# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: 13-05-2014
# last modification: 15-03-2015
# Copyright (C) 2014
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


# perform the replication of the splsda on random subsamplings of the data

bootsPLS=function(X,
        Y,
        near.zero.var = TRUE,
        many = 50,
        ncomp = 2,
        dist = c("max.dist", "centroids.dist", "mahalanobis.dist"),
        save.file,
        ratio,
        kCV = 10,
        grid,
        cpus,
        nrepeat =1,
        showProgress=TRUE)
{
    #-------------------------- bootstrapped sPLS-DA and Cross-Validation
    # X input data
    # Y factor input
    # near.zero.var: discard variables with near to zero variance, see function from mixOmics package
    # many= number of resampling. On each subsampling is performed the bootstrap sPLS-DA (with CV to choose the parameters)
    # ncomp= number of components
    # dist= which distance should be used to classify the samples?
    # save.file= name of the file with the results to save
    # ratio: proportion of samples left out in the first random subsampling
    # kCV: number of k-fold in the cross validation
    # grid: CV grid, a vector of values of keepX to test in the CV
    # nrepeat: number of replication of the Mfold process, for each of many
    # showProgress=TRUE, show the progress of the iteration
    
    if(missing(X)) stop("missing X")
    if(missing(Y)) stop("missing Y")

    check=Check.entry.bootsPLS(X,Y)
    X=check$X
    Y=check$Y

    nlevelY=nlevels(Y)
    
    if(missing(grid))
    {
        grid=1:min(40,ncol(X))
    }
    
    dist=dist[1] #if multiple entries
    
    #construct a dummy matrix
    Y.mat=matrix(0,nrow=nrow(X),ncol=nlevelY)
    for(i in 1:nlevelY)
    {
        Y.mat[which(Y==levels(Y)[i]),i]=1
    }
    colnames(Y.mat)=levels(Y)
    
    #remove variables with nearzerovar function
    nzv=list()
    if(near.zero.var == TRUE)
    {
        nzv = nearZeroVar(X)
        if (length(nzv$Position > 0))
        {
            names.remove=colnames(X)[nzv$Position]
            warning("Zero- or near-zero variance predictors have been discarded.\n See $nzv for problematic predictors.")
            X = X[, -nzv$Position,drop=FALSE]
            if(ncol(X)==0) {stop("No more predictors")}
        }else{nzv$Position=NULL}
    }
    
    #initialise some parameters that we want to record
    P=ncol(X)
    ClassifResult=array(0,c(nlevelY,nlevelY,ncomp,many))
    rownames(ClassifResult)=levels(Y)
    colnames(ClassifResult)=paste("predicted.as.",levels(Y),sep="")
    dimnames(ClassifResult)[[3]]=paste("comp.",1:ncomp,sep="")
    dimnames(ClassifResult)[[4]]=paste("iteration.",1:many,sep="")
    
    selection.variable=array(0,c(ncomp,P,many))
    dimnames(selection.variable)[[2]]=colnames(X)
    loadings.X=array(0,c(P,ncomp, many))
    dimnames(loadings.X)[[1]]=colnames(X)
    nbr.var=NULL
    
    #--------------------------------------------------------------------------------------------------------
    #---------------------		LOOP   	---------------------
    #--------------------------------------------------------------------------------------------------------
    learning.sample=matrix(0,nrow=nrow(X),ncol=many) #record which sample are in the learning set
    prediction=array(0,c(nrow(X),many,ncomp)) #record the class associated to each sample (either in learning or test set)
    rownames(learning.sample)=rownames(X)
    dimnames(prediction)[[1]]=rownames(X)
    
    for(abc in 1:many)
    {
        if(showProgress)
        cat("iteration ",abc,"\n")
        
        #--------------- 1st step: random subsampling
        
        A=suppressWarnings(random.subsampling(Y,ratio))
        A=sort(A) # to keep the same order as the data
        #A contains the sample we want to keep in the learning set, -A in the test set
        learning.sample[A,abc]=1
        
        data.learn.signature=X[A,]
        Y.learn.signature=Y[A]
        
        data.test.signature=X[-A,]
        Y.test.signature=Y[-A]
        
        #NZV
        nzv.temp=nearZeroVar(data.learn.signature)
        if(length(nzv.temp$Position)>0)
        {
            data.learn.signature=data.learn.signature[,-nzv.temp$Position]
            data.test.signature=data.test.signature[,-nzv.temp$Position]
        }
        
        #--------------- kCV on each component to find the optimal number of variables per component
        SIGN=matrix(0,nrow=ncomp,ncol=ncol(X))
        CH=NULL
        uloadings.X=NULL # the one of size p = ncol(X)
        ind.var=NULL
        

        if(TRUE) #debug: get rid of the CV part that is time-consuming
        {
            #perform a Cross validation to tune the number of variables to keep on all components
            CV = suppressWarnings(mixOmics::tune.splsda(X=data.learn.signature, Y=Y.learn.signature, ncomp=ncomp, test.keepX=grid, validation= "Mfold",
            folds = kCV, measure="BER", dist = dist, scale=TRUE, auc=FALSE, near.zero.var=near.zero.var, cpus=cpus, progressBar=showProgress,
            nrepeat=nrepeat))
            
            
            nbr.var.opt=CV$choice.keepX
            #if(nbr.var.opt==0){save(list=ls(),file=save.file)}
        }else{
            
            ind=sample(1:ncol(data.learn.scale),3)
            nbr.var.opt=length(ind)
        }
        #----------- learning the model with nbr.var.opt
        res=suppressWarnings(mixOmics::splsda(X=data.learn.signature, Y=Y.learn.signature, keepX = nbr.var.opt, ncomp=ncomp, near.zero.var=near.zero.var, scale=TRUE))
            
        #record the signature for each comp
        signature=vector("list",length=3)
        for(num.comp in 1:ncomp)
        {
            ind=which(res$loadings$X[,num.comp]!=0)
            names(ind)=colnames(data.learn.signature)[ind]
            
            signature[[num.comp]]=names(ind)

            #record the signature
            signature.value.X=matrix(0,nrow=ncol(X),ncol=1)
            a=match(names(ind),colnames(X))
            signature.value.X[a]=res$loadings$X[ind,num.comp]
            uloadings.X=cbind(uloadings.X,signature.value.X)
            SIGN[num.comp,a]=1
        }# end num.comp
        
        names(signature)=paste("comp.",1:ncomp,sep="")
        if(showProgress) {print(signature)}

        
        #----------- test of the signature on the unscaled-learning.set
        out = predict(res, newdata = data.learn.signature, dist = dist)
        predicted.learn = out$class[[1]]#only one dist, so one predicted class
        
        prediction[which(learning.sample[,abc]==1),abc,]=predicted.learn    # record the prediction
        
        #----------- test of the signature on the unscaled-test.set
        out = predict(res, newdata = data.test.signature, dist = dist)
        predicted = out$class[[1]]#only one dist, so one predicted class
        
        prediction[which(learning.sample[,abc]==0),abc,]=predicted    # record the prediction
        
        
        #save(list=ls(),file="temp.Rdata")
        #--------record of the classification accuracy for each level of Y
        for(num.comp in 1:ncomp){
            ClassifResult[,,num.comp,abc] = mixOmics::get.confusion_matrix(truth = Y.test.signature, all.levels=levels(Y), predicted = predicted[,num.comp])
        }
        selection.variable[,,abc]=SIGN
        loadings.X[,,abc]=uloadings.X
        nbr.var=rbind(nbr.var,nbr.var.opt)
    
        #calculation of the frequency of selection for each variable, on each components, after the abc replications.
        # this is done to ensured that the file saved can be re-used.
        frequency=matrix(0,nrow=ncomp,ncol=dim(loadings.X)[1])
        for(j in 1:abc)
        {
            for(k in 1:ncomp)
            {a=which(loadings.X[,k,j]!=0)
                frequency[k,a]=frequency[k,a]+1 #add 1 everytime the gene is selected
            }
        }
        frequency=frequency/abc #get the probability of selection (percentage of times each gene is selected, per component
        colnames(frequency)=colnames(X)
        
        out=list(ClassifResult=ClassifResult,loadings.X=loadings.X,selection.variable=selection.variable,frequency=frequency,
        nbr.var=nbr.var,learning.sample=learning.sample,prediction=prediction,data=list(X=X,Y=Y,dist=dist),nzv=nzv)
        structure(out,class="bootsPLS")

        data=list(X=X,Y=Y,dist=dist)
        if(!missing(save.file))
        save(ClassifResult,loadings.X,selection.variable,frequency,nbr.var,learning.sample,prediction,data,nzv,file=save.file)
        
    }
    out
    structure(out,class="bootsPLS")

}#end function

