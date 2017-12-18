# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: 28-05-2014
# last modification: 10-10-2014
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


# random.subsampling(Y)
# prediction.formula(data,ncomp,uloadings,vloadings,CH,means.X,means.Y,sigma.X,sigma.Y)
# pre.screening(X,coeff)

# random subsamplings
fit.model=function(object,auto.tune,X,Y,ncomp,signature,alpha,limit,showProgress=TRUE)
{
    # object: class bootsPLS
    # auto.tune: tune ncomp and signature. only works with object
    
    # X: data matrix
    # Y: factor of observation
    # ncomp: number of component chosen
    # signature: lisst of which variables to keep on each of the ncomp components
    # alpha: level of the test for tuning of the ncomp and/or signature
    # limit: maximal number of genes to include in the model when tuning the signature parameter. It's a vector of length ncomp
    

    if(missing(object))
    {
        if(missing(X)) stop("missing X")
        if(missing(Y)) stop("missing Y")
        
        #check the entries
        check=Check.entry.bootsPLS(X,Y)
        X=check$X
        Y=check$Y
        
        if(missing(ncomp) | missing(signature)) stop("no tuning possible for ncomp and/or signature without a bootsPLS object")

    }else{
        if(class(object)!="bootsPLS") stop("problem")
        X=object$data$X
        Y=object$data$Y #factor
        
    }

    out.temp=list()

    if(missing(auto.tune))
    {
        auto.tune=FALSE
        if(missing(ncomp)&missing(signature))
        stop("Either auto.tune=TRUE or ncomp & signature need to be arguments")
    }
    
    if(auto.tune==TRUE)
    {
        if(missing(object)) stop("no tuning possible for ncomp and/or signature without a bootsPLS object")
        if(!missing(ncomp)|!missing(signature)) message("Auto.tune has been set to TRUE, arguments ncomp and/or signature are ignored")
        
        if(showProgress)
        message("tuning number of component")
        #optimisation of the number of components
        comp.select=component.selection(object,alpha, showProgress=showProgress)
        ncomp=comp.select$opt
        out.temp$component.selection=comp.select
        
        if(showProgress)
        message("tuning number of variables on each component")
        #optimisation of the number of variables on each components
        var.select=variable.selection(object,ncomp,alpha=alpha,limit=limit,showProgress=showProgress)
        signature=var.select$signature
        out.temp$variable.selection=var.select


    }else{
        if(missing(ncomp))
        {
            if(showProgress)
            message("tuning number of component")
            #optimisation of the number of components
            comp.select=component.selection(object,alpha,showProgress=showProgress)
            ncomp=comp.select$opt
            out.temp$component.selection=comp.select

        }
        if(missing(signature))
        {
            if(showProgress)
            message("tuning number of variables on each component")
            #optimisation of the number of variables on each components
            var.select=variable.selection(object,ncomp,alpha=alpha,limit=limit,showProgress=showProgress)
            signature=var.select$signature
            out.temp$variable.selection=var.select
        }
        
        
    }
    
    #check signature and ncomp
    if(ncomp>length(signature)) stop("The number of components has to be lower than the length of signature")
    if(ncomp!=length(signature)) {signature.temp=signature; for(i in (ncomp+1):length(signature)) signature.temp[[i]]=NULL ; signature=signature.temp}
    
    if(ncomp>length(signature)) stop("The number of components has to be lower than the length of signature")
    if(ncomp!=length(signature)) {signature.temp=signature; for(i in (ncomp+1):length(signature)) signature.temp[[i]]=NULL ; signature=signature.temp}
 
 
    
    data=X[,unique(unlist(signature)),drop=FALSE] #pick the genes
    if(ncol(data)<1) stop("problem  with signature, should be at least one gene")
    
    nlevelY=nlevels(Y)
    #construct a dummy matrix
    Y.mat=matrix(0,nrow=nrow(X),ncol=nlevelY)
    for(i in 1:nlevelY)
    {
        Y.mat[which(Y==levels(Y)[i]),i]=1
    }
    
    colnames(Y.mat)=levels(Y)
    
    
    # model learnt on the scale data and Y, the prediction will be done on the unscaled matrices using the means and variances
    data.init=data
    data.scale=scale(data)
    Y.mat.scale=scale(Y.mat)
    means.X=attr(data.scale,"scaled:center")
    sigma.X=attr(data.scale,"scaled:scale")
    means.Y=attr(Y.mat.scale,"scaled:center")
    sigma.Y=attr(Y.mat.scale,"scaled:scale")
    
    # put names in signature
    signature=lapply(signature,function(x){temp=match(x,colnames(data));x=colnames(data)[temp]})

    #remove the variables with no variance to update signature
    remove=which(sigma.X==0)
    if(length(remove)>0)
    {
        names.remove=colnames(data)[remove]
    }else{names.remove=NULL}
    #match signature and data after removing `remove'
    signature=match.signature(data,names.remove=names.remove,signature)
    
    #construct a temporary signature, suited for `data'
    signature.temp= lapply(signature,function(x){match(x,colnames(data))})


    fit=spls.hybrid(X=data.scale,Y=Y.mat.scale,ncomp=ncomp,keepX.constraint=signature.temp,near.zero.var=FALSE) # removed the variable with null variance

    #calcul explained variance
    explX=explained_variance(fit$X,fit$variates$X,ncomp)
    explY=explained_variance(fit$Y,fit$variates$Y,ncomp)

    out=fit
    out$ind.mat=out$Y
    out$Y=Y # replace Y by a factor to match splsda
    out$data=list(X=X,Y=Y,Y.mat=Y.mat,signature=signature)
    out$coeff=list(means.X=means.X,sigma.X=sigma.X,means.Y=means.Y,sigma.Y=sigma.Y)
    out$component.selection=out.temp$component.selection
    out$variable.selection=out.temp$variable.selection
    out$explained_variance=list(X=explX,Y=explY)
        
    structure(out,class= c("spls.constraint","splsda","spls","DA"))
    
}
