# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: 15-12-2017
# last modification: 15-12-2017
# Copyright (C) 2017
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




plot.predictCI=function(x,ncomp=1, level=1, las=2, col, title, name.var = TRUE, abline=TRUE, ...)
{
    # x: a predictCI object, as from CI.prediction() or $out.CI from prediction()
    
    # we extract the confidence interval for the relevant component and level of Y
    if(is.numeric(level))
    {
        CI.X = x$CI[[ncomp]][[level]]
        ind=level
    }else if(is.character(level)){
        ind = which(level == names(x$CI[[ncomp]]))
        if(length(ind)==0)
        stop("`level' does not seem to be one of your outcome categories")
        
        CI.X = x$CI[[ncomp]][[paste(level)]]
    }
    if(ncol(CI.X)!=2) stop("Problem with the entry. CI[[ncomp]][[level]] does not have two columns")
    
    if(missing(col)) {
        col=rep(color.mixo(1),nrow(CI.X))
    } else {
        if(all(length(col)!=c(1,nrow(CI.X))))
        stop("'col' needs to be a vector of length 1 or ", nrow(CI.X))
    }
    if(missing(name.var)) name.var = rownames(CI.X)
    
    if(!is.logical(name.var))
    {
        if(length(name.var)!= nrow(CI.X))
        stop("'name.var' should be a vector of length ", nrow(CI.X))
        axis.label = TRUE
    } else if (name.var){
        name.var = substr(rownames(CI.X),1,20)
        axis.label = TRUE
    }
    if(missing(title))
    title = paste0("Confidence Intervals on comp ",ncomp," based on the ", names(x$CI[[ncomp]])[ind] ," level")


    plot(1:nrow(CI.X),1:nrow(CI.X),las=las,ylab="",ylim=c(min(0,min(CI.X)),max(1,max(CI.X))),
    xaxt="n",xlab="",main=title,type="n", ...)
    
    axis(at=1:nrow(CI.X),side=1,las=las,labels=name.var)

    for(i in 1:nrow(CI.X))
    {
        arrows(i,CI.X[i,1],i,CI.X[i,2],angle=90,code=3,lwd=2,length=0.1, col=col[i])
    }
    
    if(abline==TRUE)
    {
        abline(h=0.5, lty=2)
        abline(h=1/length(x$CI[[ncomp]]), lty=2)
    }
}

