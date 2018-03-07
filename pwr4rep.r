# b - beta, se - se, pval - pvalue, nd - N discovery, nr - N replication, nl - number of loci for replication

pwr.rep<-function(b,se,nd,nr,pval,nl) {
           chi<-qchisq(1-pval/nl,1)
           t2<-(b/se)^2
           ncp<-t2*(nr/nd)-1
           if(ncp<0) ncp<-0
           pwr<-pchisq(chi,1,ncp=ncp,lower.tail = F)
           #L<-as.list(c(round(chi,4),round(t2,4),round(ncp,4),round(pwr,4)))
           #names(L)<-c("chi2","t2","ncp","power")
           #return(L)
           return(round(pwr,4))
           }
