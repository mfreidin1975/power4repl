### generic
pwr.repl<-function(nl,    # number of loci to replicate
                   pvalD, # critical p-value for discovery
                   b,     # effect in discovery
                   Nd,    # N in discovery
                   Nr)    # N in replication
{
pval<-0.05/nl # p-value for replication, depends on the number of loci to replicate
chi<-qchisq(1-pval,1) # corresponding threshold chi2; same as qchisq(pval,1,low=F)
se<- abs(b/qnorm(pvalD/2)) # se in discovery, significant for given pvalue
T2<-(b/se)^2
ncp<-T2*(Nr/Nd)-1
return(pchisq(chi,1,ncp=ncp,lower.tail = F)) # power of replication
}

pwr.repl(nl=10,pvalD=0.05/200,b=0.05,Nd=1200,Nr=1300)

### need to find chi from p-value

nl<-35 # number of loci to replicate
pval<-0.05/nl # p-value for replication, depends on the number of loci to replicate
chi<-qchisq(1-pval,1) # corresponding chi2; same as qchisq(pval,1,low=F)


b<- -0.0419	 # effect in discovery
se<- 	0.00935 # se in discovery

Nd<-119823 # N in discovery
Nr<-100000 # N in replication

T2<-(b/se)^2
ncp<-T2*(Nr/Nd)-1

pchisq(chi,1,ncp=ncp,lower.tail = F) # power of replication

# b - beta, se - se, pval - pvalue, nd - N discovery, nr - N replication, nl - number of loci for replication

pval<-0.05
nl<-129

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

###### for Helena
x<-read.table("C:\\Users\\k1471250\\Desktop\\SigHAid_SNPs_HRC.txt",h=T)
x<-x[order(x$P),]

pwr.rep<-function(b,se,nd,nr,pdisc,alpha=0.05,nl=1) {
           chi<-qchisq(1-alpha/nl,1)
           t2<-(b/se)^2
           #ncp<-t2*(nr/nd)-1 
           # why is -1 here? 
           # in Yakov's paper it is (t2-1)*(nr/nd)
           # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6353874/
           ncp<-(t2-1)*(nr/nd) # this is not really correct, so ignore!!!
           if(ncp<0) ncp<-0
           pwr<-pchisq(chi,1,ncp=ncp,lower.tail = F)
           # alternative; https://github.com/kaustubhad/gwas-power
           # also example https://www.biostat.washington.edu/sites/default/files/modules/2017_SISG_14_5.pdf
           q2<-t2/nd
           ncp.alt<-nr*q2/(1-q2)
           pwr.alt<-pchisq(chi,1,ncp=ncp.alt,lower.tail = F)
           # another alternative based on p-value
           if(missing(pdisc)) pdisc<-2*pnorm(-abs(b/se))
           t2<-qchisq(1-pdisc,1)
           eff<-t2/nd
           ncp.rev<-nr*eff
           pwr.rev<-pchisq(chi,1,ncp=ncp.rev,lower.tail = F)
           return(data.frame(by.Yakov=pwr,by.Visscher=pwr.alt,by.Rev=pwr.rev))
           }

## Explanation about -1 (Sodbo)
## https://trello.com/c/x0zF1dtl
## NCP is the difference between mean of central distribution under null and under alternative
## Mean of central chi2 == 1, so we subtract it.


res<-apply(x[,c("BETA","SE","OBS_CT")],1,function(y)pwr.rep(y["BETA"],y["SE"],y["OBS_CT"],nr=100e3,pval=0.05,nl=128))

##### Another version - not correct as it give 2 times smaller power estimates???

pwr.rep.rev<-function(pdisc,nd,nr) { # pdisc - p-value in discovery
     t2<-qchisq(1-pdisc,1)
     chi<-qchisq(1-0.05,1) # chi for p=0.05
     eff<-t2/nd
     ncp<-nr*eff
     pwr<-pchisq(chi,1,ncp=ncp,lower.tail = F)
     return(pwr)
     }
   
     
#of ~5e-8 for rs10927035 with Hdiff, corresponding to a chi2 of ~30.
#The sample size at discovery is 250K,
#which leads to an effect size beta_hat^2  =30/250K = 0.00012.
#At replication, the sample size is 30K, so that the ncp at replication =N*beta^2 ~ 3.5.
#Drawing a chi2 distribution with an ncp at 3.5, I have approximately a 45% power to reach the 0.05 nominal significance threshold.
#However Table S1 reports a replication power of 0.7956.

#### Final function used in hearing paper
pwr.rep<-function(b,se,nd,nr, alpha=0.05,nl=1) { # nl=1 is for nominal alpha; change nl to the number of loci to replicate, e.g. 41
           chi<-qchisq(1-alpha/nl,1)
           t2<-(b/se)^2
           q2<-t2/nd
           ncp<-nr*q2/(1-q2)
           pwr<-pchisq(chi,1,ncp=ncp,lower.tail = F)
           return(pwr)
           }

###### Based on https://github.com/kaustubhad/gwas-power/blob/master/power_calc_functions.R
###### When using maf, it is always quite low
power_by_hsq <- function(n, b, se, maf, alpha=0.05, nl = 100) {

pval<-alpha/nl
th=qchisq(pval,df=1,lower.tail=F) # threshold chi2 for number of SNPs

# based on n and qsq
t2<-(b/se)^2
qsq<-t2/nd
ncp=n*qsq/(1-qsq)
pow.N.Q2=pchisq(th,df=1,lower.tail=F,ncp=ncp)

# based on b and maf
q2=2*maf*(1-maf)*(b^2)
ncp=n*q2/(1-q2)
pow.B.MAF=pchisq(th,df=1,lower.tail=F,ncp=ncp)

return(data.frame(pow.N.Q2,pow.B.MAF))
}


#### From GCTA
## https://cnsgenomics.shinyapps.io/gctaPower/
alpha=0.05
chi2<-qchisq(alpha,df=1,lower.tail=F)
pwr=1-pchisq(chi2,h2^2/se^2,df=1)