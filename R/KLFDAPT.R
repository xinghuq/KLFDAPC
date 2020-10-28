




KLFDAPT=function(klfdapc, gdsfile, estimation="auto",snp.id = NULL, n.RD = NULL, num.thread = 1L,
                  with.id = TRUE, outgds = NULL, verbose = TRUE)
{
 rv=SNPRelate::snpgdsPCACorr(pcaobj=klfdapc, gdsobj=gdsfile, snp.id=snp.id, eig.which=n.RD, num.thread=num.thread,
                           with.id=with.id, outgds=outgds, verbose=verbose)



    DLqvalues=function(DL_data,K,estimation)
    {
      require(robust)
      require(qvalue)
      loadings<-DL_data# [,1:as.numeric(K)]
      resscale <- apply(loadings, 2, scale)
      resmaha <- robust::covRob(resscale, distance = TRUE, na.action= na.omit, estim=estimation)$dist
      lambda <- median(resmaha)/qchisq(0.5,df=K)
      reschi2test <- stats::pchisq(resmaha/lambda,K,lower.tail=FALSE)
      qval <- qvalue::qvalue(reschi2test)
      q.values_DL<-qval$qvalues
      padj <- stats::p.adjust(reschi2test,method="bonferroni")
      return(data.frame(p.values=reschi2test, q.values=q.values_DL,padj=padj))
    }

    codata=as.data.frame(t(rv$snpcorr))
    rownames(codata)=rv$snp.id

    ### this step potentially remove the missing values
    q.values=DLqvalues(codata,K=length(n.RD),estimation=estimation)
  out=list(cor=rv,q.values=q.values)
  class(out)="KLFDAPT"
    return(out)

  }
