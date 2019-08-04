geno_pred <- function(phenocsv,genocsv,cvfcsv,cvf=1,mod = "BRR",label,phe)
{
    my_phe <- phe
    depends<- c("BGLR","doBy","doParallel",'R.utils',"BBmisc","dplyr")
    foo <- sapply(depends,
                  function(X){if(!suppressPackageStartupMessages(require(X,character.only = T))){install.packages(X)}})
    foo <- sapply(depends,function(X){suppressPackageStartupMessages(library(X,character.only=TRUE))})
    rm(foo)
    
        
    maze <- read.csv(genocsv, row.names = 1)
    phe <- read.csv(phenocsv, row.names = 1)
    cvffolds <- read.csv(cvfcsv,row.names=1)
    
    X <- scale(maze)
    y <- phe[[label]]
    if(any(is.na(y))){
        rms <- which(is.na(y))
        y <- y[-rms]
        X <- X[-rms,]
    }
    for(i in 1:50){
        cvf = i
        n=length(y)
        seed <- sample(1:100,1)
                                        #set.seed(seed)
                                        #folds=sample(1:cvf,size=n,replace=T)
        folds = cvffolds[,cvf]
        yHatCV=rep(NA,n)
        
        for(i in 1:max(folds)){
            cat("Predicting cv-fold ",i," of ", max(folds))
            tst=which(folds==i)
            yNA=y
            yNA[tst]=NA
            fm=BGLR(y=yNA,ETA=list(list(X=X,model=mod)),verbose =F ,nIter=7000,burnIn=1000)
            yHatCV[tst]=fm$yHat[tst]
            cat("   done\n")
        }
        
        my_cor <- cor(yHatCV,y,use = "complete.obs")
        
        print(c("Corrleation of GP", mod, my_cor))
        filename = paste0(my_phe,"_gp_results.csv")
        print(filename)
        if(!any(dir() == filename)){
            res <- matrix(ncol=8, nrow = 1) %>%
                setColNames(c("geno","pheno","cv_folds","seed","label", "cor","method","nmark"))
            res[1,] <- c(as.character(genocsv),as.character(phenocsv),as.character(cvf),as.character(seed),
                         as.character(label),as.character(my_cor),as.character(mod),dim(X)[2])
            print("#################")
	    print(res)
	    print("############")
            write.csv(res,filename)
        }else{
            res <- read.csv(filename,row.names = 1)
            for(i in 1:7){
                res[,i] <- as.character(res[,i])
            }
            res[dim(res)[1]+1,] <- c(as.character(genocsv),phenocsv,cvf,seed,as.character(label),my_cor,mod,dim(X)[2])
            write.csv(res,filename)
        }
    }
    
}


