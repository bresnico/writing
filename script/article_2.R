# Partie NBR----

# Phil est chou

# Chargement des libraries----

library(tidyverse) # manipulation selon les règles tidy
library(rio) # pour la conversion sav <-> csv ou/et xlsx

# Acquisition des 6 bdd et conversion----
# en xlsx pour gagner du temps dans la manipulation
# Ne run que si les fichiers ne sont pas visibles

if (!file.exists("data/bdd1.xlsx")) {
  for (i in 1:6) {
  a <- paste("data/bdd",i,".sav", sep = "") # fichiers reçus par J. Forest
  b <- paste("data/bdd",i,".xlsx", sep = "") # fichiers préparés pour être compilés à la main.
  convert(a,b)
  }
}

# Conversion du Master xlsx -> csv
if (!file.exists("data/quebec2020_jf_master.csv")) {
  convert("data/bdd_master.xlsx", "data/quebec2020_jf_master.csv")
}

# Partie PGA regroupée par NBR----

# Libraries supplémentaires----

# Sont-elles toutes nécessaires ?
library(mice)
library(psy)
library(psychometric)
library(sem)
library(moments)

# Préparation des fonctions "maison"----

# Function sumScoreDf----
# calculer, à partir des données brutes, les scores composites (moyennes ou sommes des scores par dimension) = varVD.df

"sumScoreDf.fct" <-
  function(dataConv.df, param.df) {
    # Mathieu d'Acremont
    # version 1.3
    #
    # first convert data with convScore.fct
    # to handle NA, use replaceNaDf.fct
    # finally use sumScoreDf.fct
    #
    #
    nVar <- length(dataConv.df)
    nCoef <- nrow(param.df)
    if (nVar != nCoef) {
      cat("nbr of coef in param is not equal to nbr of variables\n")
    }
    else {
      nParam <- length(param.df)
      NamParam <- names(param.df)
      nSubj <- nrow(dataConv.df)
      NamSubj <- rownames(dataConv.df)
      res.aa <- NULL
      
      calcScore.fct <- function(data.v, param.v) {
        indexParam.bi <- param.v == -1 | param.v == +1
        sum(data.v[indexParam.bi])
      }
      
      for (iParam in 1:nParam) {
        param.v <- param.df[,iParam]
        res.v <- apply(dataConv.df, 1, calcScore.fct, param.v)
        res.aa <- c(res.aa, res.v)
      }
      
      dim(res.aa) <- c(nSubj, nParam)
      rownames(res.aa) <-  NamSubj
      colnames(res.aa) <- NamParam
      as.data.frame(res.aa)
    }
  }

# Function convScoreDf----
# calculer, à partir des données brutes, les scores convertis (items à inverser) = varQConv.df


"convScoreDf.fct" <-
  function(data.df, param.v, scoreMin, scoreMax) {
    # Mathieu d'Acremont
    # version 1.2
    inv.bi <- param.v == -1
    if (sum(inv.bi) > 0) {
      scoreMax + scoreMin - data.df[,inv.bi] -> data.df[,inv.bi]
    }
    data.df
  }

# Function col2rowRF----
# pas obligatoire (doit organiser les données pour corrélations bivariées)...

"col2rowRF.fct" <- function(data.df) {
  # Mathieu d'Acremont
  # version 1.1
  dimDf <- dim(data.df)
  data.aa <- apply(data.df, 2, c)
  dim(data.aa) <- c(dimDf[1] * dimDf[2], 1)
  
  row.names(data.df) -> rowNam.v
  Subj.fo <- factor(rep(rowNam.v, dimDf[2]))
  Block.v <- as.numeric(gl(dimDf[2], dimDf[1]))
  
  data.frame(Subj=Subj.fo, Block=Block.v, X=data.aa)
}

# Function corTestDf----
# permet d'obntenir les corrélations avec les IC pour les variables choisies d'un DF

corTestDf.fct <- function(first.tmp, second.tmp, method="pearson", subset=NULL) {
  # Mathieu d'Acremont, Mathieu.Dacremont@pse.unige.ch
  # University of Geneva
  # version 1.0
  
  first.df <- as.data.frame(first.tmp)
  second.df <- as.data.frame(second.tmp)
  
  if (nrow(first.df) != nrow(second.df)) {return("Data frames have different # of rows. Please check")} # exit
  
  varNam1 <- names(first.df)
  varNam2 <- names(second.df)
  
  if (length(subset) != 0){ # a subset is given
    if(length(subset) != nrow(first.df)) {return("Subset do not equal the # of rows in the data frames. Please check")} # exit
    first.df <- as.data.frame(first.df[subset,])
    second.df <- as.data.frame(second.df[subset,])
  }
  
  
  l1 <- length(first.df)
  l2 <- length(second.df)
  
  index1 <- 1:l1
  index2 <- 1:l2
  
  cor.ma <- matrix(999, l1, l2)
  corLower.ma <- matrix(999, l1, l2)
  corUpper.ma <- matrix(999, l1, l2)
  df.ma <- matrix(999, l1, l2)
  t.ma <- matrix(999, l1, l2)
  p.ma <- matrix(999, l1, l2)
  
  for (i1 in index1) {
    for (i2 in index2) {
      v1 <- first.df[,i1]
      v2 <- second.df[,i2]
      
      nbrEqual <- sum(v1 == v2, na.rm=T)
      lv1 <- sum(v1 == v1, na.rm=T)
      lv2 <- sum(v2 == v2, na.rm=T)
      if (nbrEqual == max(lv1, lv2)) {
        cor.ma[i1,i2] <- NA
        corLower.ma[i1,i2] <- NA
        corUpper.ma[i1,i2] <- NA
        t.ma[i1,i2] <- NA
        df.ma[i1,i2] <- NA
        p.ma[i1,i2] <- NA
      }
      else {
        if (is.factor(v1)) {v1 <- as.numeric(v1)}
        if (is.factor(v2)) {v2 <- as.numeric(v2)}
        test <- cor.test(v1, v2, alternative="two.sided", method=method)
        if (method == "pearson") {
          cor.ma[i1,i2] <- test$estimate
          corLower.ma[i1,i2] <- test$conf.int[1]
          corUpper.ma[i1,i2] <- test$conf.int[2]
          df.ma[i1,i2] <- test$parameter
          t.ma[i1,i2] <- test$statistic
          p.ma[i1,i2] <- test$p.value
        }
        else {
          cor.ma[i1,i2] <- test$estimate
          corLower.ma[i1,i2] <- NA
          corUpper.ma[i1,i2] <- NA
          df.ma[i1,i2] <- NA
          t.ma[i1,i2] <- NA
          p.ma[i1,i2] <- test$p.value
        }
        
      }
    }
  }
  rownames(cor.ma) <- varNam1
  colnames(cor.ma) <- varNam2
  
  rownames(corLower.ma) <- varNam1
  colnames(corLower.ma) <- varNam2
  
  rownames(corUpper.ma) <- varNam1
  colnames(corUpper.ma) <- varNam2
  
  rownames(t.ma) <- varNam1
  colnames(t.ma) <- varNam2
  
  rownames(df.ma) <- varNam1
  colnames(df.ma) <- varNam2
  
  rownames(p.ma) <- varNam1
  colnames(p.ma) <- varNam2
  
  list(r=cor.ma, lower=corLower.ma, upper=corUpper.ma, t=t.ma, df=df.ma, p=p.ma)
}

# Function formCor----
# mise en forme selon la fonction précédente corDf

formCor.fct <- function(cor.li, digits=2) {
  # format correlations returned by cor.test.fct
  r.v <- formatC(cor.li$r, digits, format="f")
  lower.v <- formatC(cor.li$lower, digits, format="f")
  upper.v <- formatC(cor.li$upper, digits, format="f")
  
  cor.li$p -> p.v
  rep("", length(p.v)) -> sig.v
  sig.v[p.v < .05] <- "*"
  
  res.v <- paste(r.v, sig.v , " (", lower.v, " ", upper.v, ")", sep="")
  #  replace(res.v, res.v == " NA ( NA  NA)", " NA") -> res.v
  gsub("0\\.", ".", res.v) -> res.v
  dim(res.v) <- dim(r.v)
  colnames(res.v) <- colnames(r.v)
  rownames(res.v) <- rownames(r.v)
  cat("\nCorrelations within CI")
  cat("\n\n")
  print(res.v, quote=F)
}


# Function ICC.C1----
# calcule les IC pour les alpha (pas obligatoire)

ICC.C1.fct <- function(data.df, alpha) {
  # calculate the intra-class corelation, consistency version
  # data.df has the subj in row and the repeated measure in col
  # subj represent a random gouping factor
  
  data.ma <- as.matrix(data.df)
  
  n <- nrow(data.ma) # nbr of subj
  k <- ncol(data.ma) # nbr of repeated measures
  DFr <- n - 1
  DFe <- (n - 1)*(k - 1)
  
  Y.. <- mean(data.ma)
  Yi. <- apply(data.ma, 1, mean)
  Y.j <- apply(data.ma, 2, mean)
  
  # similar as aov(X ~ Subj + as.factor(Block), data=dataRF.df)
  SSr <- 0 # Sum of Square, Between row
  SSe <- 0 # Error
  for (i in 1:n) { # go through row
    Square <- (Yi.[i] - Y..)^2
    SSr <- SSr + Square
    for (j in 1:k) { # go through col
      Square <- (data.ma[i,j] - Yi.[i] - Y.j[j] + Y..)^2
      SSe <- SSe + Square
    }
  }
  SSr <- k*SSr
  
  MSr <- SSr / DFr # Mean Square, Between row
  MSe <- SSe / DFe # Error
  Fobs <- MSr / MSe
  
  FLtab <- qf(1-alpha/2, DFr, DFe)
  FUtab <- qf(1-alpha/2, DFe, DFr)
  
  FL <- Fobs / FLtab
  FU <- Fobs * FUtab
  
  ICC <- (MSr - MSe) / (MSr + (k - 1)*MSe)
  
  lowICC <- (FL - 1) / (FL + (k - 1)) 
  upICC <- (FU - 1) / (FU + (k - 1))
  
  sol <- c(lowICC, ICC, upICC)
  names(sol) <- c("lower", "ICC(C,1)", "upper")
  
  return(sol)
  
}

# Function ICC.CK----
# calcule ICC (pas obligatoire)

ICC.CK.fct <- function(data.df, alpha) {
  # calculate the intra-class corelation, consistency version
  # data.df has the subj in row and the repeated measure in col
  # subj represent a random gouping factor
  
  data.ma <- as.matrix(data.df)
  
  n <- nrow(data.ma) # nbr of subj
  k <- ncol(data.ma) # nbr of repeated measures
  DFr <- n - 1
  DFe <- (n - 1)*(k - 1)
  
  Y.. <- mean(data.ma)
  Yi. <- apply(data.ma, 1, mean)
  Y.j <- apply(data.ma, 2, mean)
  
  # similar as aov(X ~ Subj + as.factor(Block), data=dataRF.df)
  SSr <- 0 # Sum of Square, Between row
  SSe <- 0 # Error
  for (i in 1:n) { # go through row
    Square <- (Yi.[i] - Y..)^2
    SSr <- SSr + Square
    for (j in 1:k) { # go through col
      Square <- (data.ma[i,j] - Yi.[i] - Y.j[j] + Y..)^2
      SSe <- SSe + Square
    }
  }
  SSr <- k*SSr
  
  MSr <- SSr / DFr # Mean Square, Between row
  MSe <- SSe / DFe # Error
  Fobs <- MSr / MSe
  
  FLtab <- qf(1-alpha/2, DFr, DFe)
  FUtab <- qf(1-alpha/2, DFe, DFr)
  
  FL <- Fobs / FLtab
  FU <- Fobs * FUtab
  
  ICC <- (MSr - MSe) / MSr
  
  lowICC <- 1 - 1 / FL
  upICC <- 1 - 1 /FU
  
  sol <- c(lowICC, ICC, upICC)
  names(sol) <- c("lower", "ICC(C,k)", "upper")
  
  return(sol)
  
}

# Function alphaDf----
# calcule les alpha pour toutes les dimension de varQConv selon param.df

"alphaDf.fct" <-
  function(qConv.cov, param.df) {
    nCol <- length(param.df)
    alpha.v <- NULL
    for (col in 1:nCol) {
      param.v <- param.df[col]
      iPos <- which(param.v == 1)
      iNeg <- which(param.v == -1)
      i <- c(iPos, iNeg)
      tmp.cov <- qConv.cov[i,i]
      nv <- ncol(tmp.cov)
      alpha <- (nv/(nv-1)) * (1 - sum(diag(tmp.cov)) / sum(tmp.cov))
      alpha.v <- c(alpha.v, alpha)
    }
    alpha.a <- as.array(alpha.v)
    dim(alpha.a) <- c(1, nCol)
    colnames(alpha.a) <- names(param.df)
    rownames(alpha.a) <- "Alpha"
    alpha.a
  }

# Function corSemiParCI----
# calcule cor semi partielle avec IC (pas oblgatoire)

corSemiParCI.fct <- function(Y, lm.0A, lm.0B, sign=+1) {
  # calc semi partial cor with CI for a predictor of a linear regression
  # A for the full model
  # B for a model with a subset of predictors
  # Mathieu d'Acremont
  # v1.0
  
  n <- length(Y)
  
  cor(Y, predict(lmA.tmp)) -> rOA
  cor(Y, predict(lmB.tmp)) -> rOB
  
  rSemi <- sqrt(rOA^2 - rOB^2)
  
  VarInf <- (rOA^4 - 2*rOA^2 + 1 - rOB^4 + rOB^2) / n
  
  lower.tmp <- rSemi - 1.96*sqrt(VarInf)
  upper.tmp <- rSemi + 1.96*sqrt(VarInf)
  
  if (sign==1) {
    res.tmp <- c(rSemi, lower.tmp, upper.tmp)
  }
  else {
    res.tmp <- c(-rSemi, -upper.tmp, -lower.tmp)
  }
  
  names(res.tmp) <- c("r", "lower", "upper")
  return(res.tmp)
}


# Function suParam----
# créer une fonction pour inverser les items de la SUS et calculer les scores composites

suParamG <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1)
suParam.df <- data.frame(suParamG)
save(file="output/suParam.df", suParam.df)


# Load data----

read.csv(file="data/quebec2020_jf_master.csv", row.names=1) -> bdd.df
names(bdd.df)
names(bdd.df[15:28])
su.df <- bdd.df[15:28]
names(su.df)

str(su.df)
summary(su.df)

# remove all NA
su.df <- na.omit(su.df) 
summary(su.df)
boxplot(su.df, notch = T)

# compute scores
min(su.df, na.rm=T); max(su.df, na.rm=T)
convScoreDf.fct(su.df, suParam.df$suParamG, 1, 7) -> suQConv.df
sumScoreDf.fct(suQConv.df, suParam.df) -> suVD.df

# descriptives----
summary(suQConv.df)
summary(suVD.df)

# item distribution----
round(skewness(su.df), 2)
summary(skewness(su.df))

round(kurtosis(su.df), 2)
summary(kurtosis(su.df))

# score G distribution----
round(skewness(suVD.df), 2)
round(kurtosis(suVD.df), 2)

# factorial analysis----
scree.plot(suQConv.df, sim=100)
scree.plot(suQConv.df, sim=10)
factanal(covmat=cov(suQConv.df), factors=1)

factanal(covmat=cov(suQConv.df), factors=2, rotation = "varimax")
factanal(covmat=cov(suQConv.df), factors=2, rotation = "promax")

factanal(covmat=cov(suQConv.df), factors=3, rotation = "varimax")
factanal(covmat=cov(suQConv.df), factors=3, rotation = "promax")


fit <- princomp(suQConv.df, cor=TRUE)
summary(fit) # print variance accounted for
loadings(fit) # pc loadings
plot(fit,type="lines") # scree plot
fit$scores # the principal components
biplot(fit) 

# PCA Variable Factor Map
library(FactoMineR)
result <- PCA(suQConv.df) # graphs generated automatically 


# => à vue de nez: garder items 7-8-5-14-11-13-(9)



# reliability----

# alpha----
round(alphaDf.fct(cov(suQConv.df, use="pairwise"), suParam.df), 2)
round(ICC.CK.fct(su.df, .05), 2)#For the reliability of the scale, the ICC(C,K) formula given by McGraw and Wong (1996) was implemented in R . The ICC(C,K) is equal to the alpha coefficient (Cronbach, 1951) but offer the advantage to have a CI. 


# internal consistency----

col2rowRF.fct(su.df) -> suRF.df
cor(su.df) -> su.cor
su.lower <- su.cor[lower.tri(su.cor)]

round(su.cor, 2) # Correlations between all items are positive

corTestDf.fct(su.df, rowSums(su.df)) -> tmp.li
formCor.fct(tmp.li)

#icc(su.df)
mean(su.lower) 
#ICC1.CI(X, Subj, suRF.df) 
#ICC2.CI(X, Subj, suRF.df) 
round(ICC.C1.fct(su.df, .05), 2)#intra-class correlation with 95% IC
