# R Script to translate Aid Zhantikin Meyers java data to R data
# 
# Author: Thaddeus Aid
# 07/08/2013
###############################################################################

# make sure the package is installed before running this script
library("rJava")

# initialise the JVM to 10GB
.jinit(parameters="-Xmx10240m")

# make a string array and create the object
s <- .jarray("string", "mainArgs")
ftt <- .jnew("FullTreeTransform", s)

# load the full tree log-likelihood results
#fullTreeLlr <- ftt$getFullTreeLlr()
#fullTreeLlr <- sapply(fullTreeLlr, .jevalArray)

# process the log-likelihood results into p-values
#fullTreep <- pchisq(fullTreeLlr, df = 1, lower.tail=FALSE)
#fullTreep <- -log10(fullTreep * 0.5)

# save the Full Tree p-values
#save(fullTreep, file="fullTreep.rData")

# save the Full Tree log-likelihood values
#save(fullTreeLlr, file="fullTreeLlr.rData")

# load Full Tree Lambda values
#fullTreeLambda <- ftt$getFullTreeLambda()
#fullTreeLambda <- sapply(fullTreeLambda, .jevalArray)

# Save full Tree Lambda Values
#save(fullTreeLambda, file="fullTreeLambda.rData")

# get the meta data
snpIndex = ftt$getSNPindex()
snpAncestry = ftt$getSNPancestry()
snpNames = ftt$getSNPnames()
snpCorrelations = ftt$getSNPcorrelations()

save(snpIndex, file="snpIndex.rData")
save(snpAncestry, file="snpAncestry.rData")
save(snpNames, file="snpNames.rData")
save(snpCorrelations, file="snpCorrelations.rData")

# get Position information
positions <- ftt$getPositions()

save(positions, file="positions.rData")

# load the snp information
snpLlr <- ftt$getSNPllr()
snpBranches <- ftt$getSNPbranches()
snpLambda <- ftt$getSNPlambda()
snpR2 <- ftt$getSNPr2()

# translate snp log-likelihood values to p-values
snpp <- pchisq(snpLlr, df = 1, lower.tail=FALSE)
snpp <- -log10(snpp)

# save snp data
save(snpLlr, file="snpLlr.rData")
save(snpBranches, file="snpBranches.rData")
save(snpLambda, file="snpLambda.rData")
save(snpR2, file="snpR2.rData")
save(snpp, file="snpp.rData")

# load max Tree data
treeMaxLlr <- ftt$getTreeMaxLlr()
treeMaxLambda <- ftt$getTreeMaxLambda()
treeMaxBranches <- ftt$getTreeMaxBranches()

# translate max tree log-likelihood data to p-values
treeMaxp <- pchisq(treeMaxLlr, df = 1, lower.tail=FALSE)
treeMaxp <- -log10(treeMaxp * 0.5)

# save max tree data
save(treeMaxLlr, file="treeMaxLlr.rData")
save(treeMaxLambda, file="treeMaxLambda.rData")
save(treeMaxBranches, file="treeMaxBranches.rData")
save(treeMaxp, file="treeMaxp.rData")

###
### Inverse!
###

# load the full tree log-likelihood results
#fullTreeLlrI <- ftt$getFullTreeLlrI()
#fullTreeLlrI <- sapply(fullTreeLlrI, .jevalArray)

# process the log-likelihood results into p-values
#fullTreepI <- pchisq(fullTreeLlrI, df = 1, lower.tail=FALSE)
#fullTreepI <- -log10(fullTreepI)

# save the Full Tree p-values
#save(fullTreepI, file="fullTreepI.rData")

# save the Full Tree log-likelihood values
#save(fullTreeLlrI, file="fullTreeLlrI.rData")

# load Full Tree Lambda values
#fullTreeLambdaI <- ftt$getFullTreeLambdaI()
#fullTreeLambdaI <- sapply(fullTreeLambdaI, .jevalArray)

# Save full Tree Lambda Values
#save(fullTreeLambdaI, file="fullTreeLambdaI.rData")

# load the snp information
snpLlrI <- ftt$getSNPllrI()
snpBranchesI <- ftt$getSNPbranchesI()
snpLambdaI <- ftt$getSNPlambdaI()
snpR2I <- ftt$getSNPr2I()
snpCorrelationsI = ftt$getSNPcorrelationsI()

# translate snp log-likelihood values to p-values
snppI <- pchisq(snpLlrI, df = 1, lower.tail=FALSE)
snppI <- -log10(snppI)

# save snp data
save(snpLlrI, file="snpLlrI.rData")
save(snpBranchesI, file="snpBranchesI.rData")
save(snpLambdaI, file="snpLambdaI.rData")
save(snpR2I, file="snpR2I.rData")
save(snppI, file="snppI.rData")
save(snpCorrelationsI, file="snpCorrelationsI.rData")

# load max Tree data
treeMaxLlrI <- ftt$getTreeMaxLlrI()
treeMaxLambdaI <- ftt$getTreeMaxLambdaI()
treeMaxBranchesI <- ftt$getTreeMaxBranchesI()

# translate max tree log-likelihood data to p-values
treeMaxpI <- pchisq(treeMaxLlrI, df = 1, lower.tail=FALSE)
treeMaxpI <- -log10(treeMaxpI)

# save max tree data
save(treeMaxLlrI, file="treeMaxLlrI.rData")
save(treeMaxLambdaI, file="treeMaxLambdaI.rData")
save(treeMaxBranchesI, file="treeMaxBranchesI.rData")
save(treeMaxpI, file="treeMaxpI.rData")


