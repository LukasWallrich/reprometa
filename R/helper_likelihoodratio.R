# small sample AIC correction

Qc <- function(k1, k2, n) exp(k1 * (n/(n - k1 - 1)) - k2 * (n/(n - k2 -1)))


# change likelihood ratio to a variance ratio

LRtoVR <- function(lr, n, k1 = 2, k2 = 3, Correction = Qc) (lr/Correction(k1, k2, n))^(2/n)


# change variance ratio to d

VRtoD <- function(vr, n) 2 * sqrt(((n - 2) / n) * (vr - 1))


# change adjusted likelihood ratio to d

LRtoD <- function(lr, n, Correction = Qc) VRtoD(LRtoVR(lr, n, 2, 3, Correction), n)


# calculate d for the theoretically interesting effect

Dtie <- function(n_o, lr = 8, Correction = Qc) LRtoD(lr, n = n_o, Correction)


# change d to a ratio of unexplained variance

DtoVR <-  function(d, n) 1 + (n/(n - 2)) * (d/2)^2

# change variance ratio to an adjusted likelihood ratio

VRtoLR <- function(vr, n, k1 = 2, k2 = 3, Correction = Qc) Correction(k1, k2, n) * vr^(n/2)

# calculate the likelihood ratio comparing anticipated effect model to null model
#
LR.ae_0 <- function(es_r, Dtie, n_r) {
  VR0 <- DtoVR(es_r, n_r)
  VRae <- DtoVR(es_r - Dtie, n_r)
  return(VRtoLR( VR0/VRae , n = n_r, k1 = 2, k2 = 2))
}
#


###############################################################################
########### The following functions might be useful in the future ----
###############################################################################

# # AIC correction to the likelihood ratio
# #
# Qa <-
#   function(k1,k2,n) exp(k1-k2)
# #
#
# # BIC correction
# #
# Qb <-
#   function(k1,k2,n) (exp(k1-k2))^(log(n)/2)
#
# # calculate the probability of evidence for failure to replicate when the actual value is 0
# #
# pLRrep0 <-
#   function(lr,Dtie,n,dens=1000) {
#     ints<-seq(1/(2*dens),by=1/dens,length.out=dens-1)
#     (sum(LR.ae_0(2*qt(ints,n-2)/sqrt(n),Dtie,n)<lr)+1)/(dens)
#   }
# #
# # calculate the probability of evidence for replication when the actual value is Dtie
# #
# pLRrepAE <-
#   function(lr,Dtie,n,dens=1000) {
#     ints<-seq(1/(2*dens),by=1/dens,length.out=dens-1)
#     (sum(LR.ae_0(2*qt(ints,n-2,sqrt(n)*Dtie/2)/sqrt(n),Dtie,n)>lr)+1)/dens
#   }
# #
#
# # change delta AIC to an (AIC) adjusted likelihood ratio
# #
# DeltaAICtoLR <-
#   function(dAIC) exp(dAIC/2)
# #
# # change d to an adjusted likelihood ratio
# #
# DtoLR <-
#   function(d,n,Correction=Qc) VRtoLR(DtoVR(d,n),n,2,3,Correction)
# #
# # change d to t
# #
# DtoT <-
# function(d,n) sqrt(n)*d/2
# #
# # change f2 to a variance ratio
# #
# F2toVR <-
# function(f2) f2+1
# #
# # change F ratio to f2
# #
# FtoF2 <-
# function(F,df1,df2) F*df1/df2
# #
#
# # change likelihood ratio to f2
# #
# LRtoF2 <-
# function(lr,n,k1=2,k2=3,Correction=Qc) VRtoF2(LRtoVR(lr,n,k1,k2,Correction))
# #
# # calcuate power of t given an effect corresponding to a likelihood ratio
# #
# LRtoPower <-
# function(lr,n,alpha=0.05, Correction=Qc) {
# d<-LRtoD(lr,n,Correction=Correction)
# crit<-qt(alpha/2,n-2,lower.tail=F)
# ncp<-sqrt(n)*d/2
# return(pt(crit,n-2,ncp,lower.tail=F))
# }
# #
# # calculate the likelihood ratio comparing anticipated effect model to null model in terms of f2
# #
# LR.ae_0.f2 <-
# function(f2obt,f2ae,n,df1,repSign=1,Correction=Qc) {
# VR0<-F2toVR(f2obt)
# VRae<-F2toVR((sqrt(f2obt)-sign(repSign)*sqrt(f2ae))^2)
# return(VRtoLR(VR0/VRae,n,2,df1+1,Correction=Correction))
# }
# #
#
#
# #
# # change SSEs to a likelihood ratio
# #
# SSEtoLR <-
# function(SSE1,SSE2,n,k1=2,k2=3,Correction=Qc) VRtoLR(SSE1/SSE2,n,k1,k2,Correction)
# #
# # change t to d
# #
# TtoD <-
# function(t,n) 2*t/sqrt(n)
# #
#
# # change variance ratio to f2
# #
# VRtoF2 <-
# function(vr) vr-1
# #
