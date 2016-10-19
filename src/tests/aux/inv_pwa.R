library(TraME)

a = 2
B = matrix(1,1,1) + 1
C = matrix(1,1,1)
k=1

inversePWA(a,B,C,k)

.Call("inv_pwa_R", a,B,C,k, PACKAGE = "TraME")$vals
