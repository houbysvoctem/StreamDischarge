#solution of the case, when the singular matrix appear dirung ordinary kriging
#13. 5. 2018 Eliska
#update 21. 7. 2018



#-----------------
#ordinary kriging

singularita_kriging <- function(m, Ecko, EL, i, mynugget, mysill, myrange){
  #m je pocet bodu, ze kterych ma probihat interpolace
  #Ecko je setrizena vstupni matice podle vzdalenosti
  #EL je vstupni matice

  n = m
  E = Ecko
  pomCI = EL
  ndos = 0; singu = 0

  repeat{
    KVV = matrix(0, nrow = n, ncol = n)  #kovariance matrix samples/samples
    for (k in 2:(n+1)) {
      di = ((pomCI[k,1] - pomCI[2:(n+1),1])^2 + (pomCI[k,3] - pomCI[2:(n+1),3])^2) ^(1/2)
      o = k - 1
      for (j in 1:n) {
        if (di[j] != 0) {KVV[o,j] = mynugget + mysill*(1-exp(-(di[j]^2)/(myrange^2)))}
        if (di[j] == 0) {KVV[o,j] = 0}
      }
    }

    J = matrix(1, nrow = n, ncol = 1)

    L = matrix(0, nrow = n+1, ncol = n+1)
    L [1:n,1:n] = KVV
    L[n+1,1:n] = t(J)
    L[1:n,n+1] = J

    #nekdy se stane, ze L je singularni, pak nejde udelat inverzni matici a vypocitat vahy
    if (testinv(L) == TRUE) { #kdyz L neni singularni
      KSV = matrix(0, nrow = n, ncol = 1) #kovariance matrix samples/volume to be estimated
      di = ((E[i,1] - pomCI[2:(n+1),1])^2 + (E[i,3] - pomCI[2:(n+1),3])^2) ^(1/2)
      for (j in 1:n) {
        if (di[j] != 0) {KSV[j,1] = mynugget + mysill*(1-exp(-(di[j]^2)/(myrange^2)))}
        if (di[j] == 0) {KSV[j,1] = 0}
      }

      KSVv = matrix(1, nrow = n+1, ncol = 1)
      KSVv[1:n,] = KSV

      wa = solve(L) %*% KSVv
      w = wa[1:n,1]; Lagrange = wa[n+1,1]

      pon = sum ( w * pomCI[2:(n+1),4])

      singu = 0
      #print('matice neni singularni')
    }
    if (testinv(L) == FALSE) {#kdyz L je singularni
      singu = 1
      if ((n < length(E[,1]))&(ndos == 0)) {n = n + 1 #pokud je n mensi nez pocet merenych bodu, pridej n
      } else if ((n == length(E[,1]))|(ndos == 1)) {n = n - 1 ; ndos = 1 #pokud je n nebo byl roven poctu merenych bodu, uber n
      } else if (n == 0) {ndos = 2} #pokud je n nula, nedelej nic
      #print('matice je singularni'); print(n)
    }
    if ((singu == 0)|(ndos == 2)) {break}
    n = n + 1
  }

  return(list(pon = pon, podminkan = ndos))

}
