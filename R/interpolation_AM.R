#27. 5. 2018 Eliska
#interpolation with arithmetical mean

arith_mean = function(data, sit) {
  #uses function from "okrajove_podminky"

  #--parameters of interpolation:
  a = 1.71 #a = main half-axis of ellipse
  n = 2 #n = number of neighbpuring points from that the interpolation is done
  p1 = 0 #p1 = in the intersection of measured verticals and the river streambed the velocity is assumed 0
  p2 = 1 #p2 = in the intersection of measured verticals and the water surface the velocita is assumed equal to the velocity of the top measured point velocity
  p3 = 0 #p3 = the "coordinates" of measured points are transformed to normal coordinates

  b = 1 #lateral half-axis of ellipse
  #---------

  VYSL = matrix(0, ncol = 4, nrow = length(sit[,1])) #sem se posypou vysledky
  colnames(VYSL) <- c('stan', 'hlou', 'hlouBodu', 'vysl')
  VYSL[,1:3] = sit

  Loc = data$Loc
  mD = length(Loc) #pocet namerenych hodnot

  A = zjisteniP (data, p1, p2, p3)
  Stan = A$Stan
  Hlou = A$Hlou
  HlouBodu = A$HlouBodu
  Vel = A$Vel

  Ep = matrix(0, nrow = length(Vel), ncol = 5) #E je spojeni matic znamych hodnot
  Ep[,1] = Stan #1. sloupec je staniceni
  Ep[,2] = Hlou #2. sloupec jsou hloubky svislic
  Ep[,3] = HlouBodu #3. sloupec jsou hloubky bodu
  Ep[,4] = Vel #4. sloupec jsou rychlosti v bodech
  #5. sloupec budou vzdalenosti

  #odstranit hodnoty, kde nebyla merena rychlost
  #(hloubka bodu je 0 a zaroven rychlost je 0)nebo(hloubka svislice se rovna hloubce mereneho bodu)
  poc = sum(((Ep[1:mD,3] == 0)&(Ep[1:mD,4] == 0))|(Ep[1:mD,2]==Ep[1:mD,3]))
  i = 1; j = 1
  E = matrix(0, nrow = length(Vel)-poc, ncol = 5)
  while (i <= mD) {
    if (((Ep[i,3] == 0)&(Ep[i,4] == 0))|(Ep[i,2]==Ep[i,3])) {
      i = i + 1}
    else {
      E[j,] = Ep[i,]
      i = i + 1; j = j + 1
    }
  }
  #pridam zbytek hodnot do matice E - body doplnene, pokud byly podminky p1, p2
  if ((p1 == 1)|(p2==1)) { E[j:(length(Vel)-poc),] = Ep[(j+poc):length(Vel),] }

  if (p3 == 0) { #pokud nebyla podminka tranformace souradnic, musim prevest sit na "merene souradnice"
    M <- NormMeasD(sit[,2], sit[,3])
    sit[,2] <- M$Hlou
    sit[,3] <- M$HlouB
  }

  for (i in 1:(length(sit[,1])) ) { #pro kazdy bod site postupne

    DIST = 0 #budouci vektor vzdalenosti vybraneho bodu od vsech ostatnich
    DIST = ((E[,1] - sit[i,1])^2)/(a^2)+((E[,3] - sit[i,3])^2)/(b^2) #vzdalenosti v intencich elipsy
    E[,5] = DIST
    pom = E[order(E[,5],decreasing = FALSE),] #setrizeni
    pon = (sum(pom[2:(n+1),4]))/n
    VYSL[i,4] = (round(pon*10000))/10000 #do radku vypoctene rychlosti pro jednotlive body
  }

  return(D = VYSL)
}
