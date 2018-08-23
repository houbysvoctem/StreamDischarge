#creation of ortogonal regular grid in cross-secion
#26. 5. 2018 Eliska
#aktualizece 25. 7. 2018


sit_sample = function(mojedata, oko){
  library(data.table)
  library(sp)

  #vytvoreni polygonu prutocneho profilu
  #-----------------------------------

  #uprava staniceni na pocatek v nule
  if (mojedata$Loc[1] < mojedata$Loc[2]){ #pro rostouci staniceni
    Stanicenin = mojedata$Loc-mojedata$Loc[1]}
  if (mojedata$Loc[1] > mojedata$Loc[2]){ #pro klesajici staniceni
    Stanicenin = mojedata$Loc-mojedata$Loc[length(mojedata$Loc)]}

  #--osetreni pro pripad, kdy je merena hlouka u okraje
  podmhlubokybreh = 0
  if ((mojedata$Depth[1] > 0)&(mojedata$Depth[length(mojedata$Depth)] == 0) ) #pokud je v prvni svislici nejaka hloubka, pridej okraj
  {
    hloubkacs = c(0, mojedata$Depth)
    stanicenics = c(Stanicenin[1], Stanicenin)
    podmhlubokybreh = podmhlubokybreh + 1
  }
  if ((mojedata$Depth[length(mojedata$Depth)] > 0)&(mojedata$Depth[1] == 0)) #pokud je v posledni svislici nejaka hloubka, pridej okraj
  {
    hloubkacs = c(mojedata$Depth, 0)
    stanicenics = c(Stanicenin, Stanicenin[length(mojedata$Depth)])
    podmhlubokybreh = podmhlubokybreh + 1
  }
  if ((mojedata$Depth[length(mojedata$Depth)] > 0)&(mojedata$Depth[1] > 0)) #pokud je v prvni i posledni svislici nejaka hloubka, pridej okraj
  {
    hloubkacs = c(0, mojedata$Depth, 0)
    stanicenics = c(Stanicenin[1], Stanicenin, Stanicenin[length(mojedata$Depth)])
    podmhlubokybreh = podmhlubokybreh + 1
  }

  #vytvoreni polygonu pro pricny profil
  if (podmhlubokybreh == 0)
  {hloubkacs = -1*c(mojedata$Depth)
  stanicenics = Stanicenin}

  if (podmhlubokybreh > 0)
  {hloubkacs = -1*c(hloubkacs)
  stanicenics = stanicenics}

  Mat <- matrix(c(stanicenics, hloubkacs), nrow = length(hloubkacs), ncol = 2)
  Mat_unik <- unique(Mat)
  xym <- matrix(0,nrow = 1 + length(Mat_unik[,1]), ncol = 2)
  xym[1:length(Mat_unik[,1]),] <- Mat_unik
  xym[1+length(Mat_unik[,1]),] <- Mat_unik[1,]

  #aktualizace 24. 7. 2018, pouziti point.in.polygon
  x_max <- max(Mat_unik[,1]); y_min <- min(Mat_unik[,2])
  x_min <- min(Mat_unik[,1]); y_max <- max(Mat_unik[,2])
  #floor(x_max/oko)
  #ceiling(y_min/oko)

  ikska = seq(0,x_max,oko)
  ypsilonka = seq(0,y_min,(-1)*oko)
  ssbb = matrix(NA, ncol = 2, nrow = length(ikska)*length(ypsilonka))
  ind = 1; a = 1
  repeat{ #uprava 25. 7. 2018
    u = rep(ypsilonka[a],length(ikska))
    ssbb[ind : (ind + length(ikska)-1),1] = t(ikska)
    ssbb[ind : (ind + length(ikska)-1),2] = t(u)
    ind = ind + length(ikska)
    a = a + 1
    if(ind > length(ikska)*length(ypsilonka)){break}
  }

  orourke = point.in.polygon(ssbb[,1], ssbb[,2], xym[,1], xym[,2], mode.checked=FALSE)
  #zeberu jen: 1: point is strictly interior to pol
  souradnice <- matrix(NA, ncol = 2, nrow = length(orourke[orourke == 1]))
  ind = 0
  for (a in 1:length(orourke)) {
    if (orourke[a]==1){ind = ind + 1; souradnice[ind,] = ssbb[a,]}
  }
  #plot(souradnice)
  #lines(xym, col = 'red')

  return(list(SIT = souradnice, POLYGONCS = xym))

}


#zjisteni hloubky ve svislici site
HlouSvislic_sit = function(polygoncs, sit){
  SIT = matrix(0, nrow = length(sit[,1]), ncol = 3)
  SIT[,1] <- sit[,1]; SIT[,3] <- sit[,2]
  POLYGONCS <-polygoncs
  pocMerBoduDna = length(POLYGONCS[,1])-1

  #---------------------
  #prunik vertikal s useckami znameho dna pro urceni celkove hloubky ve vertikale site
  #0. udelat matici usecek dna
  MAT_USECEK = matrix(0, nrow = pocMerBoduDna-1, ncol = 4) #1.x bodu 1, 2.y bodu 1, 3.x bodu 2, 4. y bodu 2
  for (j in 1:(pocMerBoduDna-1)) {
    MAT_USECEK [j,1:2] = POLYGONCS[j,] ; MAT_USECEK [j,3:4] = POLYGONCS[j+1,]
  }

  #1. najit souradnice svislic
  svislice = unique(SIT[,1])

  #pro kazdou svislici
  for (i in 1:length(svislice)) {
    #2. najit usecku dna, kde jeji x_min < souradnice bodu < x_max
    j = 0
    repeat {
      j = j + 1
      if (((MAT_USECEK[j,1]< svislice[i])&(svislice[i] <= MAT_USECEK[j,3]))|(
        (MAT_USECEK[j,1]>= svislice[i])&(svislice[i] > MAT_USECEK[j,3]))){break} #aktualizece 21. 7. 2018
    }
    #3. najit obecny tvar primky dna a primky vertikaly
    #3.a) smerovy vektor primky
    smerovy_x = MAT_USECEK[j,3] - MAT_USECEK[j,1]
    smerovy_y = MAT_USECEK[j,4] - MAT_USECEK[j,2]
    #3.b) normalovy vektor primky
    normalovy_x = smerovy_y ; normalovy_y = (-1)*smerovy_x
    #3.c) dopocteni cecka v obecnem tvaru primky
    c = (-1) * (normalovy_x * MAT_USECEK[j,1] + normalovy_y * MAT_USECEK[j,2] )
    #3.d) obecny tvar vertikaly je x=staniceni
    prunik_x = svislice[i]
    #4. prunik usecky dna s useckou vertikaly: normalovy_x * prunik_x + normalovy_y * prunik_y + c = 0
    #do obecne rce dna doplnim vertikalu (x=staniceni), ziskam souradnici hloubky
    prunik_y = ((-1)*c + (-1)*normalovy_x * prunik_x )/ normalovy_y

    #doplnim hloubku ke vsem bodum ve svislici
    pozice = 0
    pozice <- which(SIT[,1]==svislice[i])
    SIT[pozice,2] = prunik_y
  }
  #------------------------------

  return(SIT_s_hloubkami=SIT)

}
