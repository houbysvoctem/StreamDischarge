#18. 5. 2018 Eliska

#require(data.table)


discharge_interpolace = function(mdata, interpolace){

    if ((interpolace == 1)|(interpolace == 2)|(interpolace == 3)|(interpolace == 4)){bunka = 0.12}
    if(interpolace == 5){bunka = 0.07}
   #vytvoreni site pro interpolaci
  A = sit_sample(mdata, bunka)
  SIT <- A$SIT
  POLYGONCS <- A$POLYGONCS

  #plot(SIT)
  #lines(POLYGONCS, col = 'red')

  #zjisteni hloubek ve svislicich site
  SIT_hl = HlouSvislic_sit(POLYGONCS, SIT)
  colnames(SIT_hl) <- c("stan", "hlou", "hlouBodu")

  #interpolace v bodech site
  if (interpolace == 1) {INTER = arith_mean (mdata, SIT_hl)}
  if (interpolace == 2) {INTER = inv_dis_w (mdata, SIT_hl)}
  if (interpolace == 3) {INTER = ord_kriging (mdata, SIT_hl)}

  #1. staniceni, 2. hloubka svslice, 3. hloubka bodu, 4. rychlost, 5.sirka bunky, 6.vyska bunky, 7.prutok bunkou
  CELKOVA = matrix(0, nrow = length(SIT[,1]), ncol = 7)
  CELKOVA[,1:4] = INTER

  #vyska bunek na nejvyssi horizontale:
  #1/2 bunky + hloubka bodu
  parhl = (-1) * max(CELKOVA[,3])
  pozice = 0
  pozice <- which((CELKOVA[,3])==max(CELKOVA[,3]))
  for (j in pozice) {
    CELKOVA[j,6] = parhl + bunka/2
  }

  #sirka bunek v prvni a posledni vertikale:
  #1/2 bunky + vzdalenost od pocatku staniceni(nebo vzdalenost od posledniho staniceni)
  pocatek = min(POLYGONCS[,1])
  parhl = min(CELKOVA[,1])
  pozice = 0
  pozice <- which((CELKOVA[,1])==min(CELKOVA[,1]))
  for (i in pozice) {
    CELKOVA[i,5] = (parhl-pocatek) + bunka/2
  }

  konec = max(POLYGONCS[,1])
  parhl = max(CELKOVA[,1])
  pozice = 0
  pozice <- which((CELKOVA[,1])==max(CELKOVA[,1]))
  for (i in pozice) {
    CELKOVA[i,5] = (konec-parhl) + bunka/2
  }

  pocMerBoduDna = length(POLYGONCS[,1])-1

  #vyska bunek u dna:
  #prunik vertikal s useckami znameho dna pro urceni spravne velikosti vysek bunek u dna
  #0. udelat matici usecek dna
  MAT_USECEK = matrix(0, nrow = pocMerBoduDna-1, ncol = 4) #1.x bodu 1, 2.y bodu 1, 3.x bodu 2, 4. y bodu 2
  for (j in 1:(pocMerBoduDna-1)) {
    MAT_USECEK [j,1:2] = POLYGONCS[j,] ; MAT_USECEK [j,3:4] = POLYGONCS[j+1,]
  }

  #1. najit bod u dna a jeho souradnice
  svislice = unique(SIT[,1:2])
  svisBodUDna = matrix(0,nrow = length(svislice[,1]), ncol = 2)
  for (j in 1:length(svislice[,1])) {
    pozice = 0
    pozice <- which(CELKOVA[,1]==svislice[j,1])
    POM = matrix(0,nrow = length(pozice), ncol = 2)
    for (k in 1:length(pozice)) {
      POM[k,1] = CELKOVA[pozice[k],1]
      POM[k,2] = CELKOVA[pozice[k],3]
    }
    svisBodUDna[j,] = c(POM[1,1],min(POM[,2]))
  }

  #pro kazdy bod u dna
  for (i in 1:length(svisBodUDna[,1])) {
    #zjistit hloubku ve svislici
    pozice = 0
    pozice <- which((SIT_hl[,1]==svisBodUDna[i,1])&(SIT_hl[,3]==svisBodUDna[i,2]))
    svisl_hl <- SIT_hl[pozice,2]
    #1/2 bunky + vzdalenost bodu u dna od bodu pruniku
    parhl = (-1)* (svisl_hl - svisBodUDna[i,2])
    pozice = 0
    pozice <- which((CELKOVA[,1]==svisBodUDna[i,1])&(CELKOVA[,3]==svisBodUDna[i,2]))
    CELKOVA[pozice,6] = parhl + bunka/2
  }
  #vyska a sirka ostatnich bunek = zadana velikost bunky
  for (i in (1:length(CELKOVA[,1]))) {
    if (CELKOVA[i,5] == 0) {CELKOVA[i,5] = bunka}
    if (CELKOVA[i,6] == 0) {CELKOVA[i,6] = bunka}
  }

  #prutok pro jednotlive bunky
  CELKOVA[,7] = CELKOVA[,5]*CELKOVA[,6]*CELKOVA[,4]

  #celkovy objem, prutok
  objem = sum(CELKOVA[,7])

  return(Q = objem)
}

