#depending on the boundarz conditions
#add points to known points
#transforms the coordinates

#9. 11. 2017 Eliska
#okrajove podminky---------
#aktualizace 28. 6. 2018 E

zjisteniP = function(data, p1, p2, p3) {
  #zjisteni, jake uzovatel zadal podminky a pripadne provedeni operaci pridani bodu dna nebo hladiny nebo oboje
  #a pripadna tranformace souradnic
  #p1 chcu podminku dna
  #p2 chcu podminku hladiny
  #p3 chcu podmku transformace souradnic profilu
  #----vraci matici staniceni, hloubek svislic, hloubek merenych bodu a rychlosti v bodech

  #St = data$St
  Loc = data$Loc
  MeasD = data$MeasD
  Depth = data$Depth
  Vel = data$Vel

  mD = length(Loc) #pocet namerenych hodnot

  #normalizuju staniceni tak, aby zacinalo v nule ---
  Stan = Loc
  if (Stan[1] < Stan[2]){ #pro rostouci staniceni
    Stan = Stan-Stan[1]}
  if (Stan[1] > Stan[2]){ #pro klesajici staniceni
    Stan = Stan-Stan[length(Stan)]}
  #---------------------------


  if ((p1 == 1)&(p2 == 1)) {
    BCdno = BCriverbed(data) #podminka dno
    BChlad = BCsurface(data) #podminka hladina

    Stan = c(Stan,BCdno[,1], BChlad [,1])
    Depth = c(Depth, BCdno[,2], BChlad [,2])
    MeasD = c(MeasD, BCdno[,3], BChlad [,3])
    Vel = c(Vel, BCdno[,4], BChlad [,4])
  }

  if ((p1 == 0)&(p2 == 1)) {
    BChlad = BCsurface(data) #pouze podminka hladina

    Stan = c(Stan,BChlad [,1])
    Depth = c(Depth, BChlad [,2])
    MeasD = c(MeasD,  BChlad [,3])
    Vel = c(Vel, BChlad [,4])
  }

  if ((p1 == 1)&(p2 == 0)) {
    BCdno = BCriverbed(data) #pouze podminka dno

    Stan = c(Stan,BCdno[,1])
    Depth = c(Depth, BCdno[,2])
    MeasD = c(MeasD, BCdno[,3])
    Vel = c(Vel, BCdno[,4])
  }

  if (p3 == 1) { #"transformace" hloubek merenych bodu probiha
    A = NormMeasD(Depth, MeasD)
    Hlou = A$Hlou
    HlouBodu = A$HlouB
  }
  if (p3 == 0) { #"transformace" hloubek merenych bodu NEprobiha
    Hlou = Depth
    HlouBodu  = MeasD
  }

  return(list(Stan = Stan, Hlou = Hlou, HlouBodu = HlouBodu, Vel = Vel))
}


NormMeasD = function(aa, ba) {
  ##transformace merenych hloubek tak, aby srovnavaci rovina s nulovou hloubkou nebylo dno, ale hladina ---
  #vraci vektor transformovanych merenych hloubek

  Depth = aa
  MeasD = ba
  HlouB = (-1) * (Depth - MeasD)
  Hlou = (-1) * Depth

  return(list(Hlou = Hlou, HlouB = HlouB))
}


BCriverbed = function(data) {
  ##na dne je nulova rychlost ------------
  #vraci matici, jde 1. sloupec je staniceni, 2. sloupec je hloubka svislic,
  #3. sloupec je merena hloubka (= 0), 4. sloupec rychlost (= 0)

  #St = data$St
  Loc = data$Loc
  MeasD = data$MeasD
  Depth = data$Depth
  Vel = data$Vel

  #normalizuju staniceni tak, aby zacinalo v nule ---
  Stan = Loc
  if (Stan[1] < Stan[2]){ #pro rostouci staniceni
    Stan = Stan-Stan[1]}
  if (Stan[1] > Stan[2]){ #pro klesajici staniceni
    Stan = Stan-Stan[length(Stan)]}
  #---------------------------

  POM = matrix(0, nrow = length(Loc), ncol = 2)
  POM[,1] = Stan; POM[,2] = Depth
  POMbc = unique(POM)
  BCdno =  matrix(0, nrow = length(POMbc[,1]) , ncol = 4)
  BCdno[,1:2] = POMbc

  return(BCDno = BCdno)
}


BCsurface = function(data) {
  ##na hladine je stejna rychlost jak v nejvyssim merenem bode-------------
  #vraci matici, jde 1. sloupec je staniceni, 2. sloupec je hloubka svislic,
  #3. sloupec je merena hloubka (je shodna celkovou s hloubkou svislice),
  #4. sloupec rychlost (je shodna s nejvyse merenou rychlosti)

  #St = data$St
  Loc = data$Loc
  Depth = data$Depth
  MeasD = data$MeasD
  Vel = data$Vel

  #normalizuju staniceni tak, aby zacinalo v nule ---
  Stan = Loc
  if (Stan[1] < Stan[2]){ #pro rostouci staniceni
    Stan = Stan-Stan[1]}
  if (Stan[1] > Stan[2]){ #pro klesajici staniceni
    Stan = Stan-Stan[length(Stan)]}
  #---------------------------

  i = 1; j = 1;
  BCStan = 1; BCD = 1; BCVel = 1;
  while (i < length(Loc)) {
    podminka = 0;
    if (Loc[i]!=Loc[i+1]){
      podminka = 1
      BCStan[j] = Stan[i]; BCD[j] = Depth[i]; BCVel[j] = Vel[i] }
    if (Loc[i]==Loc[i+1]) {
      if (Loc[i+1] != Loc[i+2])	{
        podminka = 2
        if (MeasD[i] > MeasD[i+1]) #mereni od hladiny ke dnu
        {BCStan[j] = Stan[i]; BCD[j] = Depth[i]; BCVel[j] = Vel[i] }
        if (MeasD[i] < MeasD[i+1]) #mereni od dna k hladine
        {BCStan[j] = Stan[i+1]; BCD[j] = Depth[i+1]; BCVel[j] = Vel[i+1] }
      }
      if (Loc[i+1] == Loc[i+2] & Loc[i+2] != Loc[i+3]) {
        podminka = 3
        if (MeasD[i] > MeasD[i+1]) #mereni od hladiny ke dnu
        {BCStan[j] = Stan[i]; BCD[j] = Depth[i]; BCVel[j] = Vel[i] }
        if (MeasD[i] < MeasD[i+1]) #mereni od dna k hladine
        {BCStan[j] = Stan[i+2]; BCD[j] = Depth[i+2]; BCVel[j] = Vel[i+2] }
      }
      if (Loc[i+1] == Loc[i+2] & Loc[i+2] == Loc[i+3] & Loc[i+3] != Loc[i+4])
        #ctyri body ve svislici by nemely nastat, neni na to ani vzorec z normy
      {podminka = 4
      if (MeasD[i] > MeasD[i+1]) #mereni od hladiny ke dnu
      {BCStan[j] = Stan[i]; BCD[j] = Depth[i]; BCVel[j] = Vel[i] }
      if (MeasD[i] < MeasD[i+1]) #mereni od dna k hladine
      {BCStan[j] = Stan[i+3]; BCD[j] = Depth[i+3]; BCVel[j] = Vel[i+3] }
      }
      if (Loc[i+1] == Loc[i+2] & Loc[i+2] == Loc[i+3] & Loc[i+3] == Loc[i+4] & Loc[i+4] != Loc[i+5]) {
        podminka = 5
        if (MeasD[i] > MeasD[i+1]) #mereni od hladiny ke dnu
        {BCStan[j] = Stan[i]; BCD[j] = Depth[i]; BCVel[j] = Vel[i] }
        if (MeasD[i] < MeasD[i+1]) #mereni od dna k hladine
        {BCStan[j] = Stan[i+4]; BCD[j] = Depth[i+4]; BCVel[j] = Vel[i+4] }
      }
      if (Loc[i+1] == Loc[i+2] & Loc[i+2] == Loc[i+3] & Loc[i+3] == Loc[i+4] & Loc[i+4] == Loc[i+5] & Loc[i+5] != Loc[i+6]) {
        podminka = 6
        if (MeasD[i] > MeasD[i+1]) #mereni od hladiny ke dnu
        {BCStan[j] = Stan[i]; BCD[j] = Depth[i]; BCVel[j] = Vel[i] }
        if (MeasD[i] < MeasD[i+1]) #mereni od dna k hladine
        {BCStan[j] = Stan[i+5]; BCD[j] = Depth[i+5]; BCVel[j] = Vel[i+5] }
      }
    }
    if (podminka == 1){i = i + 1}
    if (podminka == 2){i = i + 2}
    if (podminka == 3){i = i + 3}
    if (podminka == 4){i = i + 4}
    if (podminka == 5){i = i + 5}
    if (podminka == 6){i = i + 6}
    j = j+1
  }
  BCStan[j] = Stan[i]; BCD[j] = Depth[i]; BCVel[j] = Vel[i]

  pom = sum(BCD == 0)
  mpom = length(BCD) - pom
  BChladina =  matrix(0, nrow = mpom , ncol = 4)
  j = 1; i = 1
  while (i <= (mpom+pom)) {
    if (BCD[i] != 0) {
      BChladina[j,1] = BCStan[i]
      BChladina[j,2] = BCD[i]
      BChladina[j,3] = BCD[i]
      BChladina[j,4] = BCVel[i]
      j = j + 1; i = i + 1}
    else {
      i = i + 1}
  }

  return(BChladina = BChladina)
}
