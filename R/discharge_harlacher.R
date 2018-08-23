#16. 10. 2017 Eliska
#vypocet prutoku metodou Harlachera - "grafickou"

#calculation of discharge with graphical Harlacher method
#siutable for more points in vertical due to linear aproximation betwen points - not spline
#data must be in format: data.frame
#in data there must be columns: Loc,Depth,MeasD,Vel
#data <- read.table('my_data.txt' , header = T, sep = "", dec = ".")

Harlacher = function (data, m){

  Loc = data$Loc
  MeasD = data$MeasD
  Depth = data$Depth
  Vel = data$Vel

  #--overeni pro pripad, kdy je merena hloubka u okraje
  podmhlubokybreh = 0
  if (data$Depth[1] > 0) #pokud je v prvni svislici nejaka hloubka
  {podmhlubokybreh = 1}
  if (data$Depth[length(data$Depth)] > 0) #pokud je v posledni svislici nejaka hloubka
  {podmhlubokybreh = 2}

  #-----------------------------------
  #-----jednotlive plochy pro svislice
  i = 1; j = 1;k=0;PlochaSvislic = 0;SvislicoveH = 0;
  Staniceni = vector(length =length(unique(Loc))-1, mode = 'numeric');

  while (i < length(Loc)) {
    vyska = 0; podminka = 0;
    PlochaA = 0; PlochaB = 0; PlochaC = 0;
    Plocha = vector(length =5, mode = 'numeric')

    if (Loc[i]!=Loc[i+1]){ #mereni v 1 bode nebo v zadnem bode
      if ((i == 1)&(podmhlubokybreh == 1)) #mereni v pripade hlubokeho brehu
      {Vel[i] = Vel[i+1]*m/(m+1)}
      PlochaSvislic[j] = (MeasD[i]*Vel[i])/2 + #trojuhelnik
        (Depth[i]-MeasD[i])*Vel[i] #obdelnik
      podminka = 1
      SvislicoveH[j] = Depth[i]
    }

    if (Loc[i]==Loc[i+1]) {
      if (Loc[i+1] != Loc[i+2]) #mereni ve 2 bodech
      {if (MeasD[i] > MeasD[i+1]) #mereni od hladiny ke dnu
      {PlochaA = (Depth[i]-MeasD[i])*Vel[i] #obdelnik u hladiny
      PlochaC =  (MeasD[i+1]*Vel[i+1])/2 #trojuhlenik u dna
      vyska = MeasD[i] - MeasD[i+1]
      }
        if (MeasD[i] < MeasD[i+1]) #mereni od dna k hladine
        {PlochaA = (Depth[i+1]-MeasD[i+1])*Vel[i+1] #obdelnik u hladiny
        PlochaC =  (MeasD[i]*Vel[i])/2 #trojuhlenik u dna
        vyska = MeasD[i+1] - MeasD[i]
        }

        if (Vel[i] > Vel[i+1]) {
          PlochaB = (Vel[i] - Vel[i+1])*vyska/2  + #trojuhlenik
            (Vel[i+1]*vyska)} #obdelnik
        if (Vel[i] < Vel[i+1]) {
          PlochaB = (Vel[i+1] - Vel[i])*vyska/2  + #trojuhlenik
            (Vel[i]*vyska)} #obdelnik

        PlochaSvislic[j] = PlochaA + PlochaB + PlochaC
        podminka = 2
        SvislicoveH[j] = Depth[i]}

      if (Loc[i+1] == Loc[i+2] & Loc[i+2] != Loc[i+3]) #mereni ve 3 bodech
      {if (MeasD[i] > MeasD[i+1]) #mereni od hladiny ke dnu
      {PlochaA = (Depth[i]-MeasD[i])*Vel[i] #obdelnik u hladiny
      PlochaC =  (MeasD[i+2]*Vel[i+2])/2 #trojuhlenik u dna

      t = 1
      for (k in i:(i+1)) {vyska[t] = MeasD[k] - MeasD[k+1] ; t = t +1}
      }
        if (MeasD[i] < MeasD[i+1]) #mereni od dna k hladine
        {PlochaA = (Depth[i+2]-MeasD[i+2])*Vel[i+2] #obdelnik u hladiny
        PlochaC =  (MeasD[i]*Vel[i])/2 #trojuhlenik u dna
        t = 1
        for (k in i:(i+1)) {vyska[t] = MeasD[k+1] - MeasD[k] ; t = t +1}
        }

        t = 1
        for (k in i:(i+1)) {
          if (Vel[k] > Vel[k+1]) {
            Plocha[t] = (Vel[k] - Vel[k+1])*vyska[t]/2  + #trojuhlenik
              (Vel[k+1]*vyska[t])} #obdelnik
          if (Vel[k] < Vel[k+1]) {
            Plocha[t] = (Vel[k+1] - Vel[k])*vyska[t]/2  + #trojuhlenik
              (Vel[k]*vyska[t])} #obdelnik
          t = t +1
        }

        PlochaB = sum(Plocha)
        PlochaSvislic[j] = PlochaA + PlochaB + PlochaC
        podminka = 3
        SvislicoveH[j] = Depth[i]}

      if (Loc[i+1] == Loc[i+2] & Loc[i+2] == Loc[i+3] & Loc[i+3] != Loc[i+4]) #mereni ve 4 bodech
      {if (MeasD[i] > MeasD[i+1]) #mereni od hladiny ke dnu
      {PlochaA = (Depth[i]-MeasD[i])*Vel[i] #obdelnik u hladiny
      PlochaC =  (MeasD[i+3]*Vel[i+3])/2 #trojuhlenik u dna

      t = 1
      for (k in i:(i+2)) {vyska[t] = MeasD[k] - MeasD[k+1] ; t = t +1}
      }
        if (MeasD[i] < MeasD[i+1]) #mereni od dna k hladine
        {PlochaA = (Depth[i+3]-MeasD[i+3])*Vel[i+3] #obdelnik u hladiny
        PlochaC =  (MeasD[i]*Vel[i])/2 #trojuhlenik u dna
        t = 1
        for (k in i:(i+2)) {vyska[t] = MeasD[k+1] - MeasD[k] ; t = t +1}
        }

        t = 1
        for (k in i:(i+2)) {
          if (Vel[k] > Vel[k+1]) {
            Plocha[t] = (Vel[k] - Vel[k+1])*vyska[t]/2  + #trojuhlenik
              (Vel[k+1]*vyska[t])} #obdelnik
          if (Vel[k] < Vel[k+1]) {
            Plocha[t] = (Vel[k+1] - Vel[k])*vyska[t]/2  + #trojuhlenik
              (Vel[k]*vyska[t])} #obdelnik
          t = t +1
        }

        PlochaB = sum(Plocha)
        PlochaSvislic[j] = PlochaA + PlochaB + PlochaC
        podminka = 4
        SvislicoveH[j] = Depth[i]}

      if (Loc[i+1] == Loc[i+2] & Loc[i+2] == Loc[i+3] & Loc[i+3] == Loc[i+4] & Loc[i+4] != Loc[i+5])  #mereni ve 5 bodech
      {if (MeasD[i] > MeasD[i+1]) #mereni od hladiny ke dnu
      {PlochaA = (Depth[i]-MeasD[i])*Vel[i] #obdelnik u hladiny
      PlochaC =  (MeasD[i+4]*Vel[i+4])/2 #trojuhlenik u dna

      t = 1
      for (k in i:(i+3)) {vyska[t] = MeasD[k] - MeasD[k+1] ; t = t +1}
      }
        if (MeasD[i] < MeasD[i+1]) #mereni od dna k hladine
        {PlochaA = (Depth[i+4]-MeasD[i+4])*Vel[i+4] #obdelnik u hladiny
        PlochaC =  (MeasD[i]*Vel[i])/2 #trojuhlenik u dna
        t = 1
        for (k in i:(i+3)) {vyska[t] = MeasD[k+1] - MeasD[k] ; t = t +1}
        }

        t = 1
        for (k in i:(i+3)) {
          if (Vel[k] > Vel[k+1]) {
            Plocha[t] = (Vel[k] - Vel[k+1])*vyska[t]/2  + #trojuhlenik
              (Vel[k+1]*vyska[t])} #obdelnik
          if (Vel[k] < Vel[k+1]) {
            Plocha[t] = (Vel[k+1] - Vel[k])*vyska[t]/2  + #trojuhlenik
              (Vel[k]*vyska[t])} #obdelnik
          t = t +1
        }

        PlochaB = sum(Plocha)
        PlochaSvislic[j] = PlochaA + PlochaB + PlochaC
        podminka = 5
        SvislicoveH[j] = Depth[i]}

      if (Loc[i+1] == Loc[i+2] & Loc[i+2] == Loc[i+3] & Loc[i+3] == Loc[i+4] & Loc[i+4] == Loc[i+5] & Loc[i+5] != Loc[i+6])  #mereni ve 6 bodech
      {if (MeasD[i] > MeasD[i+1]) #mereni od hladiny ke dnu
      {PlochaA = (Depth[i]-MeasD[i])*Vel[i] #obdelnik u hladiny
      PlochaC =  (MeasD[i+5]*Vel[i+5])/2 #trojuhlenik u dna

      t = 1
      for (k in i:(i+4)) {vyska[t] = MeasD[k] - MeasD[k+1] ; t = t +1}
      }
        if (MeasD[i] < MeasD[i+1]) #mereni od dna k hladine
        {PlochaA = (Depth[i+5]-MeasD[i+5])*Vel[i+5] #obdelnik u hladiny
        PlochaC =  (MeasD[i]*Vel[i])/2 #trojuhlenik u dna

        t = 1
        for (k in i:(i+4)) {vyska[t] = MeasD[k+1] - MeasD[k] ; t = t +1}
        }

        t = 1
        for (k in i:(i+4)) {
          if (Vel[k] > Vel[k+1]) {
            Plocha[t] = (Vel[k] - Vel[k+1])*vyska[t]/2  + #trojuhlenik
              (Vel[k+1]*vyska[t])} #obdelnik
          if (Vel[k] < Vel[k+1]) {
            Plocha[t] = (Vel[k+1] - Vel[k])*vyska[t]/2  + #trojuhlenik
              (Vel[k]*vyska[t])} #obdelnik
          t = t +1
        }

        PlochaB = sum(Plocha)
        PlochaSvislic[j] = PlochaA + PlochaB + PlochaC
        podminka = 6
        SvislicoveH[j] = Depth[i]}

    }
    if (podminka == 1){i = i + 1}
    if (podminka == 2){i = i + 2}
    if (podminka == 3){i = i + 3}
    if (podminka == 4){i = i + 4}
    if (podminka == 5){i = i + 5}
    if (podminka == 6){i = i + 6}
    j = j+1
  }

  #vyskoci ven, pokud posledni svislice muze mit nula nebo jeden mereny bod - aspon doufam ;)
  #a pro ten pripad jsou tu nasledujici radky
  if ((i == length(Loc))&(podmhlubokybreh == 2)) #mereni v pripade hlubokeho brehu
  {Vel[i] = Vel[i-1]*m/(m+1)}

  PlochaSvislic[j] = (MeasD[i]*Vel[i])/2 + #trojuhelnik
    (Depth[i]-MeasD[i])*Vel[i] #obdelnik
  SvislicoveH[j] = Depth[i]
  #-----------------------------------------------------

  #------------------------------
  #---uprava staniceni na pocatek v nule
  Staniceni = unique(Loc)
  if (Staniceni[1] < Staniceni[2]){ #pro rostouci staniceni
    Staniceni = Staniceni-Staniceni[1]}
  if (Staniceni[1] > Staniceni[2]){ #pro klesajici staniceni
    Staniceni = Staniceni-Staniceni[length(Staniceni)]}
  #-------------------------

  #----------------------------------
  #----vypocet becek, neboli sirek mezisvislicovych pasu
  SBe = vector(length = length(Staniceni)-1 , mode = 'numeric')
  if (Staniceni[1] < Staniceni[2]){ #pro rostouci staniceni
    for (j in 1:(length(Staniceni)-1))
    {SBe[j] = (Staniceni[j+1]-Staniceni[j])}}
  if (Staniceni[1] > Staniceni[2]){ #pro klesajici staniceni
    for (j in 1:(length(Staniceni)-1))
    {SBe[j] = (Staniceni[j]-Staniceni[j+1])}}
  #-----------------------------

  #---------------------------------------------------
  #--vypocet prutoku,
  #tj. graficke plochy tvorene vynesenymi plochami svislic nad staniceni profilu
  Prutok = vector(length = length(SBe), mode = 'numeric')
  for (j in 1:(length(PlochaSvislic)-1)) {
    if (PlochaSvislic[j] < PlochaSvislic[j+1]) {
      Prutok[j] = (PlochaSvislic[j+1] - PlochaSvislic[j])*SBe[j]/2 + #trojuhlenik
        PlochaSvislic[j]*SBe[j] #obdelnik
    }
    if (PlochaSvislic[j] > PlochaSvislic[j+1]) {
      Prutok[j] = (PlochaSvislic[j] - PlochaSvislic[j+1])*SBe[j]/2 + #trojuhlenik
        PlochaSvislic[j+1]*SBe[j] #obdelnik
    }
  }
  Prutokcely = sum(Prutok)

  PlochaSvislic
  SvislicoveH
  SBe
  Prutok
  Prutokcely

  return(list(HarQ = Prutokcely))
}
