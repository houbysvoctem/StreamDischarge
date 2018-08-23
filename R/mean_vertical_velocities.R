#22. 6. 2015 Eliska
#update pro rozdil mereni od dna k hladine a naopak 16. 10. 2017
#vypocet strednich svislicovych rychloLoci
#last update 28. 6. 2018 E

#calculation of mean velocities in verticals of cross-section by CSN/EU norm
#data muLoc be in format: data.frame
#in data there muLoc be columns: Loc,Depth,MeasD,Vel
#data <- read.table('my_data.txt' , header = T, sep = "", dec = ".")

HorizontalV <- function(data){

  Loc = data$Loc
  MeasD = data$MeasD
  Depth = data$Depth
  Vel = data$Vel

  i = 1; j = 1;SVel = 1;SH = 1
  while (i < length(Loc)) {
    podminka = 0;

    if (Loc[i]!=Loc[i+1]){SVel[j] = Vel[i]; podminka = 1
    SH[j] = Depth[i]}
    if (Loc[i]==Loc[i+1]) {
      if (Loc[i+1] != Loc[i+2])
      {SVel[j] = 0.5*(Vel[i]+Vel[i+1])
      podminka = 2
      SH[j] = Depth[i]}
      if (Loc[i+1] == Loc[i+2] & Loc[i+2] != Loc[i+3])
      {SVel[j] = 0.25*(Vel[i]+2*Vel[i+1]+Vel[i+2])
      podminka = 3
      SH[j] = Depth[i]}
      if (Loc[i+1] == Loc[i+2] & Loc[i+2] == Loc[i+3] & Loc[i+3] != Loc[i+4])
        #ctyri body ve svislici by nemely naLocat, neni na to ani vzorec z normy
      {if (MeasD[i] > MeasD[i+1]) #mereni od hladiny ke dnu
      {SVel[j] = 0.2*(Vel[i]+2*Vel[i+1]+Vel[i+2]+Vel[i+3])}
        if (MeasD[i] < MeasD[i+1]) #mereni od dna k hladine
        {SVel[j] = 0.2*(Vel[i]+Vel[i+1]+2*Vel[i+2]+Vel[i+3])}
        podminka = 4
        SH[j] = Depth[i]}
      if (Loc[i+1] == Loc[i+2] & Loc[i+2] == Loc[i+3] & Loc[i+3] == Loc[i+4] & Loc[i+4] != Loc[i+5])
      {if (MeasD[i] > MeasD[i+1]) #mereni od hladiny ke dnu
      {SVel[j] = 0.1*(Vel[i]+3*Vel[i+1]+3*Vel[i+2]+2*Vel[i+3]+Vel[i+4])}
        if (MeasD[i] < MeasD[i+1]) #mereni od dna k hladine
        {SVel[j] = 0.1*(Vel[i]+2*Vel[i+1]+3*Vel[i+2]+3*Vel[i+3]+Vel[i+4])}
        podminka = 5
        SH[j] = Depth[i]}
      if (Loc[i+1] == Loc[i+2] & Loc[i+2] == Loc[i+3] & Loc[i+3] == Loc[i+4] & Loc[i+4] == Loc[i+5] & Loc[i+5] != Loc[i+6])
      {SVel[j] = 0.1*(Vel[i]+2*Vel[i+1]+2*Vel[i+2]+2*Vel[i+3]+2*Vel[i+4]+Vel[i+5])
      podminka = 6
      SH[j] = Depth[i]}
    }
    if (podminka == 1){i = i + 1}
    if (podminka == 2){i = i + 2}
    if (podminka == 3){i = i + 3}
    if (podminka == 4){i = i + 4}
    if (podminka == 5){i = i + 5}
    if (podminka == 6){i = i + 6}
    j = j+1
  }
  SVel[j] = Vel[i]
  SVel = round(SVel,4)
  SH[j] = Depth[i]

  return(list(SVel = SVel, SH = SH))
}

