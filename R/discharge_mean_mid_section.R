#22. 6. 2015 Eliska
#vypocet prutoku metodou svislisovych nebo mezisvislicovych pasu
#last update 28. 6. 2018

#calculation of discharge with method of mid/mean section by CNS/EU norm
#data must be in format: data.frame
#in data there must be columns: Loc,Depth,MeasD,Vel
#data <- read.table('my_data.txt' , header = T, sep = "", dec = ".")
#need of a function HorizontalV

#------------------------------------------
#-----svislicove pasy -- mid-section segment method
#------------------------------------------

SvisP = function (data){

  Loc = data$Loc
  MeasD = data$MeasD
  Depth = data$Depth

  a = HorizontalV(data)
  SVel = a$SVel
  SH = a$SH

  #---uprava staniceni na pocatek v nule
  Staniceni = unique(Loc)
  if (Staniceni[1] < Staniceni[2]){ #pro rostouci staniceni
    Staniceni = Staniceni-Staniceni[1]}
  if (Staniceni[1] > Staniceni[2]){ #pro klesajici staniceni
    Staniceni = Staniceni-Staniceni[length(Staniceni)]}
  #-------------------------
  #----vypocet becek, neboli sirek svislicovych pasu
  SBe = vector(length = length(Staniceni)-1 , mode = 'numeric')
  if (Staniceni[1] < Staniceni[2]){ #pro rostouci staniceni
    SBe[1] = 0 #(Staniceni[2]-Staniceni[1])/2 --se dle normy zanedvaba
    SBe[length(Staniceni)] = 0 #(Staniceni[length(Staniceni)] - Staniceni[length(Staniceni)-1])/2 --se dle normy zanedvaba
    for (j in 2:(length(Staniceni)-1))
    {SBe[j] = (Staniceni[j+1]-Staniceni[j-1])/2}
  }
  if (Staniceni[1] > Staniceni[2]){ #pro klesajici staniceni
    SBe[1] = 0 #(Staniceni[1]-Staniceni[2])/2 --se dle normy zanedvaba
    SBe[length(Staniceni)] = 0 #(Staniceni[length(Staniceni)-1] - Staniceni[length(Staniceni)])/2 --se dle normy zanedvaba
    for (j in 2:(length(Staniceni)-1))
    {SBe[j] = (Staniceni[j-1]-Staniceni[j+1])/2}
  }
  #-----------------------------

  #vypocet prutoku
  SQ = SVel*SH*SBe
  SQ = sum(SQ)

  return(list(SQ = SQ))
}

#-------------------------------------------
#------mezisvislicove pasy -- mean-section segment method
#-------------------------------------------
#m je pro mista, kde hloubka a rychlost u brehu nejsou nulove, hodnota obvykle 5 az 7
#m=4 pro hrube dno, m=10 pro hladke dno

MeziSvisP = function (data, m){

  Loc = data$Loc
  MeasD = data$MeasD
  Depth = data$Depth

  a = HorizontalV(data)
  MSVel = a$SVel
  MSH = a$SH

  #---uprava staniceni na pocatek v nule
  Staniceni = unique(Loc)
  if (Staniceni[1] < Staniceni[2]){ #pro rostouci staniceni
    Staniceni = Staniceni-Staniceni[1]}
  if (Staniceni[1] > Staniceni[2]){ #pro klesajici staniceni
    Staniceni = Staniceni-Staniceni[length(Staniceni)]}
  #-------------------------

  #----------------------------------
  #----vypocet becek, neboli sirek mezisvislicovych pasu
  MSBe = vector(length = length(Staniceni)-1 , mode = 'numeric')
  if (Staniceni[1] < Staniceni[2]){ #pro rostouci staniceni
    for (j in 1:(length(Staniceni)-1))
    {MSBe[j] = (Staniceni[j+1]-Staniceni[j])}}
  if (Staniceni[1] > Staniceni[2]){ #pro klesajici staniceni
    for (j in 1:(length(Staniceni)-1))
    {MSBe[j] = (Staniceni[j]-Staniceni[j+1])}}
  #-----------------------------

  #vypocet prutoku

  MSQ1 = vector(mode = 'numeric', length = length(MSBe))
  for (j in 2:(length(MSBe)-1))
  {MSQ1[j] = ((MSVel[j]*MSH[j] + MSVel[j+1]*MSH[j+1]) * (MSBe[j]))/2}
  MSQ1[1] = (MSVel[2]*m/(m+1))*MSH[1]*MSBe[1]/2
  MSQ1[length(MSBe)] = (MSVel[length(MSVel)-1]*m/(m+1))*MSH[length(MSH)]* MSBe[length(MSBe)]/2
  MSQ = sum(MSQ1)

  return(list(MSQ = MSQ))
}

