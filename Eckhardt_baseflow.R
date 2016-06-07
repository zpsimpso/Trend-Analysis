Eckhardt <- function(streamflow, filter_parameter=0.925, BFI_max = 0.80, passes=1){
  #this function is shamelessly based on the function for BaseflowSeparation in the
  #EcoHydRology package. I'm merely adjusting it according to Eckhardt (2005) by updating the filter.
  #Nathan and McMahon (1990) suggest a = 0.925, and Eck. (2005) shows it's not terribly sensitive
  #Though, it can be more precisely estimated using methods shown in Nathan and McMahon (1990)
  #Eckhardt suggests BFI_max = 0.80 for perrennial, porous aquifer streams; 0.50 for ephemeral, porous;
  #0.25 for perennial with hard rock aquifer
  
  #It seems that Eckhardt (2005) suggests using only 1 pass. That'll be the default.
  suppressWarnings(Ends<-c(1,length(streamflow))*rep(1,(passes+1))) # Start and end values for the filter function
  suppressWarnings(AddToStart<-c(1,-1)*rep(1,passes))
  btP<-streamflow##Previous pass's baseflow approximation
  qft<-vector(length=length(streamflow))
  bt<-vector(length=length(streamflow))
  bt[1]<-if(streamflow[1]<quantile(streamflow,0.25)) streamflow[1] else mean(streamflow)/1.5
  ##Guess baseflow value in first time step.  
  for(j in 1:passes){
    for (i in (Ends[j]+AddToStart[j]):Ends[j+1]){
      bk <- ((((1-BFI_max)*filter_parameter*bt[i-AddToStart[j]])+((1-filter_parameter)*BFI_max*btP[i]))/(1-(filter_parameter*BFI_max)))
      
      if (bk > btP[i]){
        bt[i]<-btP[i]
      } else bt[i]<-bk
      qft[i]<-streamflow[i]-bt[i]
    }
    if (j<passes){
      btP<-bt
      bt[Ends[j+1]]<-if(streamflow[Ends[j+1]]<mean(btP))streamflow[Ends[j+1]]/1.2 else mean(btP)
      ##Refines the approximation of end values after the first pass
    }
  }
  f <- data.frame(bt,qft)
  return(f)
}
