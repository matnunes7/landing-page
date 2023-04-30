MainFunction <- function (SCENARIO, NF, D_SO, D_EB, D_LB, LR_SO, LR_LB, D_R,
                          D_BW, LR_R, LR_EB, T, C0_t, Q_t){
  for(case in 1:NROW(SCENARIO)){
    
    # Duration  of filter cycle 
    D <- (D_R[case]+D_SO[case]+D_EB[case]+D_LB[case]+D_BW[case])*3600
    
    # Initializing matrices of log removal and effluent concentration
    LR_ti <- matrix(, nrow=T, ncol=NF[case]) #instantaneous log removal of each filter
    status_ti <- matrix(, nrow=T, ncol=NF[case]) #instantaneous status (active or not) of each filter 
    phase_ti <- matrix(, nrow=T, ncol=NF[case]) #instantaneous phase of each filter
    LR_t <- c() # cumulative log removal of set of filters
    Ce_t <- c() # cumulative effluent concentration of set of filters
   
    
    # Startup condition for each filter (time within the filter cycle in s)
    
    # for evenly staggered 2,4,8,20,40 filters
    t0_1 <- c(1*2/D,1)*D/2
    t0_2 <- c(1*4/D,1:3)*D/4
    t0_3 <- c(1*8/D,1:7)*D/8
    t0_4 <- c(1*20/D,1:19)*D/20
    t0_5 <- c(1*50/D,1:49)*D/50
    #for backwash staggering scenario
    set.seed(5)
    t0_e <- c(1/3600,15.625,31.25,46.875)*3600 # a - even
    t0_ne1 <- runif(4,1,62.5*3600) # b - random
    t0_ne2 <- c(1/3600,0.5,1,1.5)*3600 # c - sequential backwash (no overlap, but overlap poor phases)
    t0_ne3 <- c(1/3600,10,20,30)*3600  # d - no overlap backwash and breakthrough
    t0_ne4 <- c(1/3600,0.25,0.50,0.75)*3600 # e - overlap backwash and poor
    if (identical(SCENARIO, nf)){
      t0 <- cbind(t0_1, t0_2, t0_3, t0_4, t0_5)}
    else if (identical(SCENARIO, stg)){
      t0 <- cbind(t0_e, t0_ne1,t0_ne2,t0_ne3,t0_ne4)}
    else {t0 <- cbind(t0_2,t0_2,t0_2,t0_2,t0_2)}
    
    # Determining status of each filter (active/inactive) and what is the log-removal at each time step
    
    for (i in 1:NF[case]){
      nt=1
      for (t in t0[i,case]:(t0[i,case]+D-1)){
        # instantaneous log removal - depends on the phase of filter cycle
        if ((t%%D) <= D_R[case]*3600){
          LR_ti[nt,i] <- LR_R[case]
          status_ti[nt,i] <- 1 # 1 for active, 0 for inactive
          phase_ti[nt,i] <- "R"} 
        else if ((t%%D) > D_R[case]*3600 & (t%%D) <= (D_R[case]+D_SO[case])*3600){
          LR_ti[nt,i] <- LR_SO[case]
          status_ti[nt,i] <- 1
          phase_ti[nt,i] <- "SO"} 
        else if ((t%%D) > (D_R[case]+D_SO[case])*3600 & (t%%D) <= (D_R[case]+D_SO[case]+D_EB[case])*3600){
          LR_ti[nt,i] <- LR_EB
          status_ti[nt,i] <- 1
          phase_ti[nt,i] <- "EB"} 
        else if ((t%%D) > (D_R[case]+D_SO[case]+D_EB[case])*3600 & (t%%D) <= (D_R[case]+D_SO[case]+D_EB[case]+D_LB[case])*3600){
          LR_ti[nt,i] <- LR_LB[case]
          status_ti[nt,i] <- 1
          phase_ti[nt,i] <- "LB"}
        else if ((t%%D) > (D_R[case]+D_SO[case]+D_EB[case]+D_LB[case])*3600 & (t%%D) <= (D_R[case]+D_SO[case]+D_EB[case]+D_LB[case]+D_BW[case])*3600){
          LR_ti[nt,i]= 0
          status_ti[nt,i] <- 0
          phase_ti[nt,i] <- "B"}
        nt<-nt+1
      }
    }
    
    # Only for backstaggering case, checking the max number of filters being backwashed at the same time
    
    if (identical(SCENARIO,stg)){
     if(min(rowSums(status_ti)) == 3){print("No overlapping")}
        else if(min(rowSums(status_ti)) == 2){print("2 filters in backwash")}
          else if(min(rowSums(status_ti)) == 1){print("3 filters in backwash")}
    }
    
    # Calculating combined and effective log removal
    
    active <- rowSums(status_ti)
    if (identical(SCENARIO,ftw)){status_ti[which(phase_ti=="R")]<-0}    
    Q_ti <- status_ti*(Q_t/NF[case])
    R_ti <- 10^(-LR_ti)
    Ce_t <- cumsum(C0_t*rowSums(Q_ti*R_ti))/cumsum(rowSums(Q_ti)) # effective effluent concentration
    LR_t <- -log10(Ce_t/C0_t)  # effective log removal
    Ce_t2 <- C0_t*rowSums(Q_ti*R_ti)/rowSums(Q_ti) # combined effluent concentration
    LR_t2 <- -log10(Ce_t2/C0_t)  # combined log removal
    R_t2 <- 10^(-LR_t2)
    
    # Plot
    
    t <- c(1:T)/3600
    if (identical(SCENARIO, nf)){
      plot(t, LR_t2, xlim=c(3600,T)/3600, ylim=c(1,5.5),type="l",
         col = c("blue"), cex.axis=1.5, yaxt='n', cex.main=1.5,
         main= if (case == 2){paste(SCENARIO[case], "(base case)")} else {paste(SCENARIO[case])})
      title("Number of filters", outer = TRUE, line = 1, cex.main=2)}
    else if (identical(SCENARIO, d_so)){
      plot(t, LR_t2, xlim=c(3600,D)/3600, ylim=c(1,5.5),type="l",
           col = c("blue"), cex.axis=1.5, yaxt='n', cex.main=1.5,
           main= if (case == 3){paste(SCENARIO[case], "(base case)")} else {paste(SCENARIO[case])})
      title("Duration of stable operation (h)", outer = TRUE, line = 1, cex.main=2)}
    else if (identical(SCENARIO, LR_so)){
      plot(t, LR_t2, xlim=c(3600,T)/3600, ylim=c(1,7),type="l",
           col = c("blue"), cex.axis=1.5, yaxt='n', cex.main=1.5,
           main= if (case == 3){paste(formatC(SCENARIO[case], format = "f", digits=1), "(base case)")} else {paste(formatC(SCENARIO[case], format = "f", digits=1), "")})
      title("Removal of stable operation", outer = TRUE, line = 1, cex.main=2)}
    else if (identical(SCENARIO, d_eb)){
      plot(t, LR_t2, xlim=c(3600,T)/3600, ylim=c(1,5.5),type="l",
           col = c("blue"), cex.axis=1.5, yaxt='n', cex.main=1.5,
           main= if (case == 1){paste(SCENARIO[case], " / ",d_lb[case], " (base case)")} else {paste(SCENARIO[case], " / ",d_lb[case])})
      title("Duration of early / late breakthrough (h)", outer = TRUE, line = 1, cex.main=2)}
    else if (identical(SCENARIO, LR_lb)){
      plot(t, LR_t2, xlim=c(3600,T)/3600, ylim=c(0,5.5),type="l",
           col = c("blue"), cex.axis=1.5, yaxt='n', cex.main=1.5,
           main= if (case == 3){paste(formatC(SCENARIO[case], format = "f", digits=1), "(base case)")} else {paste(formatC(SCENARIO[case], format = "f", digits=1), "")})
      title("Removal of late breakthrough", outer = TRUE, line = 1, cex.main=2)}
    else if (identical(SCENARIO, stg)){
      plot(t, LR_t2, xlim=c(3600,T)/3600, ylim=c(1,5.5),type="l",
           col = c("blue"), cex.axis=1.5, yaxt='n', cex.main=1.5,
           main= if (case == 1){paste(SCENARIO[case], "(base case)")} else {paste(SCENARIO[case],"")})
      title("Backwash staggering", outer = TRUE, line = 1, cex.main=2)}
    else if (identical(SCENARIO, LR_so_coag)){
      plot(t, LR_t2, xlim=c(3600,T)/3600, ylim=c(1,5.5),type="l",
           col = c("blue"), cex.axis=1.5, yaxt='n', cex.main=1.5,
           main= if(case==1){c("No coagulation")}else if (case==2){c("Sub-optimal coagulation")}else{c("Optimal coagulation (base case)")})}
    else if (identical(SCENARIO, ftw)){
      plot(t, LR_t2, xlim=c(3600,T)/3600, ylim=c(1,5.5),type="l",
           col = c("blue"), cex.axis=1.5, yaxt='n', cex.main=1.5,
           main= if (case == 1){paste(SCENARIO[case], "")} else {paste(SCENARIO[case],"")})
      title("Removal of ripening", outer = TRUE, line = 1, cex.main=2)}
    else if (identical(SCENARIO, LR_r)){
      plot(t, LR_t2, xlim=c(3600,T)/3600, ylim=c(1,5.5),type="l",
           col = c("blue"), cex.axis=1.5, yaxt='n', cex.main=1.5,
           main= if (case == 1){paste(formatC(SCENARIO[case], format = "f", digits=1), "(base case)")} else {paste(formatC(SCENARIO[case], format = "f", digits=1), "")})}
    
    title(xlab="Time (h)", ylab="Pathogen log-reduction",cex.lab=1.5, outer = TRUE, line = 3)
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = 3)
    if (identical(SCENARIO, ftw)|identical(SCENARIO, LR_r)){
      xx <- c(D/25,0.35*D,0.35*D,D/25)/3600
      yy <- c(LR_t[D]+0.05,LR_t[D]+0.05,LR_t[D]+0.3,LR_t[D]+0.3)
      polygon(xx, yy, col="white", border=FALSE)
      text(x=D/(3600*5), y=LR_t[D]+0.175, labels =paste(formatC(LR_t[D], format = "f", digits = 3)),
         las=2, mgp=c(3,0.35,0), cex=1.5, col="darkorange")}
    else{
      xx <- c(D/25,D/4,D/4,D/25)/3600
      yy <- c(LR_t[D]+0.05,LR_t[D]+0.05,LR_t[D]+0.3,LR_t[D]+0.3)
      polygon(xx, yy, col="white", border=FALSE)
      text(x=3*D/(3600*20), y=LR_t[D]+0.175, labels =paste(formatC(LR_t[D], format = "f", digits = 2)),
           las=2, mgp=c(3,0.35,0), cex=1.5, col="darkorange")}
    
    if (identical(SCENARIO, LR_so)){
      axis(side=2, at=1:7,labels = if (case == 1) {1:7} else FALSE, cex.axis=1.5)}
    else if (identical(SCENARIO, LR_r)){}
    else {axis(side=2, at=1:5.5,labels = if (case == 1) {1:5.5} else FALSE, cex.axis=1.5)}
    
    abline(h=LR_t[D], col="darkorange", lty=1, lwd=2)
