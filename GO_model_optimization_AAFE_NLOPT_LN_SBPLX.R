#Graphene oxide model

library(deSolve)
source("/Users/eviepapakyriakopoulou/Documents/GitHub/Graphene_oxide_Nimble_model/Final model/Goodness-of-fit-metrics.R")
setwd("/Users/eviepapakyriakopoulou/Documents/GitHub/Graphene_oxide_Nimble_model/Final model")


create.params <- function(user.input){
  with(as.list(user.input),{
    
    
    #Volumes of organs as percent of BW
    #blood, kidneys, liver, stomach, Small intestine, large intestine, lungs, 
    #spleen, heart, brain, RoB  
    MW <- 124.91 #g/mol
    
    #Blood
    PVB <- 1.7e-3/0.02 #Davies et al. 1993, 1.7 for BW = 0.02 kg
    VB <- PVB * BW #blood volume kg=L
    PVplasma <- 1e-3 /0.02 #Davies et al. 1993, 1.0 for BW = 0.02 kg
    Vplasma <- PVplasma * BW #plasma volume kg=L
    VBven <- BW*0.5/20 	#volume of venous plasma (L); from doi:https://doi.org/10.1007/bf02353860
    VBart <- BW*0.22/20	#volume of arterial plasma (L); from doi:https://doi.org/10.1007/bf02353860
    
    #Kidney
    PVKi <- 1.67e-2 #Brown et al. 1997, Table 4
    VKi <- PVKi * BW #kidney volume kg=L
    
    #Liver
    PVLi <- 5.49e-2 #Brown et al. 1997, Table 4
    VLi <- PVLi * BW #liver volume kg=L
    V_macro_Li = 27.5/100*VLi # https://doi.org/10.3892/etm.2019.7450
    
    #Stomach
    PVSt <- 6e-3 #Brown et al. 1997, Table 4
    VSt <- PVSt * BW #Stomach volume kg=L
    VSt_fluid <- 0.37 #mL https://doi.org/10.1211/jpp.60.1.0008 TO CHECK
    
    #Small intestine
    PVSIn <- 2.53e-2 #Brown et al. 1997, Table 4
    VSIn <- PVSIn * BW #Small intestine volume kg=L
    
    #Large intestine
    PVLIn <- 1.09e-2 #Brown et al. 1997, Table 4
    VLIn <- PVLIn * BW #Large intestine volume kg=L
    
    #Lung
    PVLn <- 7.3e-3 #Brown et al. 1997, Table 4
    VLn <- PVLn * BW #Lung volume kg=L
    V_macro_Ln <- 7.5/100 * VLn #random
    
    #Spleen
    PVSpl <- 3.5e-3 #Brown et al. 1997, Table 4
    VSpl <- PVSpl * BW #Spleen volume kg=L
    V_macro_Spl = 6.94/100*VSpl #https://doi.org/10.1038/s41374-018-0137-1
    
    #Heart
    PVH <- 5e-3 #Brown et al. 1997, Table 4
    VH <- PVH * BW #Heart volume kg=L
    
    #Brain
    PVBr <- 1.65e-2 #Brown et al. 1997, Table 4
    VBr <- PVBr * BW #Brain volume kg=L
    
    # #Thyroid glands
    # PVTg <- 4.92e-3/24.2 #g tisse/g BW, 4.92 uL = 4.92 mg, Mancini et al., 2009 https://doi.org/10.1210/en.2009-0417
    # VTg <- PVTg * BW #Thyroid glands volume kg=L
    
    #RoB
    PVRe <- 1 - PVB - PVKi - PVLi - PVSt - PVSIn - PVLIn - PVLn - PVSpl - PVH - PVBr #- PVTg
    VRe <- PVRe * BW #volume of the rest of the body kg=L
    
    #Capillary surface area for each tissue (Ai) as percentage of body weight (m^2/kg), values from pkSim "Endothelial Surface area"
    
    PAKi <- 33.92e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    AKi <- PAKi * BW #the surface area of kidney (m^2)
    
    PALi <- 142.0e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ALi <- PALi * BW #liver surface area (m^2)
    
    PASt <- 3.34e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ASt <- PASt * VSt * 1e3 #stomach surface area (m^2)
    
    PASIn <- 9.62e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ASIn <- PASIn * VLIn * 1e3 #Small intestine surface area (m^2)
    
    PALIn <- 5.88e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ALIn <- PALIn * VLIn * 1e3 #Small intestine surface area (m^2)
    
    PLn <- 59.47e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ALn <- PLn* BW #lung surface area (m^2)
    
    PSpl <- 26.79e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ASpl <- PSpl* BW #spleen surface area (m^2), same as muscle #assumption
    
    PH <-  23.65e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    AH <- PH* VH #heart surface area (m^2)
    
    PBr <- 5.98e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ABr <- PBr* VBr *1e03 #brain surface area (m^2)
    
    ARe <- (AKi+ALi+ASt+ASIn+ALIn+ALn+ASpl+AH+ABr)/9 #assumption ???
    
    
    np_size_small <- 148 #nm,  #3/2
    np_size_large <- 556 #nm, #300/2
    
    
    #(QBi, in L/min) to different tissues (i=L, K, G, A, M, R)
    #Hall et al., 2012, https://doi.org/10.1002/jps.22811
    #PQi, fraction of cardiac output
    
    BW_ref <- 0.03 #kg
    Cardiac_output <- 11.4/1000 * BW/BW_ref #mL/min --> L/min
    
    QKi <- 1.3/1000 * BW/BW_ref # mL/min --> L/min
    PQKi <- QKi/Cardiac_output
    QLi <- 1.8/1000 * BW/BW_ref # mL/min --> L/min
    PQLi <- QLi/Cardiac_output
    QLn <- 1
    QSpl <- 0.09/1000 * BW/BW_ref # mL/min --> L/min
    PQSpl <- QSpl/Cardiac_output
    QH <- 0.28/1000 * BW/BW_ref # mL/min --> L/min
    PQH <- QH/Cardiac_output
    QBr <- 0.26/1000 * BW/BW_ref # mL/min --> L/min
    PQBr <- QBr/Cardiac_output
    PQSt <-((0.53+1.45)/2)/100 #Table 1 (mean of fore and glandular), https://doi.org/10.1002/jat.2550030607
    QSt <- PQSt*Cardiac_output 
    PSIn <- 13.5/100 #Table 1, https://doi.org/10.1002/jat.2550030607
    QSIn <- PSIn*Cardiac_output 
    PLIn <- 3.67/100 #Table 1, https://doi.org/10.1002/jat.2550030607
    QLIn <- PLIn*Cardiac_output 
    PQRe <- 1 - PQKi - PQLi - PQSpl - PQH - PQBr - PQSt - PSIn - PLIn
    QRe <- PQRe*Cardiac_output
    
    
    Qtotal <- QKi+QLi+QRe+QSpl+QH+QBr+QSIn+QSt+QLIn
    QLitot <- QLi+QSpl+QSIn+QSt+QLIn
    QGE <- 1.233/60*BW^0.25 #gastric emptying time (1/(min*BW^0.25)); https://doi.org/10.1002/mrm.10207
    
    if(sex == "M"){
      Vper <- 0.02*1e-3 #L peritoneal fluid volume https://doi.org/10.1038/2101123a0
      
    }else{
      Vper <- 0.1*1e-3 #L peritoneal fluid volume https://doi.org/10.1038/2101123a0
    }
    
    
    return(list('VB'=VB,'Vplasma'=Vplasma,'VBven'=VBven,'VBart'=VBart,
                'VKi'=VKi,'VLi'=VLi, 'V_macro_Li'=V_macro_Li,'VSt'=VSt, 'VSt_fluid'=VSt_fluid,
                'VSIn'=VSIn,'VLIn'=VLIn,'VLn'=VLn, 'V_macro_Ln'=V_macro_Ln,
                'VSpl'=VSpl,'VH'=VH,'VBr'=VBr,'VRe'=VRe,
                'V_macro_Spl'=V_macro_Spl,
                'AKi'=AKi,'ALi'=ALi,'ASt'=ASt,'ASIn'=ASIn,'ALIn'=ALIn,
                'ALn'=ALn,'ASpl'=ASpl,'AH'=AH,'ABr'=ABr,'ARe'=ARe,
                
                'QKi'=QKi, 'QLi'=QLi, 'QRe'=QRe, 'QLn'=QLn, 'QSpl'=QSpl, 'QH'=QH,
                'QBr'=QBr, 'QSt'=QSt,'QSIn'=QSIn, 'QLIn'=QLIn,
                
                "admin.time" = admin.time, "admin.dose" = admin.dose,
                "admin.type" = admin.type, "MW"=MW, 
                "np_size_small"=np_size_small, "np_size_large"=np_size_large,
                "Qtotal"=Qtotal, "QLitot"=QLitot, "QGE"=QGE,
                "sex"=sex, "estimated_params"=estimated_params,
                'Vper'=Vper
                
                
    ))
  })
}  


ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
#coef_i is the ratio of permeability coefficient (xi) between capillary blood and organ
#and partition coefficient of nanoparticles between tissue and blood (Pi), coef_i=xi/Pi
# based on the work of Li et al., 2013, https://doi.org/10.3109/17435390.2013.863406
  
      
    coef_liver <- estimated_params[1]
    coef_spleen <- estimated_params[2]
    coef_kidney <- estimated_params[3]
    coef_heart <- estimated_params[4]
    coef_lung <- estimated_params[5]
    coef_brain <- estimated_params[6]
    coef_rob <- estimated_params[7]
    coef_stomach <- estimated_params[8]
    coef_smallIn <- estimated_params[9]
    coef_largeIn <- estimated_params[10]
    CLurine <- estimated_params[11]
    CLfeces <- estimated_params[12]
    kabs_oral <- estimated_params[13]
    kabs_ip <- estimated_params[14]
    
    
    # Blood concentration
    CBven <- MBven/VBven
    CBart <- MBart/VBart
    
    # Kidney 
    CKi = MKi/VKi # tissue concentration
    
    #Liver
    CLi = MLi/VLi # tissue concentration
    
    #Stomach
    CSt = MSt/VSt # tissue concentration
    CStL = MStL/VSt_fluid
    
    #Small Intestine
    CSIn = MSIn/VSIn # tissue concentration
    
    #Large Intestine
    CLIn = MLIn/VLIn # tissue concentration
    
    #Lungs
    CLn = MLn/VLn # tissue concentration
    
    #Spleen
    CSpl = MSpl/VSpl # tissue concentration
    
    #Heart
    CH = MH/VH # tissue concentration
    
    #Brain
    CBr = MBr/VBr # tissue concentration
    
    #Rest-of-the-body
    CRe = MRe/VRe # tissue concentration
    

    
    #========================================================================================================
    #Arterial Blood
    dMBart =  QLn*(CLn/coef_lung) -QKi*CBart -QLi*CBart -QSt*CBart -QSIn*CBart -QLIn*CBart -
      QSpl*CBart -QH*CBart -QBr*CBart -QRe*CBart
    
    #Venous Blood
    dMBven = -QLn*CBven +QKi*(CKi/coef_kidney) +QLitot*(CLi/coef_liver) +
      QH*(CH/coef_heart) +QBr*(CBr/coef_brain) +QRe*(CRe/coef_rob) + kabs_ip*Mper
    
    #Kidney
    dMKi = QKi*(CBart - (CKi/coef_kidney)) - CLurine*MKi
    
    #Liver
    dMLi =  QLi*CBart - QLitot*(CLi/coef_liver) + QSt*(CSt/coef_stomach) + QSpl*(CSpl/coef_spleen) + 
      QSIn*(CSIn/coef_smallIn) + QLIn*(CLIn/coef_largeIn) 
    
    #Stomach
    dMSt = kabs_oral*CStL + QSt*(CBart - (CSt/coef_stomach))
    #Stomach lumen
    dMStL = - QGE*CStL - kabs_oral*CStL 
    
    #Small Intestine
    dMSIn = QGE*CStL + QSIn*(CBart - (CSIn/coef_smallIn)) 
    
    #Large Intestine
    dMLIn = QLIn*(CBart - (CLIn/coef_largeIn)) - CLfeces*CLIn
    
    #Lung 
    dMLn = QLn*(CBven - (CLn/coef_lung))
    
    #Spleen
    dMSpl = QSpl*(CBart - (CSpl/coef_spleen)) 
    
    #Heart
    dMH = QH*(CBart - (CH/coef_heart))
    
    #Brain
    dMBr = QBr*(CBart - (CBr/coef_brain))
    
    #Rest of body
    dMRe = QRe*(CBart - (CRe/coef_rob))
    
    # Urine
    dMurine = CLurine*MKi
    # Feces
    dMfeces = CLfeces*CLIn
    
    #Peritoneal space 
    dMper = -kabs_ip*Mper
    
    dVurine = CLurine
    dVfeces = CLfeces
    
    
    
    #Concentration calculation in each compartment 
    
    Cblood <- (MBven + MBart)/ (VBven + VBart)
    Ckidneys <- MKi/VKi
    Cliver <- MLi /VLi
    Cstomach <-  MSt/VSt
    Csmall_intestine <-  MSIn/VSIn
    Clarge_intestine <- MLIn/VLIn
    Clungs <-  MLn /VLn
    Crest <-  MRe/VRe
    Cspleen <- MSpl /VSpl
    Cheart <- MH/VH
    Cbrain <-  MBr/VBr
    
    Mblood <- (MBven + MBart)
    Mkidneys <- MKi
    Mliver <- MLi
    Mstomach <-  MSt
    Msmall_intestine <-  MSIn
    Mlarge_intestine <- MLIn
    Mintestine <- MSIn+MLIn
    Mlungs <-  MLn
    Mrest <-  MRe
    Mspleen <- MSpl
    Mheart <- MH
    Mbrain <-  MBr
    
    
    
    list(c( 
      "dMBart"=dMBart, "dMBven"=dMBven,"dMKi"=dMKi,
      "dMLi"=dMLi, "dMSt"=dMSt, "dMStL"=dMStL,
      "dMSIn"=dMSIn, "dMLIn"=dMLIn, 
      "dMLn"=dMLn, "dMSpl"=dMSpl, "dMH"=dMH, "dMBr"=dMBr,
      "dMRe"=dMRe, "dMurine"=dMurine, 
      "dMfeces"=dMfeces,  "dMper"=dMper,
      "dVurine"=dVurine, "dVfeces"=dVfeces), 
      
      'Cblood'=Cblood, 'Ckidneys'=Ckidneys, 'Cliver'=Cliver, 'Cstomach'=Cstomach,
      'Csmall_intestine'=Csmall_intestine, 'Clarge_intestine'=Clarge_intestine,
      'Clungs'=Clungs, 'Crest'=Crest,
      'Cspleen'=Cspleen, 'Cheart'=Cheart, 'Cbrain'=Cbrain,
      
      'Mblood'=Mblood, 'Mkidneys'=Mkidneys, 'Mliver'=Mliver, 'Mstomach'=Mstomach,
      'Msmall_intestine'=Msmall_intestine, 'Mlarge_intestine'=Mlarge_intestine,
      'Mlungs'=Mlungs, 'Mrest'=Mrest, 
      'Mspleen'=Mspleen, 'Mheart'=Mheart, 'Mbrain'=Mbrain, 'Mintestine'=Mintestine,
      
      'CBven'=CBven, 'CBart'=CBart,'CKi'=CKi, 
      'CLi'=CLi, 'CSt'=CSt, 'CSIn'=CSIn, 
      'CLIn'=CLIn, 'CLn'=CLn, 'CSpl'=CSpl, 
      'CH'=CH, 'CBr'=CBr,'CRe'=CRe)
            
  
  })
}


create.inits <- function(parameters){
  with(as.list(parameters),{
    
    MBart<-0; MBven<-0;MKi<-0; MLi<-0; MSt<-0; MStL <-0;
    MSIn<-0; MLIn<-0; MLn<-0; MSpl<-0; MH<-0; MBr<-0;
    MRe<-0; Murine<-0; Mfeces<-0; Mper<-0;
    Vurine <-0; Vfeces <-0
    
    return(c(
      "MBart"=MBart, "MBven"=MBven,"MKi"=MKi,
      "MLi"=MLi, "MSt"=MSt,"MStL"=MStL, "MSIn"=MSIn,
      "MLIn"=MLIn,"MLn"=MLn,"MSpl"=MSpl, 
      "MH"=MH, "MBr"=MBr,"MRe"=MRe,
      "Murine"=Murine, "Mfeces"=Mfeces,"Mper"=Mper, 
      "Vurine"=Vurine, "Vfeces"=Vfeces
    ))
  })
}


create.events <- function(parameters){
  with(as.list(parameters), {
    
    # Calculate number of administrated doses and corresponding administration time
    ldose <- length(admin.dose)
    ltimes <- length(admin.time)
    # If not equal, then stop 
    if (ltimes != ldose){
      stop("The times of administration should be equal in number to the doses")
    }else{
      if (admin.type == "iv"){
        events <- list(data = rbind(data.frame(var = c("MBven"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }else if (admin.type == "oral"){
        events <- list(data = rbind(data.frame(var = c("MStL"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }else if (admin.type == "inh"){
        events <- list(data = rbind(data.frame(var = c("MLn"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }else if (admin.type == "ip"){
        events <- list(data = rbind(data.frame(var = c("Mper"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }
    }
    return(events)
  })
}

obj.func <- function(x, dataset){
  N_data <- length(dataset)
  score <- rep(NA, N_data)
  
  # x: a vector with the values of the optimized parameters (it is not the x
  # from the odes!!!)
  estimated_params <- exp(x)

  ##########################
  #-------------------------
  # Liu et al., 2012
  #-------------------------
  ##########################
  # Set up simulations for the 1st case, i.e. Liu (2012) 1 mg/kg small_tissues
  
  BW <- 0.03  # body weight (kg) #not reported
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- c(admin.dose_per_kg*BW*1e03) #ug PFOA
  admin.time <- c(0) # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,10,0.05) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df1=========================================================
  
  exp_data <- dataset$df1 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  preds_Liu_1_small_tissues <- as.data.frame(solution[solution$time %in% unique(exp_data$time), c("Mheart",
                                                                                                  "Mliver","Mspleen", "Mstomach",
                                                                                                  "Mkidneys", "Mlungs", "Mbrain",
                                                                                                  "Msmall_intestine", "Mlarge_intestine")])
  
  
  
  obs_Liu_1_small_tissues <- list(exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Small_intestine", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Large_intestine", "mass"]*admin.dose/100)
  
  score[1] <- AAFE(predictions = preds_Liu_1_small_tissues, observations = obs_Liu_1_small_tissues)
  
  
  # Set up simulations for the 2nd case, Liu (2012) 1 mg/kg small_tissues_diff_time_points
  
  BW <- 0.03  # body weight (kg) #not reported
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 3/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,180,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df2=========================================================
  
  exp_data <- dataset$df2 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg diff_time_points
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mheart", "Mliver","Mspleen", "Mstomach","Mkidneys", "Mlungs", "Mbrain",
                    "Msmall_intestine", "Mlarge_intestine")
  preds_Liu_1_small_diftp_tissues <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Liu_1_small_diftp_tissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  obs_Liu_1_small_diftp_tissues <- list( exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Small_intestine", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Large_intestine", "mass"]*admin.dose/100) 
  
  
  score[2] <- AAFE(predictions = preds_Liu_1_small_diftp_tissues, observations = obs_Liu_1_small_diftp_tissues)
  
  
  
  # Set up simulations for the 3rd case, i.e. Liu (2012) 1 mg/kg large_tissues_diff_time_points
  
  
  BW <- 0.03  # body weight (kg) #not reported
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_large <- 300/2
  np_size <- np_size_large #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,180,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df3=========================================================
  
  exp_data <- dataset$df3 # retrieve data of Liu et al. 2012 tissues large p.s., 1 mg/kg diff_time_points
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mheart", "Mliver","Mspleen", "Mstomach","Mkidneys", "Mlungs", "Mbrain",
                    "Msmall_intestine", "Mlarge_intestine")
  preds_Liu_1_large_diftp_tissues <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Liu_1_large_diftp_tissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  obs_Liu_1_large_diftp_tissues <- list( exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Small_intestine", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Large_intestine", "mass"]*admin.dose/100) 
  
  score[3] <- AAFE(predictions = preds_Liu_1_large_diftp_tissues, observations = obs_Liu_1_large_diftp_tissues)
  
  
  # Set up simulations for the 4th case, i.e. Liu (2012) 2 mg/kg small_tissues
  
  BW <- 0.03  # body weight (kg) #not reported
  admin.dose_per_kg <- 2 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 3/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,10,0.05) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  
  #======================================df4=========================================================
  
  exp_data <- dataset$df4 # retrieve data of Liu et al. 2012 tissues Small p.s., 2 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_Liu_2_small_tissues <- as.data.frame(solution[solution$time %in% unique(exp_data$time), c("Mheart",
                                                                                                  "Mliver","Mspleen", "Mstomach",
                                                                                                  "Mkidneys", "Mlungs", "Mbrain",
                                                                                                  "Msmall_intestine", "Mlarge_intestine")])
  
  
  
  obs_Liu_2_small_tissues <- list(exp_data[exp_data$Tissue == "Heart", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Liver", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Spleen", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Stomach", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Kidneys", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Lungs", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Brain", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Small_intestine", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Large_intestine", "concentration"]*admin.dose/100)
  
  score[4] <- AAFE(predictions = preds_Liu_2_small_tissues, observations = obs_Liu_2_small_tissues)
  
  
  # Set up simulations for the 5th case, Liu (2012) 10 mg/kg small_tissues
  
  BW <- 0.03  # body weight (kg) #not reported
  admin.dose_per_kg <- 10 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 3/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,10,0.05) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df5=========================================================
  
  exp_data <- dataset$df5 # retrieve data of Liu et al. 2012 tissues Small p.s., 10 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_Liu_10_small_tissues <- as.data.frame(solution[solution$time %in% unique(exp_data$time), c("Mheart",
                                                                                                   "Mliver","Mspleen", "Mstomach",
                                                                                                   "Mkidneys", "Mlungs", "Mbrain",
                                                                                                   "Msmall_intestine", "Mlarge_intestine")])
  
  
  
  obs_Liu_10_small_tissues <- list(exp_data[exp_data$Tissue == "Heart", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Liver", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Spleen", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Stomach", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Kidneys", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Lungs", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Brain", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Small_intestine", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Large_intestine", "concentration"]*admin.dose/100)
  
  score[5] <- AAFE(predictions = preds_Liu_10_small_tissues, observations = obs_Liu_10_small_tissues)
  
  
  # Set up simulations for the 6th case, i.e. Liu (2012) 1 mg/kg small_blood diff_time_points
  
  BW <- 0.03  # body weight (kg) #not reported
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 3/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,180,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df6=========================================================
  
  exp_data <- dataset$df6 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mblood")
  
  preds_Liu_1_small_diftp_blood <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Liu_1_small_diftp_blood [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  
  obs_Liu_1_small_diftp_blood <- list(exp_data[exp_data$Tissue == "Blood", "mass"]*admin.dose/100)
  
  score[6] <- AAFE(predictions = preds_Liu_1_small_diftp_blood, observations = obs_Liu_1_small_diftp_blood)
  
  
  # Set up simulations for the 7th case, i.e. Liu (2012) 1 mg/kg large_blood_diff_time_points
  
  
  BW <- 0.03  # body weight (kg) #not reported
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_large <- 300/2
  np_size <- np_size_large #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M"
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,180,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df7=========================================================
  
  exp_data <- dataset$df7 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mblood")
  
  preds_Liu_1_large_diftp_blood <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Liu_1_large_diftp_blood [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  
  obs_Liu_1_large_diftp_blood <- list(exp_data[exp_data$Tissue == "Blood", "mass"]*admin.dose/100)
  
  score[7] <- AAFE(predictions = preds_Liu_1_large_diftp_blood, observations = obs_Liu_1_large_diftp_blood)
  
  
  
  ##########################
  #-------------------------
  # Li et al., 2013
  #-------------------------
  ##########################
  
  # Set up simulations for the 8th case, i.e. Li (2013) 10 mg/kg nanoGO_tissues_diff_time_points
  
  BW <- 0.03  # body weight (kg) #not reported
  admin.dose_per_kg <- 10 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "inh"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,720,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df8=========================================================
  
  exp_data <- dataset$df8 # retrieve data of Li (2013) 10 mg/kg nanoGO_tissues_diff_time_points
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mblood", "Mheart", "Mlungs","Mliver", "Mspleen", "Mkidneys",
                    "Mstomach","Msmall_intestine", "Mlarge_intestine")
  preds_Li_IN_tissues <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Li_IN_tissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  obs_Li_IN_tissues <- list( exp_data[exp_data$Tissue == "Blood", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Kidney", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Small_intestine", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Large_intestine", "mass"]*admin.dose/100) 
  
  score[8] <- AAFE(predictions = preds_Li_IN_tissues, observations = obs_Li_IN_tissues)
  
  
  ##########################
  #-------------------------
  # Li et al., 2014
  #-------------------------
  ##########################
  
  # Set up simulations for the 9th case, i.e. Li (2014) 5 mg/kg nanoGO_tissues_diff_time_points
  
  BW <- 0.02717  # body weight (kg) 
  PVB <- 1.7e-3/0.02
  VB <- PVB * BW
  admin.dose_per_kg <- 5 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     'VB'=VB,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,360,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df9=========================================================
  
  exp_data <- dataset$df9 # retrieve data of Li (2014) 5 mg/kg nanoGO_tissues_diff_time_points
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mblood", "Mheart", "Mlungs","Mliver", "Mspleen", "Mkidneys",
                    "Mstomach","Msmall_intestine", "Mlarge_intestine")
  preds_Li_IV_tissues <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Li_IV_tissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  obs_Li_IV_tissues <- list( exp_data[exp_data$Tissue == "Blood", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Kidney", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Small_intestine", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Large_intestine", "mass"]*admin.dose/100) 
  
  score[9] <- AAFE(predictions = preds_Li_IV_tissues, observations = obs_Li_IV_tissues)
  
  #======================================df10=========================================================
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,1440,5) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df10=========================================================
  
  exp_data <- dataset$df10 # retrieve data of Li (2014) 5 mg/kg nanoGO_tissues_diff_time_points
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mblood")
  preds_Li_IV_blood <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Li_IV_blood[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  obs_Li_IV_blood <- list( exp_data[exp_data$Tissue == "Blood", "mass"]*admin.dose*VB*1000/100) 
  
  score[10] <- AAFE(predictions = preds_Li_IV_blood, observations = obs_Li_IV_blood)
  
  
  # Set up simulations for the 11th case, i.e. Li (2014) 5 mg/kg nanoPEGGO_tissues_diff_time_points
  
  BW <- 0.02874 # body weight (kg) 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  user_input <- list('BW'=BW,
                     'VB'=VB,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,360,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df11=========================================================
  
  exp_data <- dataset$df11 # retrieve data of Li (2014) 5 mg/kg nanoGO_tissues_diff_time_points
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mlungs","Mliver", "Mspleen")
  preds_Li_IV_tissues_PEG <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Li_IV_tissues_PEG[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  obs_Li_IV_tissues_PEG <- list( exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose/100) 
  
  score[11] <- AAFE(predictions = preds_Li_IV_tissues_PEG, observations = obs_Li_IV_tissues_PEG)
  
  
  #======================================df12=========================================================
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,1440,5) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df12=========================================================
  
  exp_data <- dataset$df12 # retrieve data of Li (2014) 5 mg/kg nanoGO_PEG_blood_diff_time_points
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mblood")
  preds_Li_IV_blood_PEG <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Li_IV_blood_PEG[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  obs_Li_IV_blood_PEG <- list( exp_data[exp_data$Tissue == "Blood", "mass"]*admin.dose*VB*1000/100) 
  
  score[12] <- AAFE(predictions = preds_Li_IV_blood_PEG, observations = obs_Li_IV_blood_PEG)
  
  
  ##########################
  #-------------------------
  # Mao et al., 2016
  #-------------------------
  ##########################

  # # Observed body weight data (Graphene group) in grams, from Figure S4
  # bw_days <- c(0, 3, 7, 10, 14, 21, 28)
  # bw_values_g <- c(18.5, 22, 25.5, 27, 31.5, 36, 34.5)
  # 
  # bw_time_min <- bw_days * 1440  # 1 day = 1440 min
  # time_min <- seq(0, 40320, by = 1440)
  # 
  # # Interpolate body weight at each time point, convert to kilograms
  # BW_interp_g <- approx(x = bw_time_min, y = bw_values_g, xout = time_min, rule = 2)$y
  # BW_interp_kg <- BW_interp_g / 1000  # convert to kg
  # 
  # # Now use your original structure to simulate BW initialization and tracking
  # BW <- BW_interp_kg[1]  # initial BW at time 0 (in kg)
  # 
  # BW_init <- BW  # store the initialized BW if needed
  # 
  # # Build the full BW vector (kg), labeled by time in minutes
  # BW <- c(BW_init, rep(NA, length(time_min) - 1))
  # for (i in 2:length(BW)) {
  #   BW[i] <- BW_interp_kg[i]
  # }
  # 
  
  BW_init <- 0.03
  admin.dose <- 5 #ug PFOA
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "inh"
  sex <- "M"

  user_input <- list('BW'=BW_init,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time,
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)

  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)

  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,40320,60) #min

  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))

  # We need to keep only the predictions for the relevant compartments for the time points
  # at which we have available data.

  #======================================df13=========================================================

  exp_data <- dataset$df13 # retrieve data of Ma0 (2016) 5 ug tissues
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mliver", "Mspleen", "Mstomach","Msmall_intestine",
                    "Mlarge_intestine","Mlungs", "Mfeces")
  preds_Mao_IN_tissues <- list()


  ## loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]

    preds_Mao_IN_tissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }



  obs_Mao_IN_tissues <- list(exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Small_intestine", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Large_intestine", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose/100,
                             exp_data[exp_data$Tissue == "Feces", "mass"]*admin.dose/100)

  score[13] <- AAFE(predictions = preds_Mao_IN_tissues, observations = obs_Mao_IN_tissues)


  #======================================df14=========================================================

  BW <- 0.03  # body weight (kg) #not reported
  admin.dose <- 10 #ug PFOA
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "oral"
  sex <- "M"

  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time,
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)

  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)

  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,4320,60) #min

  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))

  #======================================df14=========================================================

  exp_data <- dataset$df14 # retrieve data of Mao et al. 2016 tissues 10 ug
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mstomach","Msmall_intestine","Mlarge_intestine", "Mfeces")
  preds_Mao_OR_tissues <-  list()
  
  
  ## loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]

    preds_Mao_OR_tissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }



  obs_Mao_OR_tissues <- list(exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Small_intestine", "mass"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Large_intestine", "mass"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Feces", "mass"]*admin.dose/100)


  
  score[14] <- AAFE(predictions = preds_Mao_OR_tissues, observations = obs_Mao_OR_tissues)
  
  
  ##########################
  #-------------------------
  # Yang et al., 2010
  #-------------------------
  ##########################
  
  # Set up simulations for the 15th case, i.e. Yang (2010) 10 mg/kg nanoGO_tissues_diff_time_points
  
  BW <- 0.022  # body weight (kg) #mean value of the BW in figure S6
  PVB <- 1.7e-3/0.02
  VB <- PVB * BW
  admin.dose_per_kg <- 4 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "F" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,1440,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df15=========================================================
  
  exp_data <- dataset$df15 # retrieve data of Yang (2010) 4 mg/kg nanographene sheets blood
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mblood")
  preds_Yang_IV_blood <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Yang_IV_blood[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  obs_Yang_IV_blood <- list( exp_data[exp_data$Tissue == "Blood", "mass"]*admin.dose*VB*1000/100) 
  
  score[15] <- AAFE(predictions = preds_Yang_IV_blood, observations = obs_Yang_IV_blood)
  
  
  #======================================df16=========================================================
  PVLi <- 5.49e-2 #Brown et al. 1997, Table 4
  VLi <- PVLi * BW #liver volume kg=L
  PVSpl <- 3.5e-3 #Brown et al. 1997, Table 4
  VSpl <- PVSpl * BW #Spleen volume kg=L
  V_macro_Spl = 6.94/100*VSpl #https://doi.org/10.1038/s41374-018-0137-1
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,86400,60) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df16=========================================================
  
  exp_data <- dataset$df16 # retrieve data of Yang (2010) 4 mg/kg nanographene sheets liver, spleen
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mliver","Mspleen")
  preds_Yang_IV_LiSpl <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Yang_IV_LiSpl[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  obs_Yang_IV_LiSpl <- list( exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose*VLi*1000/100,
                             exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose*VSpl*1000/100)
  
  
  score[16] <- AAFE(predictions = preds_Yang_IV_LiSpl, observations = obs_Yang_IV_LiSpl)
  
  
  #======================================df17=========================================================
  PVLi <- 5.49e-2 #Brown et al. 1997, Table 4
  VLi <- PVLi * BW #liver volume kg=L
  PVSpl <- 3.5e-3 #Brown et al. 1997, Table 4
  VSpl <- PVSpl * BW #Spleen volume kg=L
  V_macro_Spl = 6.94/100*VSpl #https://doi.org/10.1038/s41374-018-0137-1
  PVKi <- 1.67e-2 #Brown et al. 1997, Table 4
  VKi <- PVKi * BW #kidney volume kg=L
  PVH <- 5e-3 #Brown et al. 1997, Table 4
  VH <- PVH * BW #Heart volume kg=L
  PVLn <- 7.3e-3 #Brown et al. 1997, Table 4
  VLn <- PVLn * BW #Lung volume kg=L
  V_macro_Ln <- 7.5/100 * VLn #random
  PVBr <- 1.65e-2 #Brown et al. 1997, Table 4
  VBr <- PVBr * BW #Brain volume kg=L
  PVSt <- 6e-3 #Brown et al. 1997, Table 4
  VSt <- PVSt * BW #Stomach volume kg=L
  PVSIn <- 2.53e-2 #Brown et al. 1997, Table 4
  VSIn <- PVSIn * BW #Small intestine volume kg=L
  PVLIn <- 1.09e-2 #Brown et al. 1997, Table 4
  VLIn <- PVLIn * BW #Large intestine volume kg=L
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,86400,60) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df17=========================================================
  
  exp_data <- dataset$df17 # retrieve data of Yang (2010) 4 mg/kg nanographene sheets tissues
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mliver","Mspleen", "Mkidneys", "Mheart","Mlungs",
                    "Mstomach","Mintestine","Mbrain")
  preds_Yang_IV_tissues <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Yang_IV_tissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  obs_Yang_IV_tissues <- list( exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose*VLi*1000/100,
                               exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose*VSpl*1000/100,
                               exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose*VKi*1000/100,
                               exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose*VH*1000/100,
                               exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose*VLn*1000/100,
                               exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose*VSt*1000/100,
                               exp_data[exp_data$Tissue == "Intestine", "mass"]*admin.dose*(VSIn+VLIn)*1000/100,
                               exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose*VBr*1000/100)
  
  
  score[17] <- AAFE(predictions = preds_Yang_IV_tissues, observations = obs_Yang_IV_tissues)
  
  
  ##########################
  #-------------------------
  # Yang et al., 2013
  #-------------------------
  ##########################
  
  # Set up simulations for the 18th case, i.e. Yang (2013) 50 mg/kg tissues nanoGO
  
  BW <- 0.022  # body weight (kg) #mean value of the BW in figure S6
  
  PVLi <- 5.49e-2 #Brown et al. 1997, Table 4
  VLi <- PVLi * BW #liver volume kg=L
  PVSpl <- 3.5e-3 #Brown et al. 1997, Table 4
  VSpl <- PVSpl * BW #Spleen volume kg=L
  V_macro_Spl = 6.94/100*VSpl #https://doi.org/10.1038/s41374-018-0137-1
  PVKi <- 1.67e-2 #Brown et al. 1997, Table 4
  VKi <- PVKi * BW #kidney volume kg=L
  PVH <- 5e-3 #Brown et al. 1997, Table 4
  VH <- PVH * BW #Heart volume kg=L
  PVLn <- 7.3e-3 #Brown et al. 1997, Table 4
  VLn <- PVLn * BW #Lung volume kg=L
  V_macro_Ln <- 7.5/100 * VLn #random
  PVBr <- 1.65e-2 #Brown et al. 1997, Table 4
  VBr <- PVBr * BW #Brain volume kg=L
  PVSt <- 6e-3 #Brown et al. 1997, Table 4
  VSt <- PVSt * BW #Stomach volume kg=L
  PVSIn <- 2.53e-2 #Brown et al. 1997, Table 4
  VSIn <- PVSIn * BW #Small intestine volume kg=L
  PVLIn <- 1.09e-2 #Brown et al. 1997, Table 4
  VLIn <- PVLIn * BW #Large intestine volume kg=L
  
  admin.dose_per_kg <- 50 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "ip"
  sex <- "F" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,10080,40) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df18=========================================================
  
  exp_data <- dataset$df18 # retrieve data of  Yang (2010) 50 mg/kg tissues nanoG
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mliver","Mspleen", "Mkidneys", "Mheart","Mlungs",
                    "Mstomach","Mintestine","Mbrain")
  preds_Yang_IP_tissues <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Yang_IP_tissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  obs_Yang_IP_tissues <- list( exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose*VLi*1000/100,
                               exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose*VSpl*1000/100,
                               exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose*VKi*1000/100,
                               exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose*VH*1000/100,
                               exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose*VLn*1000/100,
                               exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose*VSt*1000/100,
                               exp_data[exp_data$Tissue == "Intestine", "mass"]*admin.dose*(VSIn+VLIn)*1000/100,
                               exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose*VBr*1000/100)
  
  
  score[18] <- AAFE(predictions = preds_Yang_IP_tissues, observations = obs_Yang_IP_tissues)
  
  
  #======================================df19=========================================================

  admin.dose_per_kg <- 100 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "oral"
  sex <- "F" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,1440,5) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df19=========================================================
  
  exp_data <- dataset$df19 # retrieve data of  Yang (2010) 100 mg/kg tissues nanoG
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mliver","Mspleen", "Mkidneys", "Mheart","Mlungs",
                    "Mstomach","Mintestine","Mbrain")
  preds_Yang_OR_tissues <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Yang_OR_tissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  obs_Yang_OR_tissues <- list( exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose*VLi*1000/100,
                               exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose*VSpl*1000/100,
                               exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose*VKi*1000/100,
                               exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose*VH*1000/100,
                               exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose*VLn*1000/100,
                               exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose*VSt*1000/100,
                               exp_data[exp_data$Tissue == "Intestine", "mass"]*admin.dose*(VSIn+VLIn)*1000/100,
                               exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose*VBr*1000/100)
  
  
  score[19] <- AAFE(predictions = preds_Yang_OR_tissues, observations = obs_Yang_OR_tissues)
  
  
  ##########################
  #-------------------------
  # Yang et al., 2013
  #-------------------------
  ##########################
  
  # Set up simulations for the 18th case, i.e. Yang (2010) 50 mg/kg tissues RGO
  
  BW <- 0.025  # body weight (kg) #mean value of the BW in figure S6
  
  PVLi <- 5.49e-2 #Brown et al. 1997, Table 4
  VLi <- PVLi * BW #liver volume kg=L
  PVSpl <- 3.5e-3 #Brown et al. 1997, Table 4
  VSpl <- PVSpl * BW #Spleen volume kg=L
  V_macro_Spl = 6.94/100*VSpl #https://doi.org/10.1038/s41374-018-0137-1
  PVKi <- 1.67e-2 #Brown et al. 1997, Table 4
  VKi <- PVKi * BW #kidney volume kg=L
  PVH <- 5e-3 #Brown et al. 1997, Table 4
  VH <- PVH * BW #Heart volume kg=L
  PVLn <- 7.3e-3 #Brown et al. 1997, Table 4
  VLn <- PVLn * BW #Lung volume kg=L
  V_macro_Ln <- 7.5/100 * VLn #random
  PVBr <- 1.65e-2 #Brown et al. 1997, Table 4
  VBr <- PVBr * BW #Brain volume kg=L
  PVSt <- 6e-3 #Brown et al. 1997, Table 4
  VSt <- PVSt * BW #Stomach volume kg=L
  PVSIn <- 2.53e-2 #Brown et al. 1997, Table 4
  VSIn <- PVSIn * BW #Small intestine volume kg=L
  PVLIn <- 1.09e-2 #Brown et al. 1997, Table 4
  VLIn <- PVLIn * BW #Large intestine volume kg=L
  
  admin.dose_per_kg <- 50 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "ip"
  sex <- "F" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,10080,40) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df20=========================================================
  
  exp_data <- dataset$df20 # retrieve data of  Yang (2010) 100 mg/kg tissues RGO
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mliver","Mspleen", "Mkidneys", "Mheart","Mlungs",
                    "Mstomach","Mintestine","Mbrain")
  preds_Yang_IP_tissues_RGO <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Yang_IP_tissues_RGO[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  obs_Yang_IP_tissues_RGO <- list( exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose*VLi*1000/100,
                                   exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose*VSpl*1000/100,
                                   exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose*VKi*1000/100,
                                   exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose*VH*1000/100,
                                   exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose*VLn*1000/100,
                                   exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose*VSt*1000/100,
                                   exp_data[exp_data$Tissue == "Intestine", "mass"]*admin.dose*(VSIn+VLIn)*1000/100,
                                   exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose*VBr*1000/100)
  
  
  score[20] <- AAFE(predictions = preds_Yang_IP_tissues_RGO, observations = obs_Yang_IP_tissues_RGO)
  
  
  #======================================df21=========================================================
  
  admin.dose_per_kg <- 100 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "oral"
  sex <- "F" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,1440,5) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df21=========================================================
  
  exp_data <- dataset$df21 # retrieve data of  Yang (2010) 100 mg/kg tissues nanoG
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mliver","Mspleen", "Mkidneys", "Mheart","Mlungs",
                    "Mstomach","Mintestine","Mbrain")
  preds_Yang_OR_tissues_RGO <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Yang_OR_tissues_RGO[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  obs_Yang_OR_tissues_RGO <- list( exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose*VLi*1000/100,
                                   exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose*VSpl*1000/100,
                                   exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose*VKi*1000/100,
                                   exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose*VH*1000/100,
                                   exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose*VLn*1000/100,
                                   exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose*VSt*1000/100,
                                   exp_data[exp_data$Tissue == "Intestine", "mass"]*admin.dose*(VSIn+VLIn)*1000/100,
                                   exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose*VBr*1000/100)
  
  
  score[21] <- AAFE(predictions = preds_Yang_OR_tissues_RGO, observations = obs_Yang_OR_tissues_RGO)
  
  
  # Set up simulations for the 22nd case, i.e. Yang (2013) nanoRGO
  
  BW <- 0.026  # body weight (kg) #mean value of the BW in figure S6
  
  PVLi <- 5.49e-2 #Brown et al. 1997, Table 4
  VLi <- PVLi * BW #liver volume kg=L
  PVSpl <- 3.5e-3 #Brown et al. 1997, Table 4
  VSpl <- PVSpl * BW #Spleen volume kg=L
  V_macro_Spl = 6.94/100*VSpl #https://doi.org/10.1038/s41374-018-0137-1
  PVKi <- 1.67e-2 #Brown et al. 1997, Table 4
  VKi <- PVKi * BW #kidney volume kg=L
  PVH <- 5e-3 #Brown et al. 1997, Table 4
  VH <- PVH * BW #Heart volume kg=L
  PVLn <- 7.3e-3 #Brown et al. 1997, Table 4
  VLn <- PVLn * BW #Lung volume kg=L
  V_macro_Ln <- 7.5/100 * VLn #random
  PVBr <- 1.65e-2 #Brown et al. 1997, Table 4
  VBr <- PVBr * BW #Brain volume kg=L
  PVSt <- 6e-3 #Brown et al. 1997, Table 4
  VSt <- PVSt * BW #Stomach volume kg=L
  PVSIn <- 2.53e-2 #Brown et al. 1997, Table 4
  VSIn <- PVSIn * BW #Small intestine volume kg=L
  PVLIn <- 1.09e-2 #Brown et al. 1997, Table 4
  VLIn <- PVLIn * BW #Large intestine volume kg=L
  
  admin.dose_per_kg <- 50 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "ip"
  sex <- "F" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,10080,40) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df22=========================================================
  
  exp_data <- dataset$df22 # retrieve data of  Yang (2010) 50 mg/kg tissues nanoGO
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mliver","Mspleen", "Mkidneys", "Mheart","Mlungs",
                    "Mstomach","Mintestine","Mbrain")
  preds_Yang_IP_tissues_nRGO <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Yang_IP_tissues_nRGO[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  obs_Yang_IP_tissues_nRGO <- list( exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose*VLi*1000/100,
                                    exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose*1000/100,
                                    exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose*VKi*1000/100,
                                    exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose*VH*1000/100,
                                    exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose*VLn*1000/100,
                                    exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose*VSt*1000/100,
                                    exp_data[exp_data$Tissue == "Intestine", "mass"]*admin.dose*(VSIn+VLIn)*1000/100,
                                    exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose*VBr*1000/100)
  
  
  score[22] <- AAFE(predictions = preds_Yang_IP_tissues_nRGO, observations = obs_Yang_IP_tissues_nRGO)
  
  
  #======================================df23=========================================================
  
  admin.dose_per_kg <- 100 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "oral"
  sex <- "F" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,1440,5) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df23=========================================================
  
  exp_data <- dataset$df23 # retrieve data of  Yang (2013) 50 mg/kg tissues nanoGO
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mliver","Mspleen", "Mkidneys", "Mheart","Mlungs",
                    "Mstomach","Mintestine","Mbrain")
  preds_Yang_OR_tissues_nRGO <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Yang_OR_tissues_nRGO[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  obs_Yang_OR_tissues_nRGO <- list( exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose*VLi*1000/100,
                                    exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose*1000/100,
                                    exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose*VKi*1000/100,
                                    exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose*VH*1000/100,
                                    exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose*VLn*1000/100,
                                    exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose*VSt*1000/100,
                                    exp_data[exp_data$Tissue == "Intestine", "mass"]*admin.dose*(VSIn+VLIn)*1000/100,
                                    exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose*VBr*1000/100)
  
  
  score[23] <- AAFE(predictions = preds_Yang_OR_tissues_nRGO, observations = obs_Yang_OR_tissues_nRGO)
  
  

  # Estimate final score
  if (sum(is.na(score))>0){
    final_score <- 100
    
  }else{
    final_score <- mean(score)
    
  }
  return(final_score)
  
}

################################################################################


MW <- 124.91 #g/mol

# Read data
Liu_1_small_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_1_tissues.xlsx")
Liu_1_small_diftp_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_1_tissues_dif_times.xlsx")
Liu_1_large_diftp_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_large_1_tissues_dif_times.xlsx")
Liu_2_small_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_2_tissues.xlsx")
Liu_10_small_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_10_tissues.xlsx")
Liu_1_small_diftp_blood <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_1_blood_dif_times.xlsx")
Liu_1_large_diftp_blood <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_large_1_blood_dif_times.xlsx")
Li_IN_tissues <- openxlsx::read.xlsx("Data/Li_2013_male_tissues_in_nanoGo.xlsx")
Li_IV_tissues <- openxlsx::read.xlsx("Data/Li_2014_male_tissues_iv_nanoGo.xlsx")
Li_IV_blood <- openxlsx::read.xlsx("Data/Li_2014_male_blood_iv_nanoGo.xlsx")
Li_IV_tissues_PEG <- openxlsx::read.xlsx("Data/Li_2014_male_tissues_iv-nanoGo-Peg.xlsx")
Li_IV_blood_PEG <- openxlsx::read.xlsx("Data/Li_2014_male_blood_iv_nanoGo-Peg.xlsx")
Mao_IN_tissues <- openxlsx::read.xlsx("Data/Mao_2016_male_tissues_in.xlsx")
Mao_OR_tissues <- openxlsx::read.xlsx("Data/Mao_2016_male_tissues_oral.xlsx")
Yang_IV_blood <- openxlsx::read.xlsx("Data/Yang_2010_blood_female_iv.xlsx")
Yang_IV_LiSpl <- openxlsx::read.xlsx("Data/Yang_2010_female_liver_spleen_iv.xlsx")
Yang_IV_tissues <- openxlsx::read.xlsx("Data/Yang_2010_female_tissues_iv.xlsx")
Yang_IP_tissues <- openxlsx::read.xlsx("Data/Yang_2013_NanoGO_female_tissues_ip.xlsx")
Yang_OR_tissues <- openxlsx::read.xlsx("Data/Yang_2013_NanoGO_female_tissues_oral.xlsx")
Yang_IP_tissues_RGO <- openxlsx::read.xlsx("Data/Yang_2013_reducedGO_female_tissues_ip.xlsx")
Yang_OR_tissues_RGO <- openxlsx::read.xlsx("Data/Yang_2013_reducedGO_female_tissues_oral.xlsx")
Yang_IP_tissues_nRGO <- openxlsx::read.xlsx("Data/Yang_2013_NanoreducedGO_female_tissues_ip.xlsx")
Yang_OR_tissues_nRGO <- openxlsx::read.xlsx("Data/Yang_2013_NanoreducedGO_female_tissues_oral.xlsx")

setwd("/Users/eviepapakyriakopoulou/Documents/GitHub/Graphene_oxide_Nimble_model/Final model/Optimization/Training/AAFE_NLOPT_LN_SBPLX")


dataset <- list("df1" = Liu_1_small_tissues,"df2" = Liu_1_small_diftp_tissues,
                "df3" = Liu_1_large_diftp_tissues,"df4" = Liu_2_small_tissues,
                "df5" = Liu_10_small_tissues,"df6" = Liu_1_small_diftp_blood,
                "df7" = Liu_1_large_diftp_blood,"df8" = Li_IN_tissues,
                "df9" = Li_IV_tissues,"df10" = Li_IV_blood, "df11" = Li_IV_tissues_PEG,
                "df12" = Li_IV_blood_PEG, "df13" = Mao_IN_tissues, "df14" = Mao_OR_tissues,
                "df15" = Yang_IV_blood, "df16" = Yang_IV_LiSpl, "df17" = Yang_IV_tissues,
                "df18" = Yang_IP_tissues,"df19" = Yang_OR_tissues,
                "df20" = Yang_IP_tissues_RGO,"df21" = Yang_OR_tissues_RGO,
                "df22" = Yang_IP_tissues_nRGO, "df23" = Yang_OR_tissues_nRGO)


#Initialise optimiser to NULL for better error handling later
opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA"
              "xtol_rel" = 1e-03,
              "ftol_rel" = 0.0,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0, 
              "maxeval" = 2000, 
              "print_level" = 1)

# Create initial conditions (zero initialisation)
#Parameter names:

N_pars <- 14 # Number of parameters to be fitted
fit <-  c(rep(log(1), 10), log(1e-5), log(1e-5), log(1), log(1))

lb = c(rep(log(1e-10),14))
ub = c(rep(log(1e10),14))
# 
# lb = c(rep(log(1e-3),15), rep(log(1e-8),5))
# ub = c(rep(log(1e3),15), rep(log(1e8),5))

# Run the optimization algorithm to estimate the parameter values
optimizer <- nloptr::nloptr( x0= fit,
                             eval_f = obj.func,
                             lb	= lb,
                             ub = ub,
                             opts = opts,
                             dataset = dataset)


estimated_params <- exp(optimizer$solution)
save.image("GO_model_optimization_rsmd_NLOPT_LN_SBPLX.RData")



# Set up simulations for the 1st case, i.e. Liu (2012) 1 mg/kg small_tissues
BW <- 0.03  # body weight (kg) #not reported
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 3/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10,0.05) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df1=========================================================

exp_data <- dataset$df1 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_small_tissues <- solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                          "Mkidneys", "Mlungs", "Mbrain",
                                          "Msmall_intestine", "Mlarge_intestine")]



# Set up simulations for the 2nd case, Liu (2012) 1 mg/kg small_tissues_diff_time_points
BW <- 0.03  # body weight (kg) #not reported
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 3/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,180,1) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))
#======================================df2=========================================================

exp_data <- dataset$df2 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg diff_time_points
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_small_diftp_tissues <-  solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                                 "Mkidneys", "Mlungs", "Mbrain",
                                                 "Msmall_intestine", "Mlarge_intestine")]



# Set up simulations for the 3rd case, i.e. Liu (2012) 1 mg/kg large_tissues_diff_time_points

BW <- 0.03  # body weight (kg)
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_large <- 300/2
np_size <- np_size_large #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M"
user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,180,1) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df3=========================================================

exp_data <- dataset$df3 # retrieve data of Liu et al. 2012 tissues large p.s., 1 mg/kg diff_time_points
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_large_diftp_tissues <-  solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                                 "Mkidneys", "Mlungs", "Mbrain",
                                                 "Msmall_intestine", "Mlarge_intestine")]


# Set up simulations for the 4th case, i.e. Liu (2012) 2 mg/kg small_tissues

BW <- 0.03  # body weight (kg) #not reported
admin.dose_per_kg <- 2 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 3/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10,0.05) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))


#======================================df4=========================================================

exp_data <- dataset$df4 # retrieve data of Liu et al. 2012 tissues Small p.s., 2 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
preds_Liu_2_small_tissues <- solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                          "Mkidneys", "Mlungs", "Mbrain",
                                          "Msmall_intestine", "Mlarge_intestine")]


# Set up simulations for the 5th case, Liu (2012) 10 mg/kg small_tissues

BW <- 0.03  # body weight (kg) #not reported
admin.dose_per_kg <- 10 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 3/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10,0.05) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

#======================================df5=========================================================

exp_data <- dataset$df5 # retrieve data of Liu et al. 2012 tissues Small p.s., 10 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
preds_Liu_10_small_tissues <-  solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                            "Mkidneys", "Mlungs", "Mbrain",
                                            "Msmall_intestine", "Mlarge_intestine")]


# Set up simulations for the 6th case, i.e. Liu (2012) 1 mg/kg small_blood diff_time_points

BW <- 0.03  # body weight (kg) #not reported
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 3/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M"


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,180,1) #h

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df6=========================================================

exp_data <- dataset$df6 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_small_diftp_blood <- solution[, c("time", "Mblood")]


# Set up simulations for the 7th case, i.e. Liu (2012) 1 mg/kg large_blood_diff_time_points

BW <- 0.03  # body weight (kg) #not reported
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_large <- 300/2
np_size <- np_size_large #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M"
user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,180,1) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df7=========================================================

exp_data <- dataset$df7 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_large_diftp_blood <- solution[, c("time", "Mblood")]



# Set up simulations for the 8th case, i.e. Li (2013) 10 mg/kg nanoGO_tissues_diff_time_points

BW <- 0.03  # body weight (kg) #not reported
admin.dose_per_kg <- 10 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "in"
sex <- "M" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,720,1) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df8=========================================================

exp_data <- dataset$df8 # retrieve data of Li (2013) 10 mg/kg nanoGO_tissues_diff_time_points
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Li_IN_tissues <- solution[, c("time","Mblood", "Mheart", "Mlungs","Mliver", "Mspleen", "Mkidneys",
                                    "Mstomach","Msmall_intestine", "Mlarge_intestine")]

# Set up simulations for the 9th case, i.e. Li (2014) 5 mg/kg nanoGO_tissues_diff_time_points

BW <- 0.02717  # body weight (kg) #not reported
admin.dose_per_kg <- 5 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,360,1) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df9=========================================================

exp_data <- dataset$df9 # retrieve data of Li (2014) 5 mg/kg nanoGO_tissues_diff_time_points
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Li_IV_tissues <- solution[, c("time","Mblood", "Mheart", "Mlungs","Mliver", "Mspleen", "Mkidneys",
                                    "Mstomach","Msmall_intestine", "Mlarge_intestine")]


# Set up simulations for the 10th case, i.e. Li (2014) 5 mg/kg nanoGO_blood_diff_time_points


# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,1440,5) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df10=========================================================

exp_data <- dataset$df10 # retrieve data of Li (2014) 5 mg/kg nanoGO_blood_diff_time_points
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Li_IV_blood <- solution[, c("time","Mblood")]


# Set up simulations for the 11th case, i.e. Li (2014) 5 mg/kg nanoGO_tissues_diff_time_points

BW <- 0.02717  # body weight (kg) #not reported

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,360,1) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df11=========================================================

exp_data <- dataset$df11 # retrieve data of Li (2014) 5 mg/kg nanoGO_PEG_tissues_diff_time_points
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Li_IV_tissues_PEG <- solution[, c("time","Mlungs","Mliver", "Mspleen")]

# Set up simulations for the 12th case, i.e. Li (2014) 5 mg/kg nanoGO_PEG_blood_diff_time_points


# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,1440,5) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df12=========================================================

exp_data <- dataset$df12 # retrieve data of Li (2014) 5 mg/kg nanoGO_PEG_blood_diff_time_points
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Li_IV_blood_PEG <- solution[, c("time","Mblood")]


# Set up simulations for the 13th case, i.e. Mao (2016) 5 mg/kg nanoGO_tissues_diff_time_points

# # Observed body weight data (Graphene group) in grams, from Figure S4
# bw_days <- c(0, 3, 7, 10, 14, 21, 28)
# bw_values_g <- c(18.5, 22, 25.5, 27, 31.5, 36, 34.5)
# 
# bw_time_min <- bw_days * 1440  # 1 day = 1440 min
# time_min <- seq(0, 40320, by = 1440)
# 
# # Interpolate body weight at each time point, convert to kilograms
# BW_interp_g <- approx(x = bw_time_min, y = bw_values_g, xout = time_min, rule = 2)$y
# BW_interp_kg <- BW_interp_g / 1000  # convert to kg
# 
# # Now use your original structure to simulate BW initialization and tracking
# BW <- BW_interp_kg[1]  # initial BW at time 0 (in kg)
# 
# BW_init <- BW  # store the initialized BW if needed
# 
# # Build the full BW vector (kg), labeled by time in minutes
# BW <- c(BW_init, rep(NA, length(time_min) - 1))
# for (i in 2:length(BW)) {
#   BW[i] <- BW_interp_kg[i]
# }

BW_init <- 0.03
admin.dose <- 5 #ug PFOA
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "inh"
sex <- "M" 

user_input <- list('BW'=BW_init,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,40320,60) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df13=========================================================

exp_data <- dataset$df13 # retrieve data of Li (2014) 5 mg/kg nanoGO_tissues_diff_time_points
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Mao_IN_tissues <- solution[, c("time","Mliver", "Mspleen", "Mstomach","Msmall_intestine",
                                     "Mlarge_intestine","Mlungs", "Mfeces")]


# Set up simulations for the 14th case, i.e. Mao (2016) 5 mg/kg nanoGO_tissues

BW <- 0.03  # body weight (kg) #not reported
admin.dose <- 10 #ug PFOA
admin.type <- "oral"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,4320,60) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df14=========================================================

exp_data <- dataset$df14 # retrieve data of Li (2014) 5 mg/kg nanoGO_tissues
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Mao_OR_tissues <- solution[, c("time","Mstomach","Msmall_intestine",
                                     "Mlarge_intestine", "Mfeces")]


# Set up simulations for the 15th case, i.e. Yang (2010) 4 mg/kg nanographene sheets

BW <- 0.022  # body weight (kg) #mean value of the BW in figure S6
PVB <- 1.7e-3/0.02
VB <- PVB * BW
admin.dose_per_kg <- 4 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "F" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,1440,1) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df15=========================================================

exp_data <- dataset$df15 # retrieve data of Yang (2010) 4 mg/kg nanographene sheets
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Yang_IV_blood <- solution[, c("time","Mblood")]


#======================================df16=========================================================
# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,86400,60) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

#======================================df16=========================================================

exp_data <- dataset$df16 # retrieve data of Yang (2010) 4 mg/kg nanographene sheets liver, spleen
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Yang_IV_LiSpl <- solution[, c("time","Mliver","Mspleen")]

#======================================df17=========================================================
PVLi <- 5.49e-2 #Brown et al. 1997, Table 4
VLi <- PVLi * BW #liver volume kg=L
PVSpl <- 3.5e-3 #Brown et al. 1997, Table 4
VSpl <- PVSpl * BW #Spleen volume kg=L
V_macro_Spl = 6.94/100*VSpl #https://doi.org/10.1038/s41374-018-0137-1
PVKi <- 1.67e-2 #Brown et al. 1997, Table 4
VKi <- PVKi * BW #kidney volume kg=L
PVH <- 5e-3 #Brown et al. 1997, Table 4
VH <- PVH * BW #Heart volume kg=L
PVLn <- 7.3e-3 #Brown et al. 1997, Table 4
VLn <- PVLn * BW #Lung volume kg=L
V_macro_Ln <- 7.5/100 * VLn #random
PVBr <- 1.65e-2 #Brown et al. 1997, Table 4
VBr <- PVBr * BW #Brain volume kg=L
PVSt <- 6e-3 #Brown et al. 1997, Table 4
VSt <- PVSt * BW #Stomach volume kg=L
PVSIn <- 2.53e-2 #Brown et al. 1997, Table 4
VSIn <- PVSIn * BW #Small intestine volume kg=L
PVLIn <- 1.09e-2 #Brown et al. 1997, Table 4
VLIn <- PVLIn * BW #Large intestine volume kg=L

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,86400,60) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

#======================================df17=========================================================

exp_data <- dataset$df17 # retrieve data of Yang (2010) 4 mg/kg nanographene sheets tissues
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Yang_IV_tissues <- solution[, c("time","Mliver","Mspleen", "Mkidneys", "Mheart","Mlungs",
                  "Mstomach","Mintestine","Mbrain")]


##########################
#-------------------------
# Yang et al., 2013
#-------------------------
##########################

# Set up simulations for the 18th case, i.e. Yang (2010) 50 mg/kg tissues nanoGO

BW <- 0.023  # body weight (kg) #mean value of the BW in figure S6

PVLi <- 5.49e-2 #Brown et al. 1997, Table 4
VLi <- PVLi * BW #liver volume kg=L
PVSpl <- 3.5e-3 #Brown et al. 1997, Table 4
VSpl <- PVSpl * BW #Spleen volume kg=L
V_macro_Spl = 6.94/100*VSpl #https://doi.org/10.1038/s41374-018-0137-1
PVKi <- 1.67e-2 #Brown et al. 1997, Table 4
VKi <- PVKi * BW #kidney volume kg=L
PVH <- 5e-3 #Brown et al. 1997, Table 4
VH <- PVH * BW #Heart volume kg=L
PVLn <- 7.3e-3 #Brown et al. 1997, Table 4
VLn <- PVLn * BW #Lung volume kg=L
V_macro_Ln <- 7.5/100 * VLn #random
PVBr <- 1.65e-2 #Brown et al. 1997, Table 4
VBr <- PVBr * BW #Brain volume kg=L
PVSt <- 6e-3 #Brown et al. 1997, Table 4
VSt <- PVSt * BW #Stomach volume kg=L
PVSIn <- 2.53e-2 #Brown et al. 1997, Table 4
VSIn <- PVSIn * BW #Small intestine volume kg=L
PVLIn <- 1.09e-2 #Brown et al. 1997, Table 4
VLIn <- PVLIn * BW #Large intestine volume kg=L

admin.dose_per_kg <- 50 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "ip"
sex <- "F" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10080,40) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df18=========================================================

exp_data <- dataset$df18 # retrieve data of  Yang (2010) 50 mg/kg tissues nanoG
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Yang_IP_tissues <- solution[, c("time","Mliver","Mspleen", "Mkidneys", "Mheart","Mlungs",
                                    "Mstomach","Mintestine","Mbrain")]


#======================================df19=========================================================

admin.dose_per_kg <- 100 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "oral"
sex <- "F" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,1440,5) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df19=========================================================

exp_data <- dataset$df19 # retrieve data of  Yang (2010) 100 mg/kg tissues nanoG
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Yang_OR_tissues <- solution[, c("time","Mliver","Mspleen", "Mkidneys", "Mheart","Mlungs",
                                      "Mstomach","Mintestine","Mbrain")]


# Set up simulations for the 20th case, i.e. Yang (2010) 50 mg/kg tissues RGO

BW <- 0.025  # body weight (kg) #mean value of the BW in figure S6

PVLi <- 5.49e-2 #Brown et al. 1997, Table 4
VLi <- PVLi * BW #liver volume kg=L
PVSpl <- 3.5e-3 #Brown et al. 1997, Table 4
VSpl <- PVSpl * BW #Spleen volume kg=L
V_macro_Spl = 6.94/100*VSpl #https://doi.org/10.1038/s41374-018-0137-1
PVKi <- 1.67e-2 #Brown et al. 1997, Table 4
VKi <- PVKi * BW #kidney volume kg=L
PVH <- 5e-3 #Brown et al. 1997, Table 4
VH <- PVH * BW #Heart volume kg=L
PVLn <- 7.3e-3 #Brown et al. 1997, Table 4
VLn <- PVLn * BW #Lung volume kg=L
V_macro_Ln <- 7.5/100 * VLn #random
PVBr <- 1.65e-2 #Brown et al. 1997, Table 4
VBr <- PVBr * BW #Brain volume kg=L
PVSt <- 6e-3 #Brown et al. 1997, Table 4
VSt <- PVSt * BW #Stomach volume kg=L
PVSIn <- 2.53e-2 #Brown et al. 1997, Table 4
VSIn <- PVSIn * BW #Small intestine volume kg=L
PVLIn <- 1.09e-2 #Brown et al. 1997, Table 4
VLIn <- PVLIn * BW #Large intestine volume kg=L

admin.dose_per_kg <- 50 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "ip"
sex <- "F" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10080,40) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df20=========================================================

exp_data <- dataset$df20 # retrieve data of  Yang (2010) 50 mg/kg tissues RGO
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Yang_IP_tissues_RGO <- solution[, c("time","Mliver","Mspleen", "Mkidneys", "Mheart","Mlungs",
                                      "Mstomach","Mintestine","Mbrain")]


#======================================df21=========================================================

admin.dose_per_kg <- 100 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "oral"
sex <- "F" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,1440,5) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df21=========================================================

exp_data <- dataset$df21 # retrieve data of  Yang (2010) 100 mg/kg tissues RGO
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Yang_OR_tissues_RGO <- solution[, c("time","Mliver","Mspleen", "Mkidneys", "Mheart","Mlungs",
                                      "Mstomach","Mintestine","Mbrain")]


# Set up simulations for the 22nd case, i.e. Yang (2013) 50 mg/kg tissues nanoGO

BW <- 0.023  # body weight (kg) #mean value of the BW in figure S6

PVLi <- 5.49e-2 #Brown et al. 1997, Table 4
VLi <- PVLi * BW #liver volume kg=L
PVSpl <- 3.5e-3 #Brown et al. 1997, Table 4
VSpl <- PVSpl * BW #Spleen volume kg=L
V_macro_Spl = 6.94/100*VSpl #https://doi.org/10.1038/s41374-018-0137-1
PVKi <- 1.67e-2 #Brown et al. 1997, Table 4
VKi <- PVKi * BW #kidney volume kg=L
PVH <- 5e-3 #Brown et al. 1997, Table 4
VH <- PVH * BW #Heart volume kg=L
PVLn <- 7.3e-3 #Brown et al. 1997, Table 4
VLn <- PVLn * BW #Lung volume kg=L
V_macro_Ln <- 7.5/100 * VLn #random
PVBr <- 1.65e-2 #Brown et al. 1997, Table 4
VBr <- PVBr * BW #Brain volume kg=L
PVSt <- 6e-3 #Brown et al. 1997, Table 4
VSt <- PVSt * BW #Stomach volume kg=L
PVSIn <- 2.53e-2 #Brown et al. 1997, Table 4
VSIn <- PVSIn * BW #Small intestine volume kg=L
PVLIn <- 1.09e-2 #Brown et al. 1997, Table 4
VLIn <- PVLIn * BW #Large intestine volume kg=L

admin.dose_per_kg <- 50 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "ip"
sex <- "F" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10080,40) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df22=========================================================

exp_data <- dataset$df22 # retrieve data of  Yang (2010) 50 mg/kg tissues nanoG ip
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Yang_IP_tissues_nRGO <- solution[, c("time","Mliver","Mspleen", "Mkidneys", "Mheart","Mlungs",
                                           "Mstomach","Mintestine","Mbrain")]


#======================================df23=========================================================

admin.dose_per_kg <- 100 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "oral"
sex <- "F" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,1440,5) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df23=========================================================

exp_data <- dataset$df23 # retrieve data of  Yang (2013) 50 mg/kg tissues nanoG oral
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Yang_OR_tissues_nRGO <- solution[, c("time","Mliver","Mspleen", "Mkidneys", "Mheart","Mlungs",
                                           "Mstomach","Mintestine","Mbrain")]



# ######################################################################################
#Plot the predictions against the observations
library(ggplot2) 

# Function that creates a plot given a compartment name and the respective predictions and observations
create.plots <- function(predictions, observations, compartment){  
  #Colours of observations and predictions
  cls <-  c("predictions" = "#56B4E9",  "Observations" = "#D55E00")
  
  ggplot(data = predictions)+
    geom_line( aes_string(x= "Time", y= rlang::expr(!!compartment), 
                          color = '"predictions"'),  size=1.5,alpha = 0.7) +
    geom_point(data=observations, aes_string(x="Time", y= rlang::expr(!!compartment), 
                                             color='"Observations"'), size=4)+
    labs(title = rlang::expr(!!compartment), 
         y = expression("GO mass (" * mu* "g)" ),
         x = "Time (min)")+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_color_manual("", values=cls,
                       guide = guide_legend(override.aes =
                                              list(shape = c(16,NA),
                                                   linetype = c(0,1))))+
    theme_light() + 
    theme(legend.position=c(1,1), 
          legend.justification=c(0, 1), 
          legend.key.size = unit(1.5, 'cm'),  
          legend.title = element_text(size=14),
          axis.title=element_text(size=14),
          legend.text = element_text(size=14)
    )
  
}


# Convert Liu 2012, male small_1_tissues from long to wide format using reshape
experiment1 <- reshape(Liu_1_small_tissues[c("Tissue" ,"Time_min", 
                                             "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment1) <- c("Time",Liu_1_small_tissues$Tissue )

# Convert Liu 2012, male small_1_tissues_dif_times from long to wide format using reshape
experiment2 <- reshape(Liu_1_small_diftp_tissues[c("Tissue" ,"Time_min", 
                                                   "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment2) <- c("Time",unique(Liu_1_small_diftp_tissues$Tissue))

# Convert Liu 2012, male large_1_tissues_dif_times from long to wide format using reshape
experiment3 <- reshape(Liu_1_large_diftp_tissues[c("Tissue" ,"Time_min", 
                                                   "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment3) <- c("Time",unique(Liu_1_large_diftp_tissues$Tissue))

# Convert Liu 2012, male small_2_tissues from long to wide format using reshape
experiment4 <- reshape(Liu_2_small_tissues[c("Tissue" ,"Time_min", 
                                             "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment4) <- c("Time",Liu_2_small_tissues$Tissue)

# Convert Liu 2012, male small_10_tissues from long to wide format using reshape
experiment5 <- reshape(Liu_10_small_tissues[c("Tissue" ,"Time_min", 
                                              "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment5) <- c("Time",Liu_10_small_tissues$Tissue)

# Convert Liu 2012, male small_1_blood_dif_times from long to wide format using reshape
experiment6 <- reshape(Liu_1_small_diftp_blood[c("Tissue" ,"Time_min", 
                                                 "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment6) <- c("Time",unique(Liu_1_small_diftp_blood$Tissue))

# Convert Liu 2012, male large_1_blood_dif_times from long to wide format using reshape
experiment7 <- reshape(Liu_1_large_diftp_blood[c("Tissue" ,"Time_min", 
                                                 "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment7) <- c("Time",unique(Liu_1_large_diftp_blood$Tissue))

# Convert Li 2013, male nanoGO_tissues_diff_time_points from long to wide format using reshape
experiment8 <- reshape(Li_IN_tissues[c("Tissue" ,"Time_min", 
                                                 "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment8) <- c("Time",unique(Li_IN_tissues$Tissue))

# Convert Li 2014, male nanoGO_tissues_diff_time_points from long to wide format using reshape
experiment9 <- reshape(Li_IV_tissues[c("Tissue" ,"Time_min", 
                                                 "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment9) <- c("Time",unique(Li_IV_tissues$Tissue))

# Convert Li 2014, male nanoGO_blood_diff_time_points from long to wide format using reshape
experiment10 <- reshape(Li_IV_blood[c("Tissue" ,"Time_min", 
                                       "%ID/gtissue")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment10) <- c("Time",unique(Li_IV_blood$Tissue))

# Convert Li 2014, male nanoGO_PEG_tissues_diff_time_points from long to wide format using reshape
experiment11 <- reshape(Li_IV_tissues_PEG[c("Tissue" ,"Time_min", 
                                      "%ID")], 
                        idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment11) <- c("Time",unique(Li_IV_tissues_PEG$Tissue))

# Convert Li 2014, male nanoGO_PEG_blood_diff_time_points from long to wide format using reshape
experiment12 <- reshape(Li_IV_blood_PEG[c("Tissue" ,"Time_min", 
                                            "%ID/gtissue")], 
                        idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment12) <- c("Time",unique(Li_IV_blood_PEG$Tissue))

# Convert Mao 2016, male nanoGO_tissues_diff_time_points from long to wide format using reshape
experiment13 <- reshape(Mao_IN_tissues[c("Tissue" ,"Time_min", 
                                          "%ID")], 
                        idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment13) <- c("Time",unique(Mao_IN_tissues$Tissue))

# Convert Mao 2016, male nanoGO_tissues from long to wide format using reshape
experiment14 <- reshape(Mao_OR_tissues[c("Tissue" ,"Time_min", 
                                         "%ID")], 
                        idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment14) <- c("Time",unique(Mao_OR_tissues$Tissue))

# Convert Yang 2010, male nanographene sheets_tissues from long to wide format using reshape
experiment15 <- reshape(Yang_IV_blood[c("Tissue" ,"Time_min", 
                                         "%ID/gtissue")], 
                        idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment15) <- c("Time",unique(Yang_IV_blood$Tissue))

# Convert Yang 2010, male nanographene sheets_tissues from long to wide format using reshape
experiment16 <- reshape(Yang_IV_LiSpl[c("Tissue" ,"Time_min", 
                                        "%ID/gtissue")], 
                        idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment16) <- c("Time",unique(Yang_IV_LiSpl$Tissue))

# Convert Yang 2010, male nanographene sheets_tissues from long to wide format using reshape
experiment17 <- reshape(Yang_IV_tissues[c("Tissue" ,"Time_min", 
                                        "%ID/gtissue")], 
                        idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment17) <- c("Time",unique(Yang_IV_tissues$Tissue))

# Convert Yang 2010, male nanographene sheets_tissues from long to wide format using reshape
experiment18 <- reshape(Yang_IP_tissues[c("Tissue" ,"Time_min", 
                                          "%ID/gtissue")], 
                        idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment18) <- c("Time",unique(Yang_IP_tissues$Tissue))

# Convert Yang 2010, male nanographene sheets_tissues from long to wide format using reshape
experiment19 <- reshape(Yang_OR_tissues[c("Tissue" ,"Time_min", 
                                          "%ID/gtissue")], 
                        idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment19) <- c("Time",unique(Yang_OR_tissues$Tissue))

# Convert Yang 2010, male nanographene sheets_tissues from long to wide format using reshape
experiment20 <- reshape(Yang_IP_tissues_RGO[c("Tissue" ,"Time_min", 
                                          "%ID/gtissue")], 
                        idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment20) <- c("Time",unique(Yang_IP_tissues_RGO$Tissue))

# Convert Yang 2010, male nanographene sheets_tissues from long to wide format using reshape
experiment21 <- reshape(Yang_OR_tissues_RGO[c("Tissue" ,"Time_min", 
                                          "%ID/gtissue")], 
                        idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment21) <- c("Time",unique(Yang_OR_tissues_RGO$Tissue))


# Convert Yang 2013, male nanographene sheets_tissues from long to wide format using reshape
experiment22 <- reshape(Yang_IP_tissues_nRGO[c("Tissue" ,"Time_min", 
                                              "%ID/gtissue")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment22) <- c("Time",unique(Yang_IP_tissues_nRGO$Tissue))

# Convert Yang 2013, male nanographene sheets_tissues from long to wide format using reshape
experiment23 <- reshape(Yang_OR_tissues_nRGO[c("Tissue" ,"Time_min", 
                                              "%ID/gtissue")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment23) <- c("Time",unique(Yang_OR_tissues_nRGO$Tissue))


# Put the experiments in a list
experiments <- list(experiment1 = experiment1, experiment2 = experiment2, experiment3 = experiment3,
                    experiment4 = experiment4, experiment5 = experiment5, experiment6 = experiment6,
                    experiment7 = experiment7, experiment8 = experiment8, experiment9=experiment9,
                    experiment10 = experiment10, experiment11 = experiment11,experiment12 = experiment12,
                    experiment13 = experiment13, experiment14 = experiment14,experiment15 = experiment15,
                    experiment16 = experiment16, experiment17 = experiment17, experiment18 = experiment18,
                    experiment19 = experiment19, experiment20 = experiment20, experiment21 = experiment21,
                    experiment22 = experiment22, experiment23 = experiment23)


# Rename predictions so that they share the same name as the names of the experimental data dataframe
colnames(preds_Liu_1_small_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                         "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine")
colnames(preds_Liu_1_small_diftp_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                               "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine") 
colnames(preds_Liu_1_large_diftp_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                               "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine") 
colnames(preds_Liu_2_small_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                         "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine") 
colnames(preds_Liu_10_small_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                          "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine") 
colnames(preds_Liu_1_small_diftp_blood) <- c("Time", "Blood")
colnames(preds_Liu_1_large_diftp_blood) <- c("Time", "Blood") 
colnames(preds_Li_IN_tissues) <- c("Time", "Blood", "Heart", "Lungs","Liver", "Spleen", "Kidney",
                                   "Stomach","Small_intestine", "Large_intestine") 
colnames(preds_Li_IV_tissues) <- c("Time", "Blood", "Heart", "Lungs","Liver", "Spleen", "Kidney",
                                   "Stomach","Small_intestine", "Large_intestine") 
colnames(preds_Li_IV_blood) <- c("Time", "Blood")
colnames(preds_Li_IV_tissues_PEG) <- c("Time", "Lungs","Liver", "Spleen")
colnames(preds_Li_IV_blood_PEG) <- c("Time", "Blood")
colnames(preds_Mao_IN_tissues) <- c("Time", "Liver", "Spleen", "Stomach","Small_intestine",
                                    "Large_intestine","Lungs", "Feces")
colnames(preds_Mao_OR_tissues) <- c("Time", "Stomach","Small_intestine",
                                    "Large_intestine","Feces")
colnames(preds_Yang_IV_blood) <- c("Time", "Blood")
colnames(preds_Yang_IV_LiSpl) <- c("Time", "Liver", "Spleen")
colnames(preds_Yang_IV_tissues) <- c("Time", "Liver","Spleen", "Kidneys", "Heart","Lungs",
                                     "Stomach","Intestine","Brain")
colnames(preds_Yang_IP_tissues) <- c("Time", "Liver","Spleen", "Kidneys", "Heart","Lungs",
                                     "Stomach","Intestine","Brain")
colnames(preds_Yang_OR_tissues) <- c("Time", "Liver","Spleen", "Kidneys", "Heart","Lungs",
                                     "Stomach","Intestine","Brain")
colnames(preds_Yang_IP_tissues_RGO) <- c("Time", "Liver","Spleen", "Kidneys", "Heart","Lungs",
                                     "Stomach","Intestine","Brain")
colnames(preds_Yang_OR_tissues_RGO) <- c("Time", "Liver","Spleen", "Kidneys", "Heart","Lungs",
                                     "Stomach","Intestine","Brain")
colnames(preds_Yang_IP_tissues_nRGO) <- c("Time", "Liver","Spleen", "Kidneys", "Heart","Lungs",
                                          "Stomach","Intestine","Brain")
colnames(preds_Yang_OR_tissues_nRGO) <- c("Time", "Liver","Spleen", "Kidneys", "Heart","Lungs",
                                          "Stomach","Intestine","Brain")

# Create a list containing the corresponding predictions
simulations <- list(predictions1 = preds_Liu_1_small_tissues, predictions2 = preds_Liu_1_small_diftp_tissues,
                    predictions3 = preds_Liu_1_large_diftp_tissues, predictions4 = preds_Liu_2_small_tissues,
                    predictions5 = preds_Liu_10_small_tissues, predictions6 = preds_Liu_1_small_diftp_blood,
                    predictions7 = preds_Liu_1_large_diftp_blood, predictions8 = preds_Li_IN_tissues,
                    predictions9 = preds_Li_IV_tissues, predictions10 = preds_Li_IV_blood,
                    predictions11 = preds_Li_IV_tissues_PEG, predictions12 = preds_Li_IV_blood_PEG,
                    predictions13 = preds_Mao_IN_tissues, predictions14 = preds_Mao_OR_tissues,
                    predictions15 = preds_Yang_IV_blood, predictions16 = preds_Yang_IV_LiSpl,
                    predictions17 = preds_Yang_IV_tissues, predictions18 = preds_Yang_IP_tissues,
                    predictions19 = preds_Yang_OR_tissues, predictions20 = preds_Yang_IP_tissues_RGO,
                    predictions21 = preds_Yang_OR_tissues_RGO, predictions22 = preds_Yang_IP_tissues_nRGO,
                    predictions23 = preds_Yang_OR_tissues_nRGO)


# Iterate over all existing experiments and create the accompanying plots
for(i in 1:length(experiments)){
  # Retrieve the corresponding observations and simulations
  observations <- experiments[[i]]
  predictions <- simulations[[i]]
  # Extract the compartment names
  compartments <- names(predictions)[2:length(predictions)]
  
  # Use lapply to iterate over the column names and create plots
  plots <- lapply(compartments, function(compartment) {
    create.plots(predictions, observations, compartment )
  })
  if(length(compartments) == 1){
    final_plot <- do.call(ggpubr::ggarrange, c(plots, ncol = 1, nrow = 1,
                                               common.legend = TRUE, legend = "right"))
    
  }else{
    final_plot <- do.call(ggpubr::ggarrange, c(plots, ncol = 3, nrow = ceiling(length(plots) / 3),
                                               common.legend = TRUE, legend = "right"))
  }
  
  
  plot.margin=unit(c(0,0,0,0), "pt")
  
  
  # Save the plot with dynamically adjusted dimensions
  ggsave(paste0("experiment", i,".png"), plot = final_plot,
         device = 'png', dpi = 300,
         width = 13,
         height = 10,
         units = "in")
}


