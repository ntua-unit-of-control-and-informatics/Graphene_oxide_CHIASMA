#Graphene oxide model

library(deSolve)
library(ggplot2)
library(jaqpotr)

#=========================
#1. Parameters of the model
#=========================

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
                "Qtotal"=Qtotal, "QLitot"=QLitot, "QGE"=QGE,
                "sex"=sex, 'Vper'=Vper,
                
                "admin.time" = admin.time, "admin.dose" = admin.dose,
                "admin.type" = admin.type, "MW"=MW, "study"=study
                
  
                
                
    ))
  })
}  

#===============================================
#2. Function to create initial values for ODEs 
#===============================================

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


#================================================================================================
#===================
#3. Events function
#===================

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


#================================================================================================
#==================
#4. Custom function 
#==================

custom.func <- function(){
  return()
}


# ============================================================================
#5. ODEs System
# ============================================================================


ode.func <- function(time, inits, params, custom.func){
  with(as.list(c(inits, params)),{
    
    if(study == "study_1"){
      coef_liver <- 1.693438e-01
      coef_spleen <- 3.455625e-01
      coef_kidney <- 7.775888e-03
      coef_heart <- 5.195568e-04
      coef_lung <- 4.960860e-01
      coef_brain <- 2.626215e-06
      coef_rob <- 1.243346e-01
      coef_stomach <- 6.284008e-02
      coef_smallIn <- 3.363185e-02
      coef_largeIn <- 3.363185e-02
      CLurine <-0
      CLfeces <- 6.585924e-07
      kabs_oral <- 0
      kabs_ip <- 0
      
    }else if(study == "study_2a") {
      coef_liver <- 2.070579e+00
      coef_spleen <- 2.409691e+00
      coef_kidney <- 5.082606e-01
      coef_heart <- 2.622628e-01
      coef_lung <- 6.892884e+02
      coef_brain <- 4.588473e-02
      coef_rob <- 8.690127e-01
      coef_stomach <- 8.585562e-01
      coef_smallIn <- 5.787620e-02
      coef_largeIn <- 5.429763e-01
      CLurine <-0
      CLfeces <- 2.003016e-05
      kabs_oral <- 0
      kabs_ip <- 0
      
    
   }else if (study == "study_2b") {
     coef_liver <- 1.917379e+00
     coef_spleen <- 4.953714e+00
     coef_kidney <- 1.664364e+00
     coef_heart <- 1.214740e+00
     coef_lung <- 6.786665e+03
     coef_brain <- 5.044096e-01
     coef_rob <- 1.530830e+00
     coef_stomach <- 5.980797e+00
     coef_smallIn <- 8.371671e-01
     coef_largeIn <- 1.003921e+00
     CLurine <-0
     CLfeces <- 2.457810e-04
     kabs_oral <- 0
     kabs_ip <- 0

  }else if (study == "study_3a"){
    coef_liver <- 8.908584e-01
    coef_spleen <- 6.548940e-01
    coef_kidney <- 6.925459e-02
    coef_heart <- 4.381825e-02
    coef_lung <- 6.580230e+00
    coef_brain <- 8.508300e-06
    coef_rob <- 4.132683e-01
    coef_stomach <- 1.670077e-01
    coef_smallIn <- 1.475720e-01
    coef_largeIn <- 4.934268e-04
    CLurine <-0
    CLfeces <- 1.552319e-14
    kabs_oral <- 0
    kabs_ip <- 3.120356e-02
    
  }else if (study == "study_3b"){
    coef_liver <- 6.636678e-02
    coef_spleen <- 8.248976e-02
    coef_kidney <- 1.096440e-01
    coef_heart <- 4.906516e-02
    coef_lung <- 1.373212e+01
    coef_brain <- 7.674925e-03
    coef_rob <- 1.185487e-02
    coef_stomach <- 6.599115e-01
    coef_smallIn <- 2.031099e-02
    coef_largeIn <- 4.433422e-01
    CLurine <-0
    CLfeces <- 2.529264e-05
    kabs_oral <- 2.162722e-03
    kabs_ip <- 0
  
  }else if (study == "study_3c"){
    coef_liver <- 2.755530e-01
    coef_spleen <- 3.011628e-01
    coef_kidney <- 2.583104e-02
    coef_heart <- 7.857248e-03
    coef_lung <- 1.397005e+00
    coef_brain <- 1.654621e-06
    coef_rob <- 4.796149e-02
    coef_stomach <- 2.644745e-02
    coef_smallIn <- 3.588826e-02
    coef_largeIn <- 1.927583e-02
    CLurine <-0
    CLfeces <- 4.845819e-10
    kabs_oral <- 0
    kabs_ip <- 7.175721e-02
    
  }else if (study == "study_3d"){
    coef_liver <- 1.530302e+00
    coef_spleen <- 1.926719e+00
    coef_kidney <- 6.191944e+00
    coef_heart <- 2.806080e+00
    coef_lung <- 7.830945e+00
    coef_brain <- 9.425943e-01
    coef_rob <- 2.134978e+09
    coef_stomach <- 8.850467e-01
    coef_smallIn <- 3.608792e+00
    coef_largeIn <- 2.437912e+01
    CLurine <-0
    CLfeces <- 4.286406e-06
    kabs_oral <- 5.194047e-01
    kabs_ip <- 0
    
  }else if (study == "study_3e"){
    coef_liver <- 3.094069e+00
    coef_spleen <- 2.838648e+02
    coef_kidney <- 2.504434e-01
    coef_heart <- 9.844240e-02
    coef_lung <- 2.621535e+01
    coef_brain <- 9.924987e-02
    coef_rob <- 8.727455e-02
    coef_stomach <- 4.967337e-01
    coef_smallIn <- 6.450172e-02
    coef_largeIn <- 8.675520e-01
    CLurine <-0
    CLfeces <- 2.356945e-09
    kabs_oral <- 0
    kabs_ip <- 5.995203e-01
    
  }else if (study == "study_3f"){
    coef_liver <- 1.016009e-03
    coef_spleen <- 7.400028e+00
    coef_kidney <- 1.105917e-03
    coef_heart <- 1.005044e-03
    coef_lung <- 1.191686e-01
    coef_brain <- 1.946798e-03 
    coef_rob <- 4.010782e-04
    coef_stomach <- 4.337152e-03
    coef_smallIn <- 5.613834e-03
    coef_largeIn <- 8.626687e-07
    CLurine <-0
    CLfeces <- 7.333796e-04
    kabs_oral <- 1.062609e-01
    kabs_ip <- 0
    
  }else if (study == "study_4"){
    coef_liver <- 5.010207e-02
    coef_spleen <- 1.658549e-01
    coef_kidney <- 2.548485e-01
    coef_heart <- 1.608659e-01
    coef_lung <- 3.547853e+02
    coef_brain <- 8.617769e-01
    coef_rob <- 1.458096e+00
    coef_stomach <- 8.606450e-01
    coef_smallIn <- 1.027909e+00
    coef_largeIn <- 1.164971e-03
    CLurine <-0
    CLfeces <- 1.321398e-06
    kabs_oral <- 0
    kabs_ip <- 0
    
  }else if (study == "study_5a"){
    coef_liver <- 1.123089e-01
    coef_spleen <- 9.805037e-01
    coef_kidney <- 1.451841e-01
    coef_heart <- 2.149089e-01
    coef_lung <- 1.566479e+03
    coef_brain <- 0
    coef_rob <- 1.357213e-02
    coef_stomach <- 6.389196e-01
    coef_smallIn <- 8.734302e-02
    coef_largeIn <- 2.516453e-01
    CLurine <-0
    CLfeces <- 6.225655e-05
    kabs_oral <- 0
    kabs_ip <- 0
    
  }else if (study == "study_5b"){
    coef_liver <- 4.394879e-01
    coef_spleen <- 1.761044e-01
    coef_kidney <- 1.783479e-01
    coef_heart <- 2.563641e-01
    coef_lung <- 1.105172e+02
    coef_brain <- 0
    coef_rob <- 3.740028e-06
    coef_stomach <- 1.710961e-02
    coef_smallIn <- 2.580722e+00
    coef_largeIn <- 9.236887e-02
    CLurine <-0
    CLfeces <- 1.335776e-05
    kabs_oral <- 0
    kabs_ip <- 0
    
  }else if (study == "study_6a"){
    coef_liver <- 1.093465e-01
    coef_spleen <- 3.138718e-01
    coef_kidney <- 0
    coef_heart <- 0
    coef_lung <- 9.792661e+03
    coef_brain <- 0
    coef_rob <- 1.453652e-01
    coef_stomach <- 1.731266e+00
    coef_smallIn <- 6.052005e-01
    coef_largeIn <- 2.063296e+00
    CLurine <-0
    CLfeces <- 7.643525e-06
    kabs_oral <- 0
    kabs_ip <- 0
    
  }else if (study == "study_6b"){
    coef_liver <- 1.477630e+00
    coef_spleen <- 4.832864e-01
    coef_kidney <- 0
    coef_heart <- 0
    coef_lung <- 7.624179e+00
    coef_brain <- 2.487163e-01
    coef_rob <- 1.147408e+01
    coef_stomach <- 0
    coef_smallIn <- 5.990568e+00
    coef_largeIn <- 2.121434e+01
    CLurine <-0
    CLfeces <- 2.373692e-05
    kabs_oral <- 2.653116e-03
    kabs_ip <- 0


   }
     
    
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

#=============
#6. User input 
#=============

BW <- 0.03  # body weight (kg) #not reported
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 
study <- "study_1"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "sex" = sex,
                   "study"=study)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10,0.05) #min

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                           y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-07, atol = 1e-07))


# ====================
# Upload on Jaqpot
# ===================


# Subset of features to be displayed on the user interface
predicted.feats <- c("MBart", "MBven","MKi",
                     "MLi", "MSt", "MStL",
                     "MSIn", "MLIn", 
                     "MLn", "MSpl", "MH", "MBr",
                     "MRe", "Murine", 
                     "Mfeces",  "Mper",
                     "CBven", "CBart","CKi", 
                     "CLi", "CSt", "CSIn", 
                     "CLIn", "CLn", "CSpl", 
                     "CH", "CBr","CRe")
                     
setwd("/Users/eviepapakyriakopoulou/Documents/GitHub/Graphene_oxide_Nimble_model/Final model")

jaqpotr::deploy.pbpk(user.input = user_input,out.vars = predicted.feats,
                     create.params = create.params,  create.inits = create.inits,
                     create.events = create.events, custom.func = custom.func,
                     envFile = "/Users/eviepapakyriakopoulou/Documents/Documents/NTUA Post-Doc/INSIGHT/jaqpotkeys/keys.env")