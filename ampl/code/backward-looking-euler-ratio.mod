  #-----------growth rate of capital as a parameter
  if s > PInf then {
    for {r in Regions, i in Sectors}{
      let GROWTH_KAP[r, i, s]
        := (KAP[r, i, s] - KAP[r, i, s - 1]) / KAP[r, i, s];
      let GROWTH_OUT[r, i, s]
        := (E_OUT[r, i, s] - E_OUT[r, i, s - 1]) / E_OUT[r, i, s - 1];
      if s > PInf + 1 then {
        #-----------Euler integrand for CES production
        let EULER_INTEGRAND[r, i, s - 1]
          := DUAL_KAP[r, i, s - 1] * (1 - DELTA[i]) 
            + DUAL_MKT_CLR[i, s - 1] * (
              SHR_KAP_OUT_CES[r, i] * KAP[r, i, s - 1] ** (RHO_OUT[i] - 1)
              * SCALE_OUT * A[i] * SHR_DOM[r, i, s - 1]
              * (DOM[r, i, s - 1] / (A[i] * SHR_DOM[r, i, s - 1]))
                  ** (1 - RHO_OUT[i] / SCALE_OUT)
            - PHI_ADJ[i]
              * (2 * GROWTH_KAP[r, i, s] + GROWTH_KAP[r, i, s] ** 2)
            );
        let EULER_RATIO[r, i, s - 1] 
          := BETA * EULER_INTEGRAND[r, i, s - 1] / DUAL_KAP[r, i, s - 2];
        #-----------Euler integrand for Cobb--Douglas production
        #let EULER_INTEGRAND[r, i, s] :=  DUAL_KAP[r, i, s] * (1 - DELTA[i]) 
        #  + DUAL_MKT_CLR[i, s] * (
        #    SHR_KAP_OUT[i]
        #    * (KAP[r, i, s] / LAB[r, i, s]) ** (SHR_KAP_OUT[i] - 1)
        #  - PHI_ADJ[i] * (2 * GROWTH_KAP[r, i, s] + GROWTH_KAP[r, i, s] ** 2)
        #  );
      };
    };
  };

