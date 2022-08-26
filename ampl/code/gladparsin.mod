param TInf 'infimum of the set of times' default 0;
param LInf 'start period/infimum of the look forward set' default 0;
param LSup 'supremum of the look forward set' > LInf default 48;
param PInf 'infimum of times on a path' default 0;
param PSup 'supremum of times on a path (T_star in CJ)' >= PInf default 3;
param TSup 'final time' >= PInf + LInf default PSup + LSup;
set Regions;
set Sectors ordered;
set LookForward 'planning set' = LInf .. LSup + LInf - 1 ordered;
set LookForwardClosure 'includes extra time for closing out the problem' =
  LInf .. LSup + LInf ordered;
set PathTimes 'times along a path' = PInf .. PSup + PInf - 1 ordered;
set PathTimesClosure 'includes extra time for closing out the problem' =
  PInf .. PSup + PInf ordered;
set PathSpace 'path space: the set of paths' default PInf .. PSup + 1;
set AllTimes 'all time periods associated with a path' = TInf .. TSup;
param StrtTm 'start time for each step on each path'{PathTimes, PathSpace}
   default time();
param EndTm 'start time for each step on each path'{PathTimes, PathSpace}
   default time();
param ALPHA 'trend growth factor' >= 1 default 1;
param ALPHA_0 'initial trend growth factor' >= 1 default 1;
param BETA 'discount factor'  in  interval(0, 1) default 0.975;
param DELTA 'rate of depreciation for kapital'{Sectors}  >= 0 default 0.025;
param PHI_ADJ 'kapital adjustment cost'{Sectors}  >= 0 default 0.5;
param GAMMA 'intertemporal elasticity of subst.'{Regions}  >= 0 default 0.5;
param EPS_LAB 'Frisch elasticity of labour supply' >= 0 default 0.5;
param EPS_CON 'elasticity of substitution for consumption flows' > 0
   default 0.1;
param EPS_OUT 'elasticity of substitution for output components' > 0
   default 0.5;
param EPS_MED 'elasticity of substitution for intermediate flows' > 0
   default 0.1;
param EPS_INV 'elasticity of substitution for kapital flows' > 0 default 0.1;
param SCALE_INV 'investment scale parameter' default 0.9;
param SCALE_CON 'consumption scale parameter' default 0.9;
param SCALE_OUT 'output scale parameter' default 0.9;
param SCALE_MED 'intermediate scale parameter' default 0.9;
param SCALE_LAB 'labour scale parameter' default 3;
param INV_MIN 'lower bound of investment' >= 0 default 0;
param ZETA1 'TFP before shock' >= 0 default 1;
param ZETA2 'TFP after shock' >= 0 default 1;
param PROB2 'one period probability of jump in TFP' >= 0 default 0.01;
param TAIL_SHR_CON 'tail consumption share (of output)' >= 0 default 0.45;
param UInf 'infimum of interval for uniform dbn'  in  interval[0, 0.5)
   default 0.4999;
param USup 'supremum of interval for uniform dbn'  in  interval(0.5, 1]
   default 0.5001;
param VInf 'infimum of interval for basic variables' default 0.0001;
param VSup 'supremum of interval for basic variables' default 1e+05;
param OInf 'infimum of interval for observed/actual values' default 1e-07;
param OSup 'supremum of interval for observed/actual values' default 1e+07;
param LabSup 'supremum of interval for labour values' default 0.66;
param RAW_CON_FLW 'raw consumption flows: table8Q1'{Regions, Sectors}  >= 0
   default Uniform(UInf, USup);
param RAW_INV_FLW 'raw investment flows: tablekapflw'{Regions, Sectors,
  Sectors}  >= 0 default Uniform(UInf, USup);
param RAW_LAB_FLW 'raw number of workers: abs'{Regions, Sectors}  >= 0
   default Uniform(UInf, USup);
param RAW_MED_FLW 'raw intermediate flows: table8'{Regions, Sectors,
  Sectors}  >= 0 default Uniform(UInf, USup);
param RAW_KAP_OUT 'raw capital primary input flows: table8P2'{Sectors}  >= 0
   default Uniform(UInf, USup);
param RAW_LAB_OUT
   'raw compensation of employees primary input flows: table8P1'{Sectors}
   >= 0 default Uniform(UInf, USup);
param RAW_MED_OUT 'raw total intermediate input flows: table8T1'{Sectors}
   >= 0 default Uniform(UInf, USup);
param RAW_DOM_CCON 'raw domestic flows to consumption: table5Q1'{Regions,
  Sectors}  >= 0 default Uniform(UInf, USup);
param RAW_YSA_CCON 'raw import flows to consumption: table(8-5)'{r in
  Regions, i in Sectors}  >= 0 default Uniform(UInf, USup);
param RAW_DOM_CINV 'raw domestic flows to investment: table5Q3+Q4+Q5'{r in
  Regions, i in Sectors, j in Sectors}  >= 0 default Uniform(UInf, USup);
param RAW_YSA_CINV 'raw import flows to investment: table(8-5)'{r in
  Regions, i in Sectors, j in Sectors}  >= 0 default Uniform(UInf, USup);
param RAW_DOM_CMED 'raw domestic flows to intermediates: table5'{r in
  Regions, i in Sectors, j in Sectors}  >= 0 default Uniform(UInf, USup);
param RAW_YSA_CMED 'raw import flows to intermediates: table(8-5)'{r in
  Regions, i in Sectors, j in Sectors}  >= 0 default Uniform(UInf, USup);
param RAW_EXO_JOUT 'raw export data: table8Q7'{Regions, Sectors}  default
  Uniform(UInf, USup)*0.1;
param RAW_DOM_JOUT 'raw total domestic uses: table8T6-Q7'{Regions, Sectors}
   default Uniform(UInf, USup)*0.9;
param CON 'observed consumption'{Regions, Sectors, PathTimes}  default 1;
param INV_SEC 'observed investment'{Regions, Sectors, PathTimes}  default 1;
param INV_SUM 'observed total investment'{Regions, Sectors, PathTimes}
   default 1;
param MED_SUM 'observed total intermediate flows'{Regions, Sectors,
  PathTimes}  default 1;
param LAB 'observed labour'{Regions, Sectors, PathTimes}  default 1;
param LAB_EXT 'observed laborforce participation'{Regions, Sectors,
  PathTimes}  default 1;
param KAP 'observed kapital'{Regions, Sectors, PathTimesClosure}  default 1;
param E_OUT 'observed Exp. output'{Regions, Sectors, PathTimes}  default 1;
param ADJ_COST_KAP 'observed adjustment costs of kapital'{Regions, Sectors,
  PathTimes}  default 0;
param MKT_CLR 'observed output'{Sectors, PathTimes}  default 0;
param DUAL_KAP 'lagrange multiplier for kapital accumulation'{Regions,
  Sectors, PathTimes}  default 1;
param DUAL_MKT_CLR 'lagrange multiplier for market clearing constraint'
  {Sectors, PathTimes}  default 1;
param GROWTH_KAP 'observed growth rate for kapital'{Regions, Sectors,
  PathTimes}  default 0.05;
param GROWTH_OUT 'observed growth rate for output'{Regions, Sectors,
  PathTimes}  default 0.05;
param EULER_INTEGRAND 'Euler error integrand'{Regions, Sectors,
  PathTimesClosure}  default 1;
param EULER_RATIO 'Expected Euler ratio'{Regions, Sectors, PathTimes}
   default 1;
param NAIRE 'Non-accelerating rate of employment'{Regions, Sectors,
  PathTimes}  default 0.95;
param EXP_LAB_EXT 'Exponent of lab_ext_sec'{Regions, Sectors, PathTimes}
   default 2;
param DOM 'actual path values for domestic output'{Regions, Sectors,
  PathTimes}  default 1;
param GAMMA_HAT 'utility parameter'{r in Regions}  = 1 - 1/GAMMA[r];
param RHO_LAB 'labour exponent parameter' = 1 + 1/EPS_LAB;
param RHO_LAB_HAT 'inverse of labour exponent parameter' = 1/RHO_LAB;
param RHO_INV 'exponent of the investment ces aggregator' = 1 - 1/EPS_INV;
param RHO_INV_HAT 'inverse of RHO_INV' = 1/RHO_INV;
param RHO_MED 'exponent of the intermediate ces aggregator' = 1 - 1/EPS_MED;
param RHO_MED_HAT 'inverse of RHO_MED' = 1/RHO_MED;
param RHO_OUT 'exponent of the output ces aggregator' = 1 - 1/EPS_OUT;
param RHO_OUT_HAT 'inverse of RHO_OUT' = 1/RHO_OUT;
param RHO_CON 'exponent of the CON ces aggregator' = 1 - 1/EPS_CON;
param RHO_CON_HAT 'inverse of RHO_CON' = 1/RHO_OUT;
param SHR_KAP_OUT 'importance of capital in production'{i in Sectors}  =
  RAW_KAP_OUT[i]/(RAW_KAP_OUT[i] + RAW_LAB_OUT[i] + RAW_MED_OUT[i]);
param SHR_LAB_OUT 'share of labour in output'{i in Sectors}  = RAW_LAB_OUT[i]/
  (RAW_KAP_OUT[i] + RAW_LAB_OUT[i] + RAW_MED_OUT[i]);
param SHR_KAPLAB_OUT 'combined importance of kapital and labour in output'
  {i in Sectors}  = SHR_KAP_OUT[i] + SHR_LAB_OUT[i];
param SHR_MED_OUT 'share of intermediates in output'{i in Sectors}  = 1 -
  SHR_KAPLAB_OUT[i];
param SHR_KAP_OUT_CES 'importance of kapital in production'{i in Sectors}  =
  SHR_KAP_OUT[i]^(1/EPS_OUT);
param SHR_LAB_OUT_CES 'importance of labour in production'{i in Sectors}  =
  SHR_LAB_OUT[i]^(1/EPS_OUT);
param SHR_KAPLAB_OUT_CES
   'combined importance of kapital and labour in prod'{i in Sectors}  =
  SHR_KAPLAB_OUT[i]^(1/EPS_OUT);
param SHR_MED_OUT_CES 'importance of intermediates in production'{i in
  Sectors}  = SHR_MED_OUT[i]^(1/EPS_OUT);
param A 'productivity trend'{i in Sectors}  default 1;
param A_LAB
   'importance of disutility of labour (weight in utility function)'{
  Regions, LookForwardClosure}  default -1;
param A_LAB_EXT 'disutility weight for labourforce deviations in utility'
   default -1;
param A_CON 'importance of consumption in utility' default 1;
param A_VAL 'Calibration factor for terminal value function' default 1;
param A_INV 'Calibration factor for investment' default 1;
param A_MED 'Calibration factor for intermediate bundle' default 1;
param REG_WGHT 'regional (population) weights'{r in Regions}  default 1/
  card(Regions);
param CON_FLW_SUM{r in Regions}  = sum{i in Sectors} RAW_CON_FLW[r,i];
param SHR_CON
   'consumption weights for each good in utility (for Cobb-Doug)'{r in
  Regions, i in Sectors}  = RAW_CON_FLW[r,i]/CON_FLW_SUM[r];
param SHR_CON_CES 'CES consumption weights for each good in utility'{r in
  Regions, i in Sectors}  = (RAW_CON_FLW[r,i]/CON_FLW_SUM[r])^(1/EPS_CON);
param LAB_FLW_SUM{r in Regions}  = sum{i in Sectors} RAW_LAB_FLW[r,i];
param SHR_LAB 'labour weights for each sector in utility'{r in Regions,
  j in Sectors}  = RAW_LAB_FLW[r,j]/LAB_FLW_SUM[r];
param INV_FLW_RSUM{r in Regions, j in Sectors}  = sum{i in Sectors}
  RAW_INV_FLW[r,i,j];
param SHR_INV_CES "sectoral share of i in j's CES investment aggregator"
  {r in Regions, i in Sectors, j in Sectors}  = (RAW_INV_FLW[r,i,j]/
  INV_FLW_RSUM[r,j])^(1/EPS_INV);
param MED_FLW_RSUM{r in Regions, j in Sectors}  = sum{i in Sectors}
  RAW_MED_FLW[r,i,j];
param SHR_MED_CES "sectoral share of i in j's CES intermediate aggregator"
  {r in Regions, i in Sectors, j in Sectors}  = (RAW_MED_FLW[r,i,j]/
  MED_FLW_RSUM[r,j])^(1/EPS_INV);
param SHR_LAB_CES 'sectoral share of i in the CES labour aggregator'{r in
  Regions, i in Sectors}  = SHR_LAB[r,i]^(-1/EPS_LAB);
param Pr_shk 'probability of SHK'{Regions, Sectors, t in LookForward}  = (1 -
  PROB2)^t;
param E_shk 'expected shock (exogenous)'{r in Regions, i in Sectors, t in
  LookForward}  = ZETA2 + Pr_shk[r,i,t]*(ZETA1 - ZETA2);
param PRC_YSA 'import prices'{Sectors, LookForward}  default 1;
param PRC_EXA 'export prices'{Sectors, LookForward}  default 1;
param A_CCON 'calibration factor for composite consumption' default 1;
param A_CMED 'scaling factor for composite intermediates' default 1;
param A_CINV 'scaling factor for composite intermediates' default 1;
param SCALE_CINV 'economies of scale for composite intermediate' default 1;
param SCALE_CMED 'economies of scale for composite intermediate' default 1;
param EPS_CINV 'elasticity of subst. for composite investment flows'{
  Regions, Sectors, Sectors}  default 2;
param EPS_CMED 'elasticity of subst. for composite intermediates'{Regions,
  Sectors, Sectors}  default 2;
param RHO_CINV 'CES exponent for composite intermediates'{r in Regions,
  i in Sectors, j in Sectors}  = 1 - 1/EPS_CINV[r,i,j];
param RHO_CINV_HAT 'inverse of CES exponent for composite intermediates'
  {r in Regions, i in Sectors, j in Sectors}  = 1/RHO_CINV[r,i,j];
param RHO_CMED 'CES exponent for composite intermediates'{r in Regions,
  i in Sectors, j in Sectors}  = 1 - 1/EPS_CMED[r,i,j];
param RHO_CMED_HAT 'inverse of CES exponent for composite intermediates'
  {r in Regions, i in Sectors, j in Sectors}  = 1/RHO_CMED[r,i,j];
param SHR_DOM_CCON 'domestic share in comp. consumption'{r in Regions,
  i in Sectors}  = RAW_DOM_CCON[r,i]/(RAW_DOM_CCON[r,i] + RAW_YSA_CCON[r,i]);
param SHR_YSA_CCON 'domestic share in comp. consumption'{r in Regions,
  i in Sectors}  = 1 - SHR_DOM_CCON[r,i];
param SHR_DOM_CINV 'domestic share in composite investment flows'{r in
  Regions, i in Sectors, j in Sectors}  = RAW_DOM_CINV[r,i,j]/(
  RAW_DOM_CINV[r,i,j] + RAW_YSA_CINV[r,i,j]);
param SHR_YSA_CINV 'import share in composite investment flows'{r in
  Regions, i in Sectors, j in Sectors}  = 1 - SHR_DOM_CINV[r,i,j];
param SHR_DOM_CMED 'domestic share in composite intermediate flows'{r in
  Regions, i in Sectors, j in Sectors}  = RAW_DOM_CMED[r,i,j]/(
  RAW_DOM_CMED[r,i,j] + RAW_YSA_CMED[r,i,j]);
param SHR_YSA_CMED 'import share in composite intermediate flows'{r in
  Regions, i in Sectors, j in Sectors}  = 1 - SHR_DOM_CMED[r,i,j];
param SHR_EXO_JOUT
   'share of exports in output (joint production CET function)'{r in
  Regions, i in Sectors}  = RAW_EXO_JOUT[r,i]/(RAW_DOM_JOUT[r,i] +
  RAW_EXO_JOUT[r,i]);
param SHR_DOM_JOUT 'share of domestic uses in output (CET function)'{r in
  Regions, i in Sectors}  = 1 - SHR_EXO_JOUT[r,i];
param PRC_EXO 'exogenous price of exports (in domestic currency units)'
  {i in Sectors, t in LookForward}  default 1;
param EPS_JOUT 'elasticity of subst. for export CET function'{Regions,
  Sectors}  = 4;
param SHR_DOM 'observed share of output for domestic uses'{Regions, Sectors,
  PathTimes}  default 1;
param aA = 0.94;
param c20A = 0.9;
param c10A = 1;
param InstanceName symbolic;
param CH_GROWTH_OUT{Regions, Sectors, PathTimes}  default 0;
var dcon 'consumption flows'{Regions, Sectors, LookForward}  in interval[
  VInf, VSup]
    default 1;
var dinv 'investment flows'{r in Regions, i in Sectors, j in Sectors,
  LookForward}  in interval[VInf, VSup]
    default DELTA[j]*KAP[r,j,PInf];
var dmed 'intermediate flows'{Regions, Sectors, Sectors, LookForward}  in
  interval[VInf, VSup]
    default 1;
var lab 'labour hours'{Regions, Sectors, LookForward}  in interval[VInf,
  LabSup]
    default 0.33;
var lab_ext 'active laborforce but can also be interpreted as effort'{r in
  Regions, i in Sectors, t in LookForward}  in interval[0, 1]
    default NAIRE[r,i,t];
var kap 'kapital stocks (dynamic: defined on LookForwardClosure)'{r in
  Regions, j in Sectors, LookForwardClosure}  in interval[VInf, VSup]
    default KAP[r,j,PInf];
var ccon 'Composite consumption good (Leontief over domestic and imports)'
  {r in Regions, i in Sectors, t in LookForward}  = A_CCON*dcon[r,i,t]/
  SHR_DOM_CCON[r,i];
var ccon_sec_CD
   'Composite consumption aggregate (Cobb-Douglas across sectors)'{r in
  Regions, t in LookForward}  = prod{i in Sectors} ccon[r,i,t]^(SHR_CON[r,i]*
  SCALE_CON);
var con_sec_CD 'Cobb--Douglas consumption aggregate (across sectors)'{r in
  Regions, t in LookForward}  = prod{i in Sectors} dcon[r,i,t]^(SHR_CON[r,i]*
  SCALE_CON);
var con_sec_CES
   'Const. Elast. Subst. consumption aggregate (across sectors)'{r in
  Regions, t in LookForward}  = (sum{i in Sectors} SHR_CON_CES[r,i]*dcon[r,i,t
  ]^RHO_CON)^(RHO_CON_HAT*SCALE_CON);
var con_sec_SumPow 'Sum of power consumption aggregate (across sectors)'
  {r in Regions, t in LookForward}  = sum{i in Sectors} SHR_CON[r,i]*dcon[r,i,
  t]^GAMMA_HAT[r]/GAMMA_HAT[r];
var con_sec_SumShr 'Sum of fractional powers from consumption shares'{r in
  Regions, t in LookForward}  = sum{i in Sectors} dcon[r,i,t]^SHR_CON[r,i];
var prc_dom_CDL
   'domestic price: Cobb-Douglas in sectors Leontief in imports'{j in
  Sectors, t in LookForward}  = SCALE_CON*SHR_CON['GLD',j]*(A_CCON/
  SHR_DOM_CCON['GLD',j])*(A_CON*ccon_sec_CD['GLD',t]/dcon['GLD',j,t]);
var prc_dom 'domestic price / the Lagrange multiplier'{j in Sectors, t in
  LookForward}  = prc_dom_CDL[j,t];
var ymed 'intermediate imports'{r in Regions, i in Sectors, j in Sectors,
  t in LookForward}  = SHR_YSA_CMED[r,i,j]/SHR_DOM_CMED[r,i,j]*(prc_dom[j,t]/
  PRC_YSA[j,t])^EPS_CMED[r,i,j]*dmed[r,i,j,t];
var cmed 'composite intermediate flows'{r in Regions, i in Sectors, j in
  Sectors, t in LookForward}  = A_CMED*(SHR_DOM_CMED[r,i,j]*dmed[r,i,j,t]^
  RHO_CMED[r,i,j] + SHR_YSA_CMED[r,i,j]*ymed[r,i,j,t]^RHO_CMED[r,i,j])^(
  RHO_CMED_HAT[r,i,j]*SCALE_CMED);
var yinv 'investment imports'{r in Regions, i in Sectors, j in Sectors,
  t in LookForward}  = SHR_YSA_CINV[r,i,j]/SHR_DOM_CINV[r,i,j]*(prc_dom[j,t]/
  PRC_YSA[j,t])^EPS_CINV[r,i,j]*dinv[r,i,j,t];
var cinv 'composite investment flows'{r in Regions, i in Sectors, j in
  Sectors, t in LookForward}  = A_CINV*(SHR_DOM_CINV[r,i,j]*dinv[r,i,j,t]^
  RHO_CINV[r,i,j] + SHR_YSA_CINV[r,i,j]*yinv[r,i,j,t]^RHO_CINV[r,i,j])^(
  RHO_CINV_HAT[r,i,j]*SCALE_CINV);
var shr_dom
   'share of output going to domestic uses (p_dom / (p_dom + P_EXO))'{r in
  Regions, i in Sectors, t in LookForward}  = 1/(1 + SHR_EXO_JOUT[r,i]/
  SHR_DOM_JOUT[r,i]*(PRC_EXO[i,t]/prc_dom[i,t])^EPS_JOUT[r,i]);
var inv_sec_CD 'Cobb--Douglas investment aggregate (across sectors)'{r in
  Regions, j in Sectors, t in LookForward}  = prod{i in Sectors} dinv[r,i,j,t]
  ^(SHR_INV_CES[r,i,j]*SCALE_INV);
var inv_sec_CES
   'Const. Elast. Subst. investment aggregate (across sectors)'{r in
  Regions, j in Sectors, t in LookForward}  = (sum{i in Sectors}
  SHR_INV_CES[r,i,j]*dinv[r,i,j,t]^RHO_INV)^(RHO_INV_HAT*SCALE_INV);
var cinv_sec_CES 'composite CES investment aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward}  = (sum{i in Sectors}
  SHR_INV_CES[r,i,j]*cinv[r,i,j,t]^RHO_INV)^(RHO_INV_HAT*SCALE_INV);
var med_sec_CES
   'Const. Elast. Subst. intermediate aggregate (across sectors)'{r in
  Regions, j in Sectors, t in LookForward}  = (sum{i in Sectors}
  SHR_MED_CES[r,i,j]*dmed[r,i,j,t]^RHO_MED)^(RHO_MED_HAT*SCALE_MED);
var cmed_sec_CES
   'Const. Elast. Subst. intermediate aggregate (across sectors)'{r in
  Regions, j in Sectors, t in LookForward}  = (sum{i in Sectors}
  SHR_MED_CES[r,i,j]*cmed[r,i,j,t]^RHO_MED)^(RHO_MED_HAT*SCALE_MED);
var lab_sec_CD 'Cobb--Douglas labour aggregate (across sectors)'{r in
  Regions, t in LookForward}  = prod{j in Sectors} A_LAB[r,t]*lab[r,j,t]^
  SHR_LAB[r,j];
var lab_sec_Q 'quadratic labour aggregate (across sectors)'{r in Regions,
  t in LookForward}  = sum{j in Sectors} lab[r,j,t]^2;
var lab_sec_bangF
   'Frisch labour aggregate (across sectors) flat level sets'{r in Regions,
  t in LookForward}  = (sum{j in Sectors} lab[r,j,t])^RHO_LAB;
var lab_sec_caveF
   'Frisch labour aggregate (across sectors) concave level sets'{r in
  Regions, t in LookForward}  = (sum{j in Sectors} SHR_LAB_CES[r,j]*lab[r,j,t]
  ^RHO_LAB)^(RHO_LAB_HAT*SCALE_LAB);
var lab_ext_sec 'labourforce aggregator'{r in Regions, t in LookForward}  =
  sum{j in Sectors} (NAIRE[r,j,t] - lab_ext[r,j,t])^EXP_LAB_EXT[r,j,t];
var adj_cost_kap_Q 'quadratic adjustment costs for kapital'{r in Regions,
  i in Sectors, t in LookForward}  = PHI_ADJ[i]*kap[r,i,t]*(kap[r,i,t + 1]/
  kap[r,i,t] - 1)^2;
var E_out_CD 'Cobb--Douglas output transformation'{r in Regions, i in
  Sectors, t in LookForward}  = E_shk[r,i,t]*A[i]*(kap[r,i,t]^SHR_KAP_OUT[i]*
  lab[r,i,t]^(1 - SHR_KAP_OUT[i]))^SCALE_OUT;
var E_out_ATA 'Atalay output transformation'{r in Regions, i in Sectors,
  t in LookForward}  = E_shk[r,i,t]*A[i]*(SHR_KAPLAB_OUT_CES[i]*(kap[r,i,t]^
  SHR_KAP_OUT[i]*lab[r,i,t]^SHR_LAB_OUT[i])^RHO_OUT + SHR_MED_OUT_CES[i]*
  med_sec_CES[r,i,t]^RHO_OUT)^(RHO_OUT_HAT*SCALE_OUT);
var E_out_CES 'Constant Elasticity of Substitution output transformation'
  {r in Regions, i in Sectors, t in LookForward}  = E_shk[r,i,t]*A[i]*(
  SHR_KAP_OUT_CES[i]*(kap[r,i,t]*lab_ext[r,i,t])^RHO_OUT + SHR_LAB_OUT_CES[i]*
  (lab[r,i,t]*ALPHA*ALPHA_0^(t + 1))^RHO_OUT + SHR_MED_OUT_CES[i]*(
  med_sec_CES[r,i,t]*lab_ext[r,i,t])^RHO_OUT)^(RHO_OUT_HAT*SCALE_OUT);
var E_cout_CES 'Constant Elasticity of Substitution output transformation'
  {r in Regions, i in Sectors, t in LookForward}  = E_shk[r,i,t]*A[i]*(
  SHR_KAP_OUT_CES[i]*(kap[r,i,t]*lab_ext[r,i,t])^RHO_OUT + SHR_LAB_OUT_CES[i]*
  (lab[r,i,t]*ALPHA*ALPHA_0^(t + 1))^RHO_OUT + SHR_MED_OUT_CES[i]*(
  cmed_sec_CES[r,i,t]*lab_ext[r,i,t])^RHO_OUT)^(RHO_OUT_HAT*SCALE_OUT);
var utility_CD 'Cobb--Douglas instantaneous utility'{t in LookForward}  = sum
  {r in Regions} REG_WGHT[r]*(con_sec_CD[r,t] - lab_sec_CD[r,t]);
var utility_CD_Q 'Cobb--Douglas instantaneous utility'{t in LookForward}  =
  sum{r in Regions} REG_WGHT[r]*(con_sec_CD[r,t] - lab_sec_Q[r,t]);
var utility_pow_CD_Q
   'Power of Cobb-Douglas and Quadratic instantaneous utility'{t in
  LookForward}  = sum{r in Regions} (REG_WGHT[r]*(con_sec_CD[r,t] -
  lab_sec_Q[r,t]))^GAMMA_HAT[r]/GAMMA_HAT[r];
var utility_CES_Q
   'Const. Elast. Subst. and Quadratic instantaneous utility'{t in
  LookForward}  = sum{r in Regions} REG_WGHT[r]*(con_sec_CES[r,t] -
  lab_sec_Q[r,t]);
var utility_powCES_Q
   'Const. Elast. Subst. and Quadratic instantaneous utility'{t in
  LookForward}  = sum{r in Regions} 0.05*REG_WGHT[r]*(con_sec_CES[r,t]^
  GAMMA_HAT[r]/GAMMA_HAT[r] - lab_sec_Q[r,t]);
var utility_SumShr_Q
   'utility: SumShr for consumption, quadratic for labour'{t in LookForward}
   = sum{r in Regions} REG_WGHT[r]*(con_sec_SumShr[r,t] - lab_sec_Q[r,t]);
var utility_SumPow_Q
   'utility: SumPow for consumption and quadratic for labour'{t in
  LookForward}  = sum{r in Regions} REG_WGHT[r]*(con_sec_SumPow[r,t] -
  lab_sec_Q[r,t]);
var utility_CD_F 'Cobb--Douglas and Frisch instantaneous utility'{t in
  LookForward}  = sum{r in Regions} REG_WGHT[r]*(con_sec_CD[r,t] -
  lab_sec_caveF[r,t]);
var utility_CES_bangF
   'Const. Elast. Subst. and conc. Frisch instant. utility'{t in
  LookForward}  = sum{r in Regions} REG_WGHT[r]*(con_sec_CES[r,t] -
  lab_sec_bangF[r,t]);
var utility_CES_caveF
   'Const. Elast. Subst. and conc. Frisch instant. utility'{t in
  LookForward}  = sum{r in Regions} REG_WGHT[r]*(A_CON*con_sec_CES[r,t] +
  A_LAB_EXT*lab_ext_sec[r,t] + A_LAB[r,t]*lab_sec_caveF[r,t]);
var utility_CD_caveF 'Cobb-Douglas and concave Frisch inst. utility'{t in
  LookForward}  = sum{r in Regions} REG_WGHT[r]*(A_CON*con_sec_CD[r,t] +
  A_LAB_EXT*lab_ext_sec[r,t] + A_LAB[r,t]*lab_sec_caveF[r,t]);
var cutility_CD_caveF
   'Cobb-Douglas-Leontief and concave Frisch inst. utility'{t in
  LookForward}  = sum{r in Regions} REG_WGHT[r]*(A_CON*ccon_sec_CD[r,t] +
  A_LAB_EXT*lab_ext_sec[r,t] + A_LAB[r,t]*lab_sec_caveF[r,t]);
var tail_val_CD_F 'continuation value from time LSup + LInf onwards' = (sum
  {r in Regions} REG_WGHT[r]*(prod{i in Sectors} (TAIL_SHR_CON*kap[r,i,LSup +
  LInf]^(SHR_KAP_OUT[i]*SCALE_OUT))^(SHR_CON[r,i]*SCALE_CON) - sum{i in
  Sectors} A_LAB[r,LSup + LInf]*1^RHO_LAB/RHO_LAB))/(1 - BETA);
var tail_val_CD_Q 'SumShr continuation value from time LSup + LInf onwards' =
  (sum{r in Regions} REG_WGHT[r]*(prod{i in Sectors} (TAIL_SHR_CON*kap[r,i,
  LSup + LInf]^(SHR_KAP_OUT[i]*SCALE_OUT))^(SHR_CON[r,i]*SCALE_CON) - 1^2))/(1
   - BETA);
var tail_val_CESutl_Q_CDout
   'continuation value from time LSup + LInf onwards' = (sum{r in Regions}
  REG_WGHT[r]*((sum{i in Sectors} SHR_CON_CES[r,i]*(TAIL_SHR_CON*A[i]*kap[r,i,
  LSup + LInf]^(SHR_KAP_OUT[i]*SCALE_OUT))^RHO_CON)^(RHO_CON_HAT*SCALE_CON) -
  1^2))/(1 - BETA);
var tail_val_CDutl_F_CESout
   'continuation value from time LSup + LInf onwards' = (sum{r in Regions}
  REG_WGHT[r]*(prod{i in Sectors} (TAIL_SHR_CON*A[i]*(SHR_KAP_OUT_CES[i]*
  kap[r,i,LSup + LInf]^RHO_OUT + SHR_MED_OUT_CES[i]*1^RHO_OUT +
  SHR_LAB_OUT_CES[i]*1^RHO_OUT)^(RHO_OUT_HAT*SCALE_OUT))^(SHR_CON[r,i]*
  SCALE_CON) - sum{i in Sectors} A_LAB[r,LSup + LInf]*1^RHO_LAB/RHO_LAB))/(1
   - BETA);
var tail_val_CESutl_bangF_CESout
   'continuation value from time LSup + LInf onwards' = (sum{r in Regions}
  REG_WGHT[r]*(A_CON*(sum{i in Sectors} SHR_CON_CES[r,i]*(TAIL_SHR_CON*A[i]*(
  SHR_KAP_OUT_CES[i]*kap[r,i,LSup + LInf]^RHO_OUT + SHR_MED_OUT_CES[i]*1^
  RHO_OUT + SHR_LAB_OUT_CES[i]*1^RHO_OUT)^(RHO_OUT_HAT*SCALE_OUT))^RHO_CON)^(
  RHO_CON_HAT*SCALE_CON) - A_LAB[r,LSup + LInf]*(sum{i in Sectors} 1)^
  RHO_LAB))/(1 - BETA);
var tail_val_CESutl_caveF_CESout
   'continuation value from time LSup + LInf onwards' = (sum{r in Regions}
  REG_WGHT[r]*(A_CON*(sum{i in Sectors} SHR_CON_CES[r,i]*(TAIL_SHR_CON*A[i]*(
  SHR_KAP_OUT_CES[i]*kap[r,i,LSup + LInf]^RHO_OUT + SHR_MED_OUT_CES[i]*1^
  RHO_OUT + SHR_LAB_OUT_CES[i]*(1*ALPHA*ALPHA_0^(LSup + LInf))^RHO_OUT)^(
  RHO_OUT_HAT*SCALE_OUT))^RHO_CON)^(RHO_CON_HAT*SCALE_CON) + A_LAB_EXT*0 +
  A_LAB[r,LSup + LInf]*(sum{i in Sectors} SHR_LAB_CES[r,i]*0.33^RHO_LAB)^(
  RHO_LAB_HAT*SCALE_LAB)))/(1 - BETA);
var tail_val_CDutl_caveF_CESout
   'continuation value from time LSup + LInf onwards' = (sum{r in Regions}
  REG_WGHT[r]*(A_CON*(prod{i in Sectors} (TAIL_SHR_CON*A[i]*(
  SHR_KAP_OUT_CES[i]*kap[r,i,LSup + LInf]^RHO_OUT + SHR_MED_OUT_CES[i]*1^
  RHO_OUT + SHR_LAB_OUT_CES[i]*(1*ALPHA*ALPHA_0^(LSup + LInf))^RHO_OUT)^(
  RHO_OUT_HAT*SCALE_OUT))^(SHR_CON[r,i]*SCALE_CON)) + A_LAB_EXT*0 + A_LAB[r,
  LSup + LInf]*(sum{i in Sectors} SHR_LAB_CES[r,i]*0.33^RHO_LAB)^(
  RHO_LAB_HAT*SCALE_LAB)))/(1 - BETA);
var con 'consumption'{r in Regions, i in Sectors, t in LookForward}  =
  ccon[r,i,t];
var inv 'domestic investment flows'{r in Regions, i in Sectors, j in
  Sectors, t in LookForward}  = dinv[r,i,j,t];
var med 'domestic intermedate flows'{r in Regions, i in Sectors, j in
  Sectors, t in LookForward}  = cmed[r,i,j,t];
var inv_sec 'current intermediate variable for aggregated investment'{r in
  Regions, j in Sectors, t in LookForward}  = A_INV*inv_sec_CES[r,j,t];
var E_out 'current intermediate variable for output'{r in Regions, i in
  Sectors, t in LookForward}  = E_cout_CES[r,i,t];
var dom 'domestic uses'{r in Regions, i in Sectors, t in LookForward}  =
  shr_dom[r,i,t]*E_out[r,i,t];
var exo 'exports'{r in Regions, i in Sectors, t in LookForward}  = (1 -
  shr_dom[r,i,t])*E_out[r,i,t];
var utility 'current intermediate variable for utility'{t in LookForward}  =
  cutility_CD_caveF[t];
var adj_cost_kap 'current adjustment costs for kapital'{r in Regions, i in
  Sectors, t in LookForward}  = adj_cost_kap_Q[r,i,t];
var tail_val 'current intermediate variable for tail value function' =
  A_VAL*tail_val_CDutl_caveF_CESout;
subject to kap_transition 'equation for the accumulation of kapital'{r in
  Regions, j in Sectors, t in LookForward} : kap[r,j,t + 1] - (1 - DELTA[j])*
  kap[r,j,t] - inv_sec[r,j,t] == 0;
subject to market_clearing 'market clearing for each sector and time'{i in
  Sectors, t in LookForward} : sum{r in Regions} (con[r,i,t] + sum{j in
  Sectors} inv[r,i,j,t] + sum{j in Sectors} med[r,i,j,t] + adj_cost_kap[r,i,t]
   - dom[r,i,t]) == 0;
maximize pres_disc_val 'present discounted value of utilities': sum{t in
  LookForward} BETA^(t - LInf)*utility[t] + BETA^(LSup - LInf)*tail_val;
