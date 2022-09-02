/*=============================================================================
This is NAIU (New Australian Inter* Model with Uncertainty, where * stands for
Sectoral, Time, Age and Regional). It is an extension of
Cai--Judd (2021) to include multiple sectors among among other things. If any
part of this code is re-used, please cite OCallaghan (2022).
=============================================================================*/

#==============================================================================
# THE MODEL starts here
#==============================================================================
#==============================================================================
# Parameters for generating sets
#==============================================================================
param TInf 'infimum of the set of times' default 0;
param LInf 'start period/infimum of the look forward set', default 0;
param LSup 'supremum of the look forward set' > LInf, default 48e+0; 
param PInf 'infimum of times on a path', default 0;
param PSup 'supremum of times on a path (T_star in CJ)'# eg 2050 - 2022 = 28 
  default 3 >= PInf; 
param TSup 'final time' >= PInf + LInf, default PSup + LSup;
#==============================================================================
# Sets
#============================================================================== 
set Regions; 
set Sectors ordered;
#-----------the planning horizon and corresponding set:
#-----------lower values of LFwd are more myopic; longer satisfy Euler Eqn.
set LookForward 'planning set' = {LInf .. LSup + LInf - 1} ordered;
set LookForwardClosure 'includes extra time for closing out the problem'
  = {LInf .. LSup + LInf} ordered;
#-----------Each path is a state of the world and has the following structure:
set PathTimes 'times along a path' = {PInf .. PSup + PInf - 1} ordered;
set PathTimesClosure 'includes extra time for closing out the problem'
  = {PInf .. PSup + PInf} ordered;
#-----------The set of paths are indexed by
set PathSpace 'path space: the set of paths'
  # for default we adopt a two-state Markov chain with unique absorbing state
    default {PInf .. PSup + 1};  
set AllTimes 'all time periods associated with a path' = {TInf .. TSup};                      
#==============================================================================
#-----------Parameters
#==============================================================================
param KAP 'observed kapital' {Regions, Sectors, PathTimesClosure}
  default 1e+0; # in (OInf, OSup);
param StrtTm 'start time for each step on each path' {PathTimes, PathSpace}
  default time();
param EndTm 'start time for each step on each path' {PathTimes, PathSpace}
  default time();
param ALPHA 'trend growth factor', >= 1, default 100e-2; 
param ALPHA_0 'initial trend growth factor', >= 1, default 100e-2; 
param BETA 'discount factor', in (0, 1), default .975; 
param DELTA 'rate of depreciation for kapital' {Sectors} default .025 >= 0;
param PHI_ADJ 'kapital adjustment cost' {Sectors} default 5e-1 >= 0;
param GAMMA 'intertemporal elasticity of subst.' {Regions} default 5e-1 >= 0;
param EPS_LAB 'Frisch elasticity of labour supply' default 5e-1 >= 0;
param EPS_CON 'elasticity of substitution for consumption flows'
  default 1e-1 > 0; 
param EPS_OUT 'elasticity of substitution for output components'
  default 5e-1 > 0; 
param EPS_MED 'elasticity of substitution for intermediate flows'
  default 1e-1 > 0; 
param EPS_INV 'elasticity of substitution for kapital flows'
  default 1e-1 > 0; 
param SCALE_INV 'investment scale parameter' default 90e-2; # in (0, 1);
param SCALE_CON 'consumption scale parameter' default 90e-2; # in (0, 1);
param SCALE_OUT 'output scale parameter' default 90e-2; # in (0, 1);
param SCALE_MED 'intermediate scale parameter' default 90e-2; # in (0, 1);
param SCALE_LAB 'labour scale parameter' default 300e-2; # in (0, 1);

param INV_MIN 'lower bound of investment' default 0 >= 0;
#param kmin 'smallest capital' default 1e-1 >= 0;
#param kmax 'largest capital' default 1e+1 >= 0;
param ZETA1 'TFP before shock' default 1 >= 0;
#param ZETA2 'TFP after shock' default .95 >= 0;
param ZETA2 'TFP after shock' default 1 >= 0;
#param ZETA2 'TFP after shock' default __ZETA2__ >= 0;
param PROB2 'one period probability of jump in TFP' default 0.01 >= 0;
param TAIL_SHR_CON 'tail consumption share (of output)' default 0.45 >= 0;
#
param UInf 'infimum of interval for uniform dbn' default 0.4999 in [0, .5);
param USup 'supremum of interval for uniform dbn' default 0.5001 in (.5, 1];
param VInf 'infimum of interval for basic variables' default 1e-4;
param VSup 'supremum of interval for basic variables' default 1e+6;
param OInf 'infimum of interval for observed/actual values' default 1e-7;
param OSup 'supremum of interval for observed/actual values' default 1e+7;
param LabSup 'supremum of interval for labour values' default 66e-2;
param EXP_LAB_EXT 'Exponent of lab_ext_sec'
  {Regions, Sectors, PathTimes} default 200e-2;
param NAIRE 'Non-accelerating rate of employment'
  {Regions, Sectors, PathTimes} default 100e-2;
#-----------productivity and relative importance of labour in utility
param A 'productivity trend'
  {i in Sectors}
  default 1; #(1 - (1 - DELTA[i]) * BETA) / (SHR_KAP_OUT_CES[r, i] * BETA) >= 0;
param A_LAB 'importance of disutility of labour (weight in utility function)' 
  {Regions, LookForwardClosure}
    default -1;
param A_LAB_EXT 'disutility weight for labourforce deviations in utility'
  default -1;
  #(1 - SHR_KAP_OUT_CES[r, i]) * A[i] * (A[i] - DELTA[i]) ** (-1 / GAMMA[r]) >= 0;
param A_CON 'importance of consumption in utility'
  default 1;
param A_VAL 'Calibration factor for terminal value function'
  default 1;
param A_INV 'Calibration factor for investment'
  default 1;
param A_MED 'Calibration factor for intermediate bundle'
  default 1;
#==============================================================================
# raw flow data parameters
#==============================================================================
#-----------set the seed for the random number generator for weights
option randseed 12345;
param RAW_CON_FLW "raw consumption flows: table8Q1"
  {Regions, Sectors}
  default Uniform(UInf, USup) >= 0;
param RAW_INV_FLW "raw investment flows: tablekapflw"
  {Regions, Sectors, Sectors} 
  default Uniform(UInf, USup) >= 0;
param REG_LAB  "raw number of workers in region"
  {Regions, Sectors}
  default Uniform(UInf, USup) >= 0;
#param REF_LAB  "raw number of workers in reference region(s)"
#  {RefRegions, Sectors}
#  default Uniform(UInf, USup) >= 0;
param RAW_MED_FLW "raw intermediate flows: table8"
  {Regions, Sectors, Sectors} 
  default Uniform(UInf, USup) >= 0;
#-----------raw shares of the three components of output
param RAW_KAP_OUT 'raw capital primary input flows: table8P2'
  {Regions, Sectors}
  default Uniform(UInf, USup) >= 0;
param RAW_LAB_OUT 'raw compensation of employees primary input flows: table8P1'
  {Regions, Sectors}
  default Uniform(UInf, USup) >= 0;
param RAW_MED_OUT 'raw total intermediate input flows: table8T1'
  {Regions, Sectors}
  default Uniform(UInf, USup) >= 0;
#-----------raw Armington data
param RAW_DOM_CCON 'raw domestic flows to consumption: table5Q1'
  {Regions, Sectors}
  default Uniform(UInf, USup) >= 0;
param RAW_YPO_CCON 'raw import flows to consumption: table(8-5)'
  {r in Regions, i in Sectors}
  default Uniform(UInf, USup) >= 0;
param RAW_DOM_CINV 'raw domestic flows to investment: kapflows per table5'
  {r in Regions, i in Sectors, j in Sectors}
  default Uniform(UInf, USup) >= 0;
param RAW_YPO_CINV 'raw import flows to investment: kapflows per tablediff'
  {r in Regions, i in Sectors, j in Sectors}
  default Uniform(UInf, USup) >= 0;
param RAW_DOM_CMED 'raw domestic flows to intermediates: table5'
  {r in Regions, i in Sectors, j in Sectors}
  default Uniform(UInf, USup) >= 0;
param RAW_YPO_CMED 'raw import flows to intermediates: table(8-5)'
  {r in Regions, i in Sectors, j in Sectors}
  default Uniform(UInf, USup) >= 0;
param RAW_XPO_JOUT 'raw export data: table8Q7'
  {Regions, Sectors}
  default Uniform(UInf, USup) * 10e-2;
param RAW_DOM_JOUT 'raw total domestic uses: table8T6-Q7'
  {Regions, Sectors}
  default Uniform(UInf, USup) * 90e-2;
param SHR_MED_ROW 'intermediate flows as a share of total demand (quasi-raw)'
  {Regions, Sectors, Sectors}
  default 100e-2;
param SHR_MED_COL 'intermediate flows as a share of total cost (quasi-raw)'
  {Regions, Sectors, Sectors}
  default 100e-2;
param RAW_OUT_REG_SEC 'raw output per region and sector'
  {Regions, Sectors}
  default 1;
#==============================================================================
# Computed parameters
#==============================================================================
param GAMMA_HAT 'utility parameter' {r in Regions} = 1 - 1 / GAMMA[r];
param RHO_LAB 'labour exponent parameter' = 1 + 1 / EPS_LAB;
param RHO_LAB_HAT 'inverse of labour exponent parameter'
  = 1 / RHO_LAB;
param RHO_INV 'exponent of the investment ces aggregator'
  = 1 - 1 / EPS_INV; 
param RHO_INV_HAT 'inverse of RHO_INV', = 1 / RHO_INV;
param RHO_MED 'exponent of the intermediate ces aggregator'
  = 1 - 1 / EPS_MED; 
param RHO_MED_HAT 'inverse of RHO_MED', = 1 / RHO_MED;
param RHO_OUT 'exponent of the output ces aggregator'
  = 1 - 1 / EPS_OUT; 
param RHO_OUT_HAT 'inverse of RHO_OUT', = 1 / RHO_OUT;
param RHO_CON 'exponent of the CON ces aggregator'
  = 1 - 1 / EPS_CON; 
param RHO_CON_HAT 'inverse of RHO_CON', = 1 / RHO_OUT;
param  REG_WGHT 'regional (population) weights' {r in Regions}
  default 1 / card(Regions);
#-----------Share parameters in the aggregators for utility and production.
#-----------For default values, we draw from the uniform distribution.
#-----------input shares for output
param SHR_LAB_OUT 'share of labour in output'
  {r in Regions, i in Sectors}
  = RAW_LAB_OUT[r, i] / (RAW_KAP_OUT[r, i] + RAW_LAB_OUT[r, i] + RAW_MED_OUT[r, i]); 
param SHR_KAP_OUT 'importance of capital in production'
  {r in Regions, i in Sectors}
  = RAW_KAP_OUT[r, i] / (RAW_KAP_OUT[r, i] + RAW_LAB_OUT[r, i] + RAW_MED_OUT[r, i]);
param SHR_KAPLAB_OUT 'combined importance of kapital and labour in output'
  {r in Regions, i in Sectors}
  = SHR_KAP_OUT[r, i] + SHR_LAB_OUT[r, i];
param SHR_MED_OUT 'share of intermediates in output'
  {r in Regions, i in Sectors}
  = 1 - SHR_KAPLAB_OUT[r, i];
#-----------output component shares for CES functions
param SHR_KAP_OUT_CES 'importance of kapital in production'
  {r in Regions, i in Sectors}
  = SHR_KAP_OUT[r, i] ** (1 / EPS_OUT);
param SHR_LAB_OUT_CES 'importance of labour in production'
  {r in Regions, i in Sectors}
  = SHR_LAB_OUT[r, i] ** (1 / EPS_OUT);
param SHR_KAPLAB_OUT_CES 'combined importance of kapital and labour in prod'
  {r in Regions, i in Sectors}
  = SHR_KAPLAB_OUT[r, i] ** (1 / EPS_OUT);
param SHR_MED_OUT_CES 'importance of intermediates in production'
  {r in Regions, i in Sectors}
  = SHR_MED_OUT[r, i] ** (1 / EPS_OUT);
#-----------consumption
param CON_FLW_SUM
  {r in Regions}
  = sum{i in Sectors} RAW_CON_FLW[r, i];
param SHR_CON 'consumption weights for each good in utility (for Cobb-Doug)'
  {r in Regions, i in Sectors} = (RAW_CON_FLW[r, i] / CON_FLW_SUM[r]);
  #{r in Regions, i in Sectors} = RAW_CON_FLW[r, i];
param SHR_CON_CES 'CES consumption weights for each good in utility'
  {r in Regions, i in Sectors}
  = (RAW_CON_FLW[r, i] / CON_FLW_SUM[r]) ** (1 / EPS_CON);
#-----------labour
param LAB_FLW_SUM {r in Regions}
  = sum{i in Sectors} REG_LAB[r, i];
param SHR_LAB 'labour weights for each sector in utility'
  {r in Regions, j in Sectors}
  = REG_LAB[r, j] / LAB_FLW_SUM[r];
  #{r in Regions, j in Sectors} = REG_LAB[r, j];
#-----------inputs of sector i into j's investment).
param INV_FLW_RSUM {r in Regions, j in Sectors} 
  = sum{i in Sectors} RAW_INV_FLW[r, i, j];
param SHR_INV_CES "sectoral share of i in j's CES investment aggregator"
  {r in Regions, i in Sectors, j in Sectors}
    = (RAW_INV_FLW[r, i, j] / INV_FLW_RSUM[r, j]) ** (1 / EPS_INV);
    #= RAW_INV_FLW[r, i, j] ** (1 / EPS_INV);
#-----------inputs of sector i into j's intermediate bundle).
param MED_FLW_RSUM {r in Regions, j in Sectors} 
  = sum{i in Sectors} RAW_MED_FLW[r, i, j];
#  = 1;
param SHR_MED_CES "sectoral share of i in j's CES intermediate aggregator"
  {r in Regions, i in Sectors, j in Sectors}
    = (RAW_MED_FLW[r, i, j] / MED_FLW_RSUM[r, j]) ** (1 / EPS_INV);
    #= RAW_MED_FLW[r, i, j] ** (1 / EPS_INV);
param SHR_LAB_CES "sectoral share of i in the CES labour aggregator"
  {r in Regions, i in Sectors}
    = SHR_LAB[r, i];
#-----------Armington share parameters
param SHR_DOM_CCON 'domestic share in comp. consumption'
  {r in Regions, i in Sectors}
  = RAW_DOM_CCON[r, i] / (RAW_DOM_CCON[r, i] + RAW_YPO_CCON[r, i]);
param SHR_YPO_CCON 'domestic share in comp. consumption'
  {r in Regions, i in Sectors}
  = 1 - SHR_DOM_CCON[r, i];
param SHR_DOM_CINV 'domestic share in composite investment flows'
  {r in Regions, i in Sectors, j in Sectors}
  = RAW_DOM_CINV[r, i, j] / (RAW_DOM_CINV[r, i, j] + RAW_YPO_CINV[r, i, j]);
param SHR_YPO_CINV 'import share in composite investment flows'
  {r in Regions, i in Sectors, j in Sectors}
  = 1 - SHR_DOM_CINV[r, i, j];
param SHR_DOM_CMED 'domestic share in composite intermediate flows'
  {r in Regions, i in Sectors, j in Sectors}
  = RAW_DOM_CMED[r, i, j] / (RAW_DOM_CMED[r, i, j] + RAW_YPO_CMED[r, i, j]);
param SHR_YPO_CMED 'import share in composite intermediate flows'
  {r in Regions, i in Sectors, j in Sectors}
  = 1 - SHR_DOM_CMED[r, i, j];
/*-----------------------------------------------------------------------------
uncertainty parameters
-----------------------------------------------------------------------------*/
param Pr_shk 'probability of SHK'
  {Regions, Sectors, t in LookForward} = (1 - PROB2) ** t;
param E_shk 'expected shock (exogenous)'
  {r in Regions, i in Sectors, t in LookForward}
    = ZETA2 + Pr_shk[r, i, t] * (ZETA1 - ZETA2);
/*-----------------------------------------------------------------------------
Computed Armington parameters
-----------------------------------------------------------------------------*/
#-----------foreign prices
param PRC_YPO 'import prices' {Sectors, LookForward} default 100e-2;
param PRC_XPO 'export prices' {Sectors, LookForward} default 100e-2;
#-----------calibration factors
param A_CCON 'calibration factor for composite consumption'
  default 100e-2;
param A_CMED 'scaling factor for composite intermediates'
  default 100e-2;
param A_CINV 'scaling factor for composite intermediates'
  default 100e-2;
#-----------economies of scale exponents
param SCALE_CINV 'economies of scale for composite intermediate'
  default 100e-2;
param SCALE_CMED 'economies of scale for composite intermediate'
  default 100e-2;
#-----------elasticities
param EPS_CINV 'elasticity of subst. for composite investment flows'
  {Regions, Sectors, Sectors}
  default 0200e-2;
param EPS_CMED 'elasticity of subst. for composite intermediates'
  {Regions, Sectors, Sectors}
  default 0200e-2;
#-----------shorthand elasticity parameters
param RHO_CINV 'CES exponent for composite intermediates'
  {r in Regions, i in Sectors, j in Sectors}
  = 1 - 1 / EPS_CINV[r, i, j];
param RHO_CINV_HAT 'inverse of CES exponent for composite intermediates'
  {r in Regions, i in Sectors, j in Sectors}
  = 1 / RHO_CINV[r, i, j];
param RHO_CMED 'CES exponent for composite intermediates'
  {r in Regions, i in Sectors, j in Sectors}
  = 1 - 1 / EPS_CMED[r, i, j];
param RHO_CMED_HAT 'inverse of CES exponent for composite intermediates'
  {r in Regions, i in Sectors, j in Sectors}
  = 1 / RHO_CMED[r, i, j];
/*=============================================================================
===============================================================================
The variables
===============================================================================
===============================================================================
Basic (economic) variables
=============================================================================*/
var dcon 'consumption flows' {Regions, Sectors, LookForward}
  in [VInf, VSup] default 1e-0;
var dinv 'investment flows'
  {r in Regions, i in Sectors, j in Sectors, LookForward}
  in [VInf, VSup] default DELTA[j] * KAP[r, j, PInf];
var dmed 'intermediate flows'
  {Regions, Sectors, Sectors, LookForward}
    in [VInf, VSup] default 1e-0;
/*-----------------------------------------------------------------------------
labour split into time spent working and laborforce both normalised
-----------------------------------------------------------------------------*/
var lab 'labour hours' {Regions, Sectors, LookForward}
  in [VInf, LabSup] default 33e-2;
#-----------to do: add occupations as a set in the following
var lab_ext 'active laborforce but can also be interpreted as effort'
  {r in Regions, i in Sectors, t in LookForward}
  in [0, 100e-2] default NAIRE[r, i, t];
/*-----------------------------------------------------------------------------
kapital, the dynamic variable, is defined on LookForwardClosure
-----------------------------------------------------------------------------*/
var kap 'kapital stocks (dynamic: defined on LookForwardClosure)'
  {r in Regions, j in Sectors, LookForwardClosure}
  in [VInf, VSup] default KAP[r, j, PInf]; 
/*=============================================================================
Computed variables for paths
=============================================================================*/
#var inv_sum 'sum over destinations for saving'
# {r in Regions, i in Sectors, s in PathTimes}
#   = sum{j in Sectors} inv[r, i, j, LFwd, s];
#var kap_growth {r in Regions, i in Sectors, s in PathTimes}
#  = (kap[r, i, LFwd + 1] - KAP[r, i, LFwd]) / KAP[r, i, s];
##-----------Euler integrand for CES production
#    let EULER_INTEGRAND[r, i, s] := DUAL_KAP[r, i, s] * (1 - DELTA[i]) 
#      + DUAL_MKT_CLR[i, s] * (
#        SHR_KAP_OUT_CES[r, i] * (KAP[r, i, s] / E_OUT[r, i, s]) ** (RHO_INV - 1)
#        - PHI_ADJ[i] * (2 * GROWTH_KAP[r, i, s] + GROWTH_KAP[r, i, s] ** 2)
#      );
#    if s > PInf then 
#      let EULER_RATIO[r, i, s] 
#        := BETA * EULER_INTEGRAND[r, i, s] / DUAL_KAP[r, i, s - 1];
#  };
/*=============================================================================
Potential intermediate variables (substituted out during pre-solving)
=============================================================================*/
#-----------variety of consumption aggregator functions 
var ccon 'Composite consumption good (Leontief over domestic and imports)'
  {r in Regions, i in Sectors, t in LookForward}
  = A_CCON * dcon[r, i, t] / SHR_DOM_CCON[r, i];
var ccon_sec_CD 'Composite consumption aggregate (Cobb-Douglas across sectors)'
  {r in Regions, t in LookForward}
  = prod{i in Sectors} ccon[r, i, t] ** (SHR_CON[r, i] * SCALE_CON);
var con_sec_CD 'Cobb--Douglas consumption aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = prod{i in Sectors} dcon[r, i, t] ** (SHR_CON[r, i] * SCALE_CON);
var con_sec_CES 'Const. Elast. Subst. consumption aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = (sum{i in Sectors} SHR_CON_CES[r, i] * dcon[r, i, t] ** RHO_CON)
    ** (RHO_CON_HAT * SCALE_CON); 
var con_sec_SumPow 'Sum of power consumption aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = sum{i in Sectors}
    SHR_CON[r, i] * dcon[r, i, t] ** GAMMA_HAT[r] / GAMMA_HAT[r];
var con_sec_SumShr 'Sum of fractional powers from consumption shares'
  {r in Regions, t in LookForward}
  = sum{i in Sectors} dcon[r, i, t] ** SHR_CON[r, i];
/*-----------------------------------------------------------------------------
Armington variables
-----------------------------------------------------------------------------*/
#-----------domestic prices
var prc_dom_CDL 'domestic price: Cobb-Douglas in sectors Leontief in imports'
  {j in Sectors, t in LookForward}
  = SCALE_CON * SHR_CON['GLD', j] * (A_CCON / SHR_DOM_CCON['GLD', j])
    * (A_CON * ccon_sec_CD['GLD', t] / dcon['GLD', j, t]);
var prc_dom  'domestic price / the Lagrange multiplier'
  {j in Sectors, t in LookForward}
  = prc_dom_CDL[j, t];
#-----------imports
var ymed 'intermediate imports'
  {r in Regions, i in Sectors, j in Sectors, t in LookForward}
  = (SHR_YPO_CMED[r, i, j] / SHR_DOM_CMED[r, i, j])
    * (prc_dom[j, t] / PRC_YPO[j, t]) ** EPS_CMED[r, i, j]
    * dmed[r, i, j, t];
var cmed 'composite intermediate flows'
  {r in Regions, i in Sectors, j in Sectors, t in LookForward}
  = A_CMED * (
    SHR_DOM_CMED[r, i, j] * dmed[r, i, j, t] ** RHO_CMED[r, i, j]
    + SHR_YPO_CMED[r, i, j] * ymed[r, i, j, t] ** RHO_CMED[r, i, j]
  ) ** (RHO_CMED_HAT[r, i, j] * SCALE_CMED);
var yinv 'investment imports'
  {r in Regions, i in Sectors, j in Sectors, t in LookForward}
  = (SHR_YPO_CINV[r, i, j] / SHR_DOM_CINV[r, i, j])
    * (prc_dom[j, t] / PRC_YPO[j, t]) ** EPS_CINV[r, i, j]
    * dinv[r, i, j, t];
var cinv 'composite investment flows' 
  {r in Regions, i in Sectors, j in Sectors, t in LookForward}
  = A_CINV * (
    SHR_DOM_CINV[r, i, j] * dinv[r, i, j, t] ** RHO_CINV[r, i, j]
    + SHR_YPO_CINV[r, i, j] * yinv[r, i, j, t] ** RHO_CINV[r, i, j]
  ) ** (RHO_CINV_HAT[r, i, j] * SCALE_CINV);
#-----------exports (other variables appear after output)
param SHR_XPO_JOUT 'share of exports in output (joint production CET function)'
  {r in Regions, i in Sectors}
  = RAW_XPO_JOUT[r, i] / (RAW_DOM_JOUT[r, i] + RAW_XPO_JOUT[r, i]);
param SHR_DOM_JOUT 'share of domestic uses in output (CET function)'
  {r in Regions, i in Sectors}
  = 1 - SHR_XPO_JOUT[r, i];
param EPS_JOUT 'elasticity of subst. for export CET function'
  {Regions, Sectors}
  default 4;
param SHR_DOM 'observed share of output for domestic uses'
  {Regions, Sectors, PathTimes}
  default 100e-2;
var shr_dom 'share of output going to domestic uses (p_dom / (p_dom + P_XPO))'
  {r in Regions, i in Sectors, t in LookForward}
    = 1 / (
      1 + SHR_XPO_JOUT[r, i] / SHR_DOM_JOUT[r, i]
        * (PRC_XPO[i, t] / prc_dom[i, t]) ** EPS_JOUT[r, i]
    );
#var dom 'domestic production'
#  = (A_YPO[j] * SHR_DOM_OUT[j] * mkt_clr[j, t] / prc_dom[j, t])
#    ** (1 - 1 / RHO_YPO[j])
#    * out[r, j, t];
#var gdp 'gross domestic product'
#  {r in Regions, j in Sectors, t in LookForward}
#  = (
#    (A_XPO[j] ** RHO_XPO[j] * SHR_DOM_GDP[j] * (1 + TAX_GDP[j]) * mkt_clr[j, t])
#      / (A_YPO[j] ** RHO_YPO[j] * SHR_YPO_OUT[j] * mkt_clr[j, t])
#  ) ** (1  - 1 / RHO_YPO[j])
#    * out[r, j, t];
#-----------variety of investment aggregator functions 
var inv_sec_CD 'Cobb--Douglas investment aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward}
  = prod{i in Sectors}
    dinv[r, i, j, t] ** (SHR_INV_CES[r, i, j] * SCALE_INV);
var inv_sec_CES 'Const. Elast. Subst. investment aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward}
  = (sum{i in Sectors} SHR_INV_CES[r, i, j] * dinv[r, i, j, t] ** RHO_INV)
    ** (RHO_INV_HAT * SCALE_INV);
var cinv_sec_CES 'composite CES investment aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward}
  = (sum{i in Sectors} SHR_INV_CES[r, i, j] * cinv[r, i, j, t] ** RHO_INV)
    ** (RHO_INV_HAT * SCALE_INV);
#-----------variety of intermediate aggregator functions 
var med_sec_CES 'Const. Elast. Subst. intermediate aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward}
  = (sum{i in Sectors} SHR_MED_CES[r, i, j] * dmed[r, i, j, t] ** RHO_MED)
    ** (RHO_MED_HAT * SCALE_MED);
var cmed_sec_CES 'Const. Elast. Subst. intermediate aggregate (across sectors)'
  {r in Regions, j in Sectors, t in LookForward}
  = (sum{i in Sectors} SHR_MED_CES[r, i, j] * cmed[r, i, j, t] ** RHO_MED)
    ** (RHO_MED_HAT * SCALE_MED);
#-----------variety of labour aggregator functions 
var lab_sec_CD 'Cobb--Douglas labour aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = prod{j in Sectors} A_LAB[r, t] * lab[r, j , t] ** SHR_LAB[r, j];
var lab_sec_Q 'quadratic labour aggregate (across sectors)'
  {r in Regions, t in LookForward}
  = sum{j in Sectors} (lab[r, j , t] ** 2);
var lab_sec_bangF 'Frisch labour aggregate (across sectors) flat level sets'
  {r in Regions, t in LookForward}
  #= (sum{j in Sectors} lab[r, j , t]) ** RHO_LAB / RHO_LAB;
  = (sum{j in Sectors} lab[r, j , t]) ** RHO_LAB;
var lab_sec_caveF 'Frisch labour aggregate (across sectors) concave level sets'
  {r in Regions, t in LookForward}
  = (sum{j in Sectors} SHR_LAB_CES[r, j] * lab[r, j , t] ** RHO_LAB)
    ** (RHO_LAB_HAT * SCALE_LAB);
var lab_ext_sec 'labourforce aggregator'
  {r in Regions, t in LookForward}
  = sum{j in Sectors}
    (NAIRE[r, j, t] - lab_ext[r, j, t]) ** EXP_LAB_EXT[r, j, t];
#-----------variety of adjustment cost functions 
var adj_cost_kap_Q 'quadratic adjustment costs for kapital'
  {r in Regions, i in Sectors, t in LookForward}
  = PHI_ADJ[i] * kap[r, i, t]
      * (kap[r, i, t + 1] / kap[r, i, t] - 1) ** 2;
#-----------variety of production functions 
#var E_out_CD 'Cobb--Douglas output transformation'
#  {r in Regions, i in Sectors, t in LookForward}
#  = E_shk[r, i, t] * A[i]
#      * (kap[r, i, t] ** SHR_KAP_OUT[r, i] * lab[r, i, t] ** (1 - SHR_KAP_OUT[r, i]))
#        ** SCALE_OUT;
#var E_out_ATA 'Atalay output transformation'
#  {r in Regions, i in Sectors, t in LookForward}
#  = E_shk[r, i, t] * A[i] * (
#    SHR_KAPLAB_OUT_CES[r, i]
#      * (kap[r, i, t] ** SHR_KAP_OUT[r, i] * lab[r, i, t] ** SHR_LAB_OUT[r, i])
#        ** RHO_OUT
#    + SHR_MED_OUT_CES[r, i] * med_sec_CES[r, i, t] ** RHO_OUT
#    ) ** (RHO_OUT_HAT * SCALE_OUT);
var E_out_CES 'Constant Elasticity of Substitution output transformation'
  {r in Regions, i in Sectors, t in LookForward}
  = E_shk[r, i, t] * A[i] * (
    SHR_KAP_OUT_CES[r, i]
      * (kap[r, i, t] * lab_ext[r, i, t]) ** RHO_OUT
    + SHR_LAB_OUT_CES[r, i]
      * (lab[r, i, t] ) ** RHO_OUT 
    + SHR_MED_OUT_CES[r, i]
      * (med_sec_CES[r, i, t] * lab_ext[r, i, t]) ** RHO_OUT
    + 1
      * (33e-2 * ALPHA * ALPHA_0 ** t) ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT);
var E_cout_CES 'Constant Elasticity of Substitution output transformation'
  {r in Regions, i in Sectors, t in LookForward}
  = E_shk[r, i, t] * A[i] * (
    SHR_KAP_OUT_CES[r, i]
      * (kap[r, i, t] * lab_ext[r, i, t]) ** RHO_OUT
    + SHR_LAB_OUT_CES[r, i]
      * (lab[r, i, t] ) ** RHO_OUT 
    + SHR_MED_OUT_CES[r, i]
      * (cmed_sec_CES[r, i, t] * lab_ext[r, i, t]) ** RHO_OUT
    + 1
      * (33e-2 * ALPHA * ALPHA_0 ** t) ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT);
#-----------variety of utility functions
var utility_CD 'Cobb--Douglas instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CD[r, t] - lab_sec_CD[r, t]);
var utility_CD_Q 'Cobb--Douglas instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CD[r, t] - lab_sec_Q[r, t]);
var utility_pow_CD_Q 'Power of Cobb-Douglas and Quadratic instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      (REG_WGHT[r] * (con_sec_CD[r, t] - lab_sec_Q[r, t]))
        ** GAMMA_HAT[r] / GAMMA_HAT[r];
var utility_CES_Q 'Const. Elast. Subst. and Quadratic instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CES[r, t] - lab_sec_Q[r, t]);
var utility_powCES_Q 'Const. Elast. Subst. and Quadratic instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      5e-2 * REG_WGHT[r] * ((con_sec_CES[r, t]) ** GAMMA_HAT[r] / GAMMA_HAT[r]
        - lab_sec_Q[r, t]);
var utility_SumShr_Q 'utility: SumShr for consumption, quadratic for labour'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_SumShr[r, t] - lab_sec_Q[r, t]);
var utility_SumPow_Q 'utility: SumPow for consumption and quadratic for labour'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_SumPow[r, t] - lab_sec_Q[r, t]);
var utility_CD_F 'Cobb--Douglas and Frisch instantaneous utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CD[r, t] - lab_sec_caveF[r, t]);
var utility_CES_bangF 'Const. Elast. Subst. and conc. Frisch instant. utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (con_sec_CES[r, t] - lab_sec_bangF[r, t]);
var utility_CES_caveF 'Const. Elast. Subst. and conc. Frisch instant. utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (
        A_CON * con_sec_CES[r, t]
        + A_LAB_EXT * lab_ext_sec[r, t]
        + A_LAB[r, t] * lab_sec_caveF[r, t]
      );
var utility_CD_caveF 'Cobb-Douglas and concave Frisch inst. utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (
        A_CON * con_sec_CD[r, t]
        + A_LAB_EXT * lab_ext_sec[r, t]
        + A_LAB[r, t] * lab_sec_caveF[r, t]
      );
var cutility_CD_caveF 'Cobb-Douglas-Leontief and concave Frisch inst. utility'
  {t in LookForward}
  = sum{r in Regions}
      REG_WGHT[r] * (
        A_CON * ccon_sec_CD[r, t]
        + A_LAB_EXT * lab_ext_sec[r, t]
      #  + A_LAB[r, t] * lab_sec_caveF[r, t]
      );
#-----------variety of tail or terminal value functions
var tail_val_CD_F 'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions}
    REG_WGHT[r] * (
      prod{i in Sectors}(
        (TAIL_SHR_CON * kap[r, i, LSup + LInf] ** (SHR_KAP_OUT[r, i] * SCALE_OUT))
          ** (SHR_CON[r, i] * SCALE_CON)
      )
      - sum{i in Sectors} A_LAB[r, LSup + LInf] * 1 ** RHO_LAB / RHO_LAB
    )) / (1 - BETA);
var tail_val_CD_Q 'SumShr continuation value from time LSup + LInf onwards'
  = (sum{r in Regions}
    REG_WGHT[r] * (
      prod{i in Sectors}(
        (TAIL_SHR_CON * kap[r, i, LSup + LInf] ** (SHR_KAP_OUT[r, i] * SCALE_OUT))
          ** (SHR_CON[r, i] * SCALE_CON)
      )
      - 1 ** 2
    )) / (1 - BETA);
var tail_val_CESutl_Q_CDout 'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * ((
    sum{i in Sectors} SHR_CON_CES[r, i] * (TAIL_SHR_CON 
      * A[i] * kap[r, i, LSup + LInf] ** (SHR_KAP_OUT[r, i] * SCALE_OUT)
    ) ** RHO_CON
    ) ** (RHO_CON_HAT * SCALE_CON) #(RHO_INV_HAT + GAMMA_HAT[r]) / GAMMA_HAT[r]
    - 1 ** 2
  )) / (1 - BETA);
var tail_val_CDutl_F_CESout
  'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * (
    prod{i in Sectors} (TAIL_SHR_CON * A[i] * (
        SHR_KAP_OUT_CES[r, i] * kap[r, i, LSup + LInf] ** RHO_OUT
        + SHR_MED_OUT_CES[r, i] * 1 ** RHO_OUT
        + SHR_LAB_OUT_CES[r, i] * 1 ** RHO_OUT
      ) ** (RHO_OUT_HAT * SCALE_OUT) 
    )  ** (SHR_CON[r, i] * SCALE_CON)
    - sum{i in Sectors} A_LAB[r, LSup + LInf] * 1 ** RHO_LAB / RHO_LAB
  )) / (1 - BETA);
var tail_val_CESutl_bangF_CESout
  'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * ( A_CON * (
    sum{i in Sectors} SHR_CON_CES[r, i] * (TAIL_SHR_CON * A[i] * (
      SHR_KAP_OUT_CES[r, i] * kap[r, i, LSup + LInf] ** RHO_OUT
      + SHR_MED_OUT_CES[r, i] * 1 ** RHO_OUT
      + SHR_LAB_OUT_CES[r, i] * 1 ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT)) ** RHO_CON
  ) ** (RHO_CON_HAT * SCALE_CON) 
    #- (sum{i in Sectors} 1) ** RHO_LAB / RHO_LAB
    - A_LAB[r, LSup + LInf] * (sum{i in Sectors} 1) ** RHO_LAB
  )) / (1 - BETA);
var tail_val_CESutl_caveF_CESout
  'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * ( A_CON * (
    sum{i in Sectors} SHR_CON_CES[r, i] * (TAIL_SHR_CON * A[i] * (
      SHR_KAP_OUT_CES[r, i] * kap[r, i, LSup + LInf] ** RHO_OUT
      + SHR_MED_OUT_CES[r, i] * 1 ** RHO_OUT
      + SHR_LAB_OUT_CES[r, i] * (1 * ALPHA * ALPHA_0 ** (LSup + LInf)) ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT)) ** RHO_CON
  ) ** (RHO_CON_HAT * SCALE_CON) 
    + A_LAB_EXT * 0
    + A_LAB[r, LSup + LInf] * (
      sum{i in Sectors} SHR_LAB_CES[r, i] * 33e-2 ** RHO_LAB
      ) ** (RHO_LAB_HAT * SCALE_LAB)
  )) / (1 - BETA);
var tail_val_CDutl_caveF_CESout
  'continuation value from time LSup + LInf onwards'
  = (sum{r in Regions} REG_WGHT[r] * ( A_CON * (
    prod{i in Sectors} (TAIL_SHR_CON * A[i] * (
      SHR_KAP_OUT_CES[r, i] * kap[r, i, LSup + LInf] ** RHO_OUT
      + SHR_MED_OUT_CES[r, i] * 1 ** RHO_OUT
      + SHR_LAB_OUT_CES[r, i] * (1 ) ** RHO_OUT
      + 1 * (33e-2 * ALPHA * ALPHA_0 ** (LSup + LInf)) ** RHO_OUT
    ) ** (RHO_OUT_HAT * SCALE_OUT)
    ) ** (SHR_CON[r, i] * SCALE_CON)
    )
    + A_LAB_EXT * 0
    + A_LAB[r, LSup + LInf] * (
      sum{i in Sectors} SHR_LAB_CES[r, i] * 33e-2 ** RHO_LAB
      ) ** (RHO_LAB_HAT * SCALE_LAB)
  )) / (1 - BETA);
/*=============================================================================
Current intermediate variables (substituted out during pre-solving)
=============================================================================*/
var con 'consumption'
  {r in Regions, i in Sectors, t in LookForward}
  = ccon[r, i, t];
var inv 'domestic investment flows'
  {r in Regions, i in Sectors, j in Sectors, t in LookForward}
  = dinv[r, i, j, t];
var med 'domestic intermedate flows'
  {r in Regions, i in Sectors, j in Sectors, t in LookForward}
  = cmed[r, i, j, t];
var inv_sec 'current intermediate variable for aggregated investment'
  {r in Regions, j in Sectors, t in LookForward}
    = inv_sec_CES[r, j, t] * A_INV;
# for nimrod:
    #= __INVSEC__[r, j, t] * A_INV;
#var med_sec 'current intermediate variable for aggregated investment'
#  {r in Regions, j in Sectors, t in LookForward}
#    = cmed_sec_CES[r, j, t] * A_INV;
var E_out 'current intermediate variable for output'
  {r in Regions, i in Sectors, t in LookForward}
    = E_cout_CES[r, i, t];
var dom 'domestic uses'
  {r in Regions, i in Sectors, t in LookForward}
  = shr_dom[r, i, t] * E_out[r, i, t];
var xpo 'exports'
  {r in Regions, i in Sectors, t in LookForward}
  = (1 - shr_dom[r, i, t]) * E_out[r, i, t];
var utility 'current intermediate variable for utility'
  {t in LookForward} = cutility_CD_caveF[t]; 
var adj_cost_kap 'current adjustment costs for kapital'
  {r in Regions, i in Sectors, t in LookForward}
    = adj_cost_kap_Q[r, i, t];
var tail_val 'current intermediate variable for tail value function'
  = tail_val_CDutl_caveF_CESout * A_VAL;
/*=============================================================================
The objectives and constraints
=============================================================================*/
maximize pres_disc_val 'present discounted value of utilities':
    sum{t in LookForward} BETA ** (t - LInf) * utility[t]
      + BETA ** (LSup - LInf) * tail_val;
subject to kap_transition 'equation for the accumulation of kapital'
  {r in Regions, j in Sectors, t in LookForward}:
    kap[r, j, t + 1] - (1 - DELTA[j]) * kap[r, j, t] - inv_sec[r, j, t] = 0;
subject to market_clearing 'market clearing for each sector and time'
  {i in Sectors, t in LookForward}:
    sum{r in Regions}(
      con[r, i, t]
      + sum{j in Sectors}(inv[r, i, j, t])
      + sum{j in Sectors}(med[r, i, j, t])
      + adj_cost_kap[r, i, t]
      - dom[r, i , t] 
      ) = 0;
var mprod_fac
  {r in Regions, i in Sectors, t in LookForward}
  = SCALE_OUT * A[i] * (A[i] / E_out[r, i, t]) ** (RHO_OUT / SCALE_OUT - 1);
var mpkk
  {r in Regions, i in Sectors, t in LookForward}
  = mprod_fac[r, i, t] * SHR_KAP_OUT[r, i] * kap[r, i, t] ** (RHO_OUT );
var mpmm
  {r in Regions, i in Sectors, t in LookForward}
  = mprod_fac[r, i, t] * SHR_MED_OUT[r, i]
    * cmed_sec_CES[r, i, t] ** (RHO_OUT );
var mpll
  {r in Regions, i in Sectors, t in LookForward}
  = mprod_fac[r, i, t] * SHR_LAB_OUT[r, i] * lab[r, i, t] ** (RHO_OUT );
subject to vertical_balance
  {r in Regions, i in Sectors, t in LookForward}:
  E_out[r, i, t] >= mpkk[r, i, t] + mpmm[r, i, t] + mpll[r, i, t];

#subject to jacobi_id 'Intertemporal constraints on investment'
#  {r in Regions, i in Sectors, j in Sectors, t in LookForward
#    ii in Sectors: 1 < ord(j) and i <> j}:
#      inv[r, i, j, t]
#        = (inv[r, i, ii, t] / SHR_INV_CES[r, i, ii])
#          / (inv[r, ii, ii, t] / SHR_INV_CES[r, ii, ii])
#          * (inv[r, ii, j, t] / SHR_INV_CES[r, ii, j])
#          * SHR_INV_CES[r, i, j];
/*=============================================================================
The_data section
=============================================================================*/
param datadir symbolic; #"the directory where parameter-input data is stored" symbolic;
param outputdir symbolic;
data;
#-----------1xANZdiv model
set Regions := GLD;
set Sectors := A B C D E F G H I J K L M N O P Q R S;
#set Sectors := Agrc Frst;
#-----------2x2 model
#set Regions := GLD RoQ;
#set Sectors := Agrc Frst;
#-----------2x3 model
#set Regions := GLD RoQ;
#set Sectors := A Mnfc O;
#-----------2x4 model
#set Regions := GLD RoQ;
#set Sectors := A Frst Mnfc O;
#-----------3x4 model
#set Regions := GLD RoQ RoA;
#set Sectors := Agrc Frst Mnfc Srvc;
#-----------3x5 model
#set Regions := GLD RoQ RoA;
#set Sectors := A Frst Mnfc Srvc O;
#-----------3x6 model
#set Regions := GLD RoQ RoA;
#set Sectors := A Frst Mnfc Srvc O Utlt;
##-----------4x4
#set Regions := GLD RoQ RoA RoW;
#set Sectors := A Frst Mnfc O;# Trns Utlt;
##-----------4x5
#set Regions := GLD RoQ RoA RoW;
#set Sectors := A Frst Mnfc O Utlt;
##-----------4x6
#set Regions := GLD RoQ RoA RoW;
#set Sectors := A Elct Frst Mnfc Srvc O Utlt;
##-----------4x7
#set Regions := GLD RoQ RoA RoW;
#set Sectors := A B C D E F O;
##-----------4x8
#set Regions := GLD RoQ RoA RoW;
#set Sectors := A B C D E F G O;
##-----------4x10
#set Regions := GLD RoQ RoA RoW;
#set Sectors := A B C D E F G H I O;
#-----------7x18
#set Regions := CnQ FNQ Mck NrQ GLD WBB RoA;
#set Sectors := A B C D E F G H I J K L M N O P Q R;
#-----------7x20
#set Regions := CnQ FNQ Mck NrQ GLD WBB RoA;
#set Sectors := A B C D E F G H I J K L M N O P Q R T U;
#-----------7x21
#set Regions := CnQ FNQ Mck NrQ GLD WBB RoA;
#set Sectors := A B C D E F G H I J K L M N O P Q R T U Al;
##-----------7x40
#set Regions := CnQ FNQ Mck NrQ GLD WBB RoA;
#set Sectors := A B C D E F G H I J K L M N O P Q R T U
#  A1 B1 C1 D1 E1 F1 G1 H1 I1 J1 K1 L1 M1 N1 O1 P1 Q1 R1 T1 U1;
##-----------7x60
#set Regions := CnQ FNQ Mck NrQ GLD WBB RoA;
#set Sectors := A B C D E F G H I J K L M N O P Q R T U
#  A1 B1 C1 D1 E1 F1 G1 H1 I1 J1 K1 L1 M1 N1 O1 P1 Q1 R1 T1 U1
#  A2 B2 C2 D2 E2 F2 G2 H2 I2 J2 K2 L2 M2 N2 O2 P2 Q2 R2 T2 U2;
#=============================================================================#
# read in the data
#=============================================================================#
load amplcsv.dll;
let datadir := "julia/output/";
#-----------------------------------------------------------------------------#
# three indices: regions sectors and sectors: flows between sectors
#-----------------------------------------------------------------------------#
table rif IN "amplcsv"
(datadir & "RAW_INV_FLW.csv"): [index0, index1, index2], RAW_INV_FLW;
read table rif;
table rmf IN "amplcsv"
(datadir & "RAW_MED_FLW.csv"): [index0, index1, index2], RAW_MED_FLW;
read table rmf;
table rdcm IN "amplcsv"
(datadir & "RAW_DOM_CMED.csv"): [index0, index1, index2], RAW_DOM_CMED;
read table rdcm;
table rycm IN "amplcsv"
(datadir & "RAW_YPO_CMED.csv"): [index0, index1, index2], RAW_YPO_CMED;
read table rycm;
table smr IN "amplcsv"
(datadir & "SHR_MED_ROW.csv"): [index0, index1, index2], SHR_MED_ROW;
read table smr;
table smc IN "amplcsv"
(datadir & "SHR_MED_COL.csv"): [index0, index1, index2], SHR_MED_COL;
read table smc;
#-----------------------------------------------------------------------------#
# two indices: regions and sectors
#-----------------------------------------------------------------------------#
table rcf IN "amplcsv"
(datadir & "RAW_CON_FLW.csv"): [index0, index1], RAW_CON_FLW;
read table rcf;
table rdcc IN "amplcsv"
(datadir & "RAW_DOM_CCON.csv"): [index0, index1], RAW_DOM_CCON;
read table rdcc;
table rycc IN "amplcsv"
(datadir & "RAW_YPO_CCON.csv"): [index0, index1], RAW_YPO_CCON;
read table rycc;
table rko IN "amplcsv"
(datadir & "RAW_KAP_OUT.csv"): [index0, index1], RAW_KAP_OUT;
read table rko;
table rlo IN "amplcsv"
(datadir & "RAW_LAB_OUT.csv"): [index0, index1], RAW_LAB_OUT;
read table rlo;
table rmo IN "amplcsv"
(datadir & "RAW_MED_OUT.csv"): [index0, index1], RAW_MED_OUT;
read table rmo;
table rl IN "amplcsv"
(datadir & "REG_LAB.csv"): [index0, index1], REG_LAB;
read table rl;
#table fl IN "amplcsv"
#(datadir & "REF_LAB.csv"): [index0, index1], REF_LAB;
#read table fl;
table rxjo IN "amplcsv"
(datadir & "RAW_XPO_JOUT.csv"): [index0, index1], RAW_XPO_JOUT;
read table rxjo;
table rdjo IN "amplcsv"
(datadir & "RAW_DOM_JOUT.csv"): [index0, index1], RAW_DOM_JOUT;
read table rdjo;
table rors IN "amplcsv"
(datadir & "RAW_OUT_REG_SEC.csv"): [index0, index1], RAW_OUT_REG_SEC;
read table rors;
display RAW_MED_OUT;
#------------------------------------------------------------------------------
#-----------set the horizon and length of paths
#------------------------------------------------------------------------------
let LSup := 15;
let PSup := 105;
#------------------------------------------------------------------------------
#-----------opportunity to tune the calibration factors (still part of data)
#------------------------------------------------------------------------------
let BETA := 950e-3;
#param aA := 094e-2;
#param c20A := 090e-2;
#param c10A := 100e-2;
#let A["P"] := 40e-2;
#let A["H"] := 20e-2;
#let A["N"] := 20e-2;
#let A["R"] := 10e-2;
#let A["S"] := 20e-2;
for {i in Sectors}{
  let DELTA[i] := 05e-2;
  let PHI_ADJ[i] := 400e-2;
#  let A[i] := 015e-2;
  let A[i] := 2e+1 +
    (SHR_MED_ROW['GLD', i, i]
      / SHR_MED_COL['GLD', i, i]) ** ((1 - RHO_OUT) / RHO_OUT)
      * SCALE_OUT ** (1 / RHO_OUT)
    + 8 * RAW_OUT_REG_SEC['GLD', i]
      / sum{j in Sectors} RAW_OUT_REG_SEC['GLD', j];
  let EPS_JOUT['GLD', i] := 10e-2;
  for {j in Sectors}{
  let EPS_CMED['GLD', i, j] := 10e-2;
  };
};
#-----------------------------------------------------------------------------#
# initial kap and growth factor if starting from a future point in time
#-----------------------------------------------------------------------------#
let ALPHA := 1;#271828182846e-11;
let ALPHA_0 := 1;#271828182846e-11;
let ALPHA := 101828182846e-11 ** 105;
let ALPHA_0 := 101828182846e-11;
load amplxl.dll;
let datadir := "ampl/data/";
table kapxl IN "amplxl"
  (datadir & "KAP-105-fin" & ".xlsx") "Sheet1":
#  (datadir & "KAP-71-A15-PHI0-ALPHA1058" & ".xlsx") "Sheet1":
  [Regions, Sectors, PathTimes], KAP;
read table kapxl; #table inkap IN "amplcsv"
# to continue the previous loop
#for {r in Regions, i in Sectors} let KAP[r, i, 0] := KAP[r, i, 36];
#(datadir & "KAP.csv"): [Regions, Sectors, PathTimes], KAP;
#read table inkap;

#-----------------------------------------------------------------------------#
# initial kap and growth factor if starting from scratch (comment out otherwise)
#-----------------------------------------------------------------------------#
let ALPHA := 101828182846e-11;
let ALPHA_0 := 101828182846e-11;
for {r in Regions, i in Sectors}{
  let KAP[r, i, PInf]
    := 1;
#    := RAW_KAP_OUT[r, i] / sum{j in Sectors} RAW_KAP_OUT[r, j]
#      + RAW_OUT_REG_SEC[r, i] / sum{j in Sectors} RAW_OUT_REG_SEC[r, j]
#      + 50e-2
#      ;
};
display KAP;
#-----------------------------------------------------------------------------#
# regionalisation
#-----------------------------------------------------------------------------#
# no shock
#-----------------------------------------------------------------------------#
let A['C'] := 130e-2 * A['C'];
#let A['M'] := 130e-2;
let A['L'] := 070e-2 * A['L'];
let A['K'] := 070e-2 * A['K'];
let RAW_XPO_JOUT['GLD', 'C'] := RAW_DOM_JOUT['GLD', 'C'] * 3; 
let RAW_MED_FLW['GLD', 'C', 'C'] := RAW_MED_FLW['GLD', 'C', 'C'] * 130e-2;
let RAW_MED_FLW['GLD', 'D', 'C'] := RAW_MED_FLW['GLD', 'D', 'C'] * 700e-2;
let RAW_YPO_CMED['GLD', 'B', 'C'] := RAW_YPO_CMED['GLD', 'B', 'C'] * 130e-2;
#-----------------------------------------------------------------------------#
# the shock
#-----------------------------------------------------------------------------#
#let A['C'] := A['C'] * 080e-2;
#let RAW_XPO_JOUT['GLD', 'C'] := RAW_XPO_JOUT['GLD', 'C'] * 3 / 4; 
#let RAW_MED_FLW['GLD', 'C', 'C'] := RAW_MED_FLW['GLD', 'C', 'C'] * 040e-2;
#let RAW_MED_FLW['GLD', 'D', 'C'] := RAW_MED_FLW['GLD', 'D', 'C'] * 15e-2;
#let KAP['GLD', 'C', 0] := KAP['GLD', 'C', 0] * 3 / 4;
display A;
#let A['B'] := 50e-2;
#let KAP['GLD', 'B', 0] := 200e-2;
#let KAP['GLD', 'R', 0] := 20e-2;
#let KAP['GLD', 'S', 0] := 40e-2;
#let KAP['GLD', 'N', 0] := 40e-2;
#let KAP['GLD', 'H', 0] := 50e-2;
#let KAP['GLD', 'P', 0] := 50e-2;
let A_CON := 00100e-2; #increase this to increase labour and consumption
let A_INV := 0090e-2;
let A_MED := 0010e-2;
let A_VAL := 0001e-2;
let A_LAB_EXT := -0545e-2;
#let A_CMED := 1;
let TAIL_SHR_CON := 050e-2;

let EPS_INV := 0300e-3;
let EPS_MED := 0400e-3;
let EPS_CON := 0999e-3;
let EPS_OUT := 0800e-3;
let EPS_LAB := 050e-2;

let SCALE_CON := 300e-3;
let SCALE_INV := 999e-3;
let SCALE_MED := 999e-3;
let SCALE_OUT := 999e-3;
let SCALE_LAB := 600e-2;
let SCALE_CMED := 990e-3;
let SCALE_CINV := 990e-3;
for {r in Regions, i in Sectors, j in Sectors, t in LookForward}{
#  fix lab[r, j, t] := 33e-2;
  fix lab_ext[r, j, t] := 100e-2;
  let NAIRE[r, j, t] := 95e-2;
  let EXP_LAB_EXT[r, j, t] := 2;
  if SHR_INV_CES[r, i, j] < 1e-12 then
    fix dinv[r, i, j, t] := 0;
  if SHR_MED_CES[r, i, j] < 1e-12 then
    fix dmed[r, i, j, t] := 0;
};
update data;
param CH_GROWTH_OUT {Regions, Sectors, PathTimes} default 0;
param nwshr 'Alternative share parameter'
  {r in Regions, i in Sectors}
  = A[i] ** (SCALE_OUT - RHO_OUT - 1)
          * RAW_KAP_OUT[r, i] ** (1 - RHO_OUT)
          / (RAW_LAB_OUT[r, i]
             + RAW_KAP_OUT[r, i]
             + RAW_MED_OUT[r, i]) ** (SCALE_OUT - RHO_OUT);
param nwshr2 'Alternative share parameter for square parameters'
  {r in Regions, i in Sectors, j in Sectors}
  = A[i] ** (SCALE_INV - RHO_INV - 1)
          * RAW_INV_FLW[r, i, j] ** (1 - RHO_INV)
          / (sum{ii in Sectors} RAW_INV_FLW[r, ii, j]
            ) ** (SCALE_INV - RHO_INV);
/*=============================================================================
Solving the model
=============================================================================*/
param InstanceName symbolic;
option solver knitro;
let outputdir := "ampl/output/";
#option solver ipopt;
option show_stats 1;
#=============================================================================#
# setup for output of results
#=============================================================================#
#==============================================================================
#-----------parameters for storing (observable) path values
#==============================================================================
param CON 'observed consumption' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup);
param INV_SEC 'observed investment' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup);
param INV_SUM 'observed total investment' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup); 
param MED_SUM 'observed total intermediate flows' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup); 
param LAB 'observed labour' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup);
param LAB_EXT 'observed laborforce participation' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup);
param E_OUT 'observed Exp. output' {Regions, Sectors, PathTimes}
  default 1e+0; # in (OInf, OSup); 
param ADJ_COST_KAP 'observed adjustment costs of kapital'
  {Regions, Sectors, PathTimes} default 0; # in [0, OSup);
param MKT_CLR 'observed output' {Sectors, PathTimes}
  default 0; # in (-1e-4, 1e-4); 
param DUAL_KAP 'lagrange multiplier for kapital accumulation'
  {Regions, Sectors, PathTimes} default 1e+0; # in (-OSup, OSup);
param DUAL_MKT_CLR 'lagrange multiplier for market clearing constraint'
  {Sectors, PathTimes} default 1e+0;# in [0, OSup);
param GROWTH_KAP 'observed growth rate for kapital'
  {Regions, Sectors, PathTimes} default 5e-2; # in (-1, 1);
param GROWTH_OUT 'observed growth rate for output'
  {Regions, Sectors, PathTimes} default 5e-2; # in (-1, 1);
param EULER_INTEGRAND 'Euler error integrand'
  {Regions, Sectors, PathTimesClosure} default 1; # in (-OSup, OSup);
param EULER_RATIO 'Expected Euler ratio'
  {Regions, Sectors, PathTimes} default 1; # in (-1e+2, 1e+2);
param DOM 'actual path values for domestic output'
  {Regions, Sectors, PathTimes}
  default 100e-2;
param XPO 'actual path values for exports'
  {Regions, Sectors, PathTimes}
  default 100e-2;
param YMED_CSUM 'actual path values for sum over a row of intermediate imports'
  {Regions, Sectors, PathTimes}
  default 100e-2;
param CMED_SEC 'actual path values for intermediate input aggregator'
  {Regions, Sectors, PathTimes}
  default 100e-2;
param MPROD_FAC 'marginal product factor (just output in the Cobb--Doug case)'
  {r in Regions, i in Sectors, s in PathTimes}
  = SCALE_OUT * A[i] * (A[i] / E_OUT[r, i, s]) ** (RHO_OUT / SCALE_OUT - 1);
param MPKK 'marginal product of kapital'
  {r in Regions, i in Sectors, s in PathTimes}
  = MPROD_FAC[r, i, s] * SHR_KAP_OUT[r, i] * KAP[r, i, s] ** (RHO_OUT );
param MPMM 'marginal product of intermediate (the aggregator)'
  {r in Regions, i in Sectors, s in PathTimes}
  = MPROD_FAC[r, i, s] * SHR_MED_OUT[r, i]
    * CMED_SEC[r, i, s] ** (RHO_OUT );
param MPLL 'marginal product of labour'
  {r in Regions, i in Sectors, s in PathTimes}
  = E_OUT[r, i, s] - MPKK[r, i, s] - MPMM[r, i, s];
param GVA 'gross value added'
  {r in Regions, i in Sectors, s in PathTimes}
  = MPKK[r, i, s] + MPLL[r, i, s];
param VERT_BAL 'vertical balance (should be zero)'
  {r in Regions, i in Sectors, s in PathTimes}
  = E_OUT[r, i, s] - GVA[r, i, s] - MPMM[r, i, s];
param AGG_OUT 'aggregate output per period'
  {r in Regions, s in PathTimes}
  = sum{i in Sectors} E_OUT[r, i, s];
param AGG_KAP 'aggregate kapital per period'
  {r in Regions, s in PathTimes}
  = sum{i in Sectors} KAP[r, i, s];
param AGG_CON 'aggregate consumption per period'
  {r in Regions, s in PathTimes}
  = sum{i in Sectors} CON[r, i, s];
param AGG_XPO 'aggregate exports per period'
  {r in Regions, s in PathTimes}
  = sum{i in Sectors} XPO[r, i, s];
param AGG_YMED_CSUM 'aggregate imports indirect alloc. (column sum) per period'
  {r in Regions, s in PathTimes}
  = sum{i in Sectors} YMED_CSUM[r, i, s];
param experiment symbolic;
param shock symbolic;
param reg symbolic;
#option solver conopt;
let shock := "historic";
let reg := "aus";
let outputdir := ("ampl/output/" & reg & shock & "/");
# declare output tables
table res OUT "amplxl" (outputdir & "burn-in-8-C130.xlsx") "Results":
  [Regions, Sectors, PathTimes],
  KAP, E_OUT, CON, XPO, YMED_CSUM, GROWTH_KAP, GROWTH_OUT,
  EULER_INTEGRAND, EULER_RATIO;
table aggres OUT "amplxl" (outputdir & "agg-burn-in-8-C130.xlsx") "Results":
  [Regions, PathTimes],
  AGG_KAP, AGG_CON, AGG_XPO, AGG_YMED_CSUM;
#=============================================================================#
# the above model may be solved in isolation, but to solve a path issue the
# command: "include gladpath.run" after instantiating this model
# or, on the command line: "ampl <this file> path.run"
#=============================================================================#
