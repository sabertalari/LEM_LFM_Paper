$CALL GDXXRW dataexcelWind_new3_twostage_newDR_ECM1n.xlsx Index=ParamsAndSets!A1
$GDXIN dataexcelWind_new3_twostage_newDR_ECM1n.gdx

Sets
t           Time interval
w
S
N_PV
PEVscen     Scenarios of PEVS
*For 20 PEV scenario add /1*20/ to the excel in set section and change the probability!A3 to probability!D3
Windscen    Wind scenarios
PVscen      PhotoVoltaic scenarios


N                                  Indices for the bus system                   /1*15/
N_DG(N)                            Conventional generation node position
N_REN(N)                           Renewable generation node position
N_REN2(N)                           Renewable generation node position
N_TN(N)
SL(N)                              Slack bar indice
R                                  Number of power section P y Q
SWITCH(N,N)                        Indice for switches
CONEX(N,N)                         Indice for line connection
CONEX_SIN(N,N)                     Indice for lines with switch
PARL                               Line parameters
PDQ                                Set of active and reactive demand
Prob(s,w)
CONDR(N,N)
N_DR(N)
N_ES(N)
;
$load t,PVscen,Windscen,w,N_PV,S, R SL conex conex_sin N_DG N_REN N_REN2 PARL PDQ Prob SWITCH N_TN CONDR N_DR N_ES
alias(N,NP,i,j)  ;
alias(t,tAR,tDEP);
*alias(t,k);

Scalars
Wind               To be or not to be!
PV                 To be or not to be!
VMIN                               Minimum value of the grid voltage            [kV]
VMAX                               Maximum value of the grid voltage            [kV]
IMAX                               Maximum value of the grid current            [A]
IMIN                               Minimum value of the grid current            [A]
VNOM                               Nominal voltage of the grid                  [kV]
VMIN2                              Minimum value of the grid voltage square     [kV]
VMAX2                              Maximum value of the grid voltage square     [kV]
IMIN2                              Valor de corriente m?nimo en la red          [A]
VNOM2                              Voltage nominal al-cuadrado                  [kV]
W_MAX2                             Maximum limit for switches voltage drop      [V]
*P_DG_MAX                           Maximum conventional active power            [kW]
Q_DG_MAX                           Maximum conventional reactive power          [kVAr]
P_RN_MAX                           Maximum wind active power                    [kW]
Q_RN_MAX                           Maximum wind reactive power                  [kW]
P_RN_MAX2                           Maximum PV active power                    [kW]
Q_RN_MAX2                           Maximum PV reactive power                  [kW]
Sbase                                                                               /2300/
SOC0                               initial state of the charge              [pu]    /0.04347/
SOCmax                             maximum state of the charge              [kW]    /0.13/
eta_c                              efficiency of charging                    [%]    /0.95/
eta_d                              efficiency of discharging                 [%]    /0.9/;

$Load Wind,PV

Parameters
Demand_P(t)           Electric demand without PEVs consideration
Windexp(t,w)
PVexp(t,w)

LineData(N,NP,PARL)                Line branch parameters
demand(N,t,PDQ)                    Demand at each node and hour
Pdem(N,PDQ,t)                      Active demand at node and for time indice    [kW]
Qdem(N,PDQ,t)                      Reactive demand at node and for time indice  [kVAr]


NumEolica                          Number of wind parks                         /1/
P_RN(N,t,w)                        Wind active power                       [kW]
P_RN2(N,t,w)                       PV active power                       [kW]
P_DG_MAX(N,w)
P_DG_MIN(N,w)
PV_S(t,S)
Wind_S(t,S)
P_RNS(N,t,w,s)
P_RNS2(N,t,w,s)

Price_Energy(t)       Tariff of DSO OR energy market(Day-ahead) Euro per kWh
Price_Reserve(t)      Price of reserve
P_del_Res(t)          Probability of being called by DSO OR ISO to generate its offered reserve
Price_Reg_Up(t)       Regulation market (Notice: typical data have been used!!)
P_del_Reg_Up(t)       Probability of being called by DSO OR ISO to generate its offered regulation (Notice: typical data have been used!!)
Price_Reg_Down(t)     Regulation market (Notice: typical data have been used!!)
P_del_Reg_Down(t)     Probability of being called by DSO OR ISO to generate its offered regulation (Notice: typical data have been used!!)
Cap_payment_price(t)  Capacity reserve supply OR ancillary services market
Cap_Reg_Up(t)         Capacity payment regulation up
Cap_Reg_Down(t)       Capacity payment regulation down

V0_2(N,t,w)                        Initialization of the voltage at node n      [kV]
DELTA_S(N,NP,t,w)                  Variation of aparent power between nodes     [kW]
testsecondstage(N,t,w,s)
*DR-------------
Lambda_DR(t,N)
;

$Load Demand_P,Windexp, PVexp,  LineData,  demand,   P_DG_MAX, P_DG_MIN, PV_S, Wind_S  Price_Energy,Price_Reserve,P_del_Res,Price_Reg_Up,Price_Reg_Down,P_del_Reg_Down,P_del_Reg_Up,Cap_payment_price,Cap_Reg_Down,Cap_Reg_Up  ,Lambda_DR


P_RN(N,t,w)$ N_REN(N) = Windexp(t,w)/2300 ;
P_RN2(N,t,w)$ N_REN2(N)= PVexp(t,w)/2300;
P_RNS(N,t,w,s)$ N_REN(N)= Wind_S(t,S)/2300;
P_RNS2(N,t,w,s)$ N_REN2(N)= PV_S(t,S)/2300;
Pdem(N,PDQ,t) = demand(N,t,'Pt') ;
Qdem(N,PDQ,t) = demand(N,t,'Qt') ;
display Pdem,Qdem,P_RN, P_RN2,P_RNs, P_RNs2 ;

***-----------------------------------------------------------------------------------------------------------------**
variables

**-------------------------------------------------------------------------------------------------------------------**
fobj                               Funcion objetivo
WV(N,NP,t,w)                       Auxiliar variable for voltage drop           [V]
WV_S(N,NP,t,w,s)
P_DG(N,t,w)                        Conventional generation active power         [kW]
Q_DG(N,t,w)                        Conventional generation reactive power       [kVAr]
Q_RN(N,t,w)                        Wind generation reactive power               [kVAr]
Q_RN2(N,t,w)
Q_RNs(N,t,w,s)                        Wind generation reactive power               [kVAr]
Q_RN2s(N,t,w,s)
Pss(N,t,w)                         Substation active power                      [kW]
Qss(N,t,w)                         Substation reactive power                    [kVAr]
PSS_S(N,t,w,s)
Qss_S(N,t,w,s)
ProfitTest


DR_final(t)
load_final(t)
Each_node_load(t,n)
TOTALDR_COST_T(t)
TOTALDR_COST_N(N)

SOC(t,N)
Pd(t,N)
Pc(t,N)

SOC_s(t,N,s)
Pd_s(t,N,s)
Pc_s(t,N,s)
;

SOC.up(t,N)=SOCmax;
SOC.lo(t,N)=0.2*SOCmax;
SOC.fx('24',N)=SOC0;
Pd.up(t,N)=0.2*SOCmax;
Pd.lo(t,N)=0;
Pc.up(t,N)=0.2*SOCmax;
Pc.lo(t,N)=0;

SOC_s.up(t,N,s)=SOCmax;
SOC_s.lo(t,N,s)=0.2*SOCmax;
SOC_s.fx('24',N,s)=SOC0;
Pd_s.up(t,N,s)=0.2*SOCmax;
Pd_s.lo(t,N,s)=0;
Pc_s.up(t,N,s)=0.2*SOCmax;
Pc_s.lo(t,N,s)=0;


positive variable
p_En_G2PL(t,Windscen,PVscen)
,SOC_PL(t,Windscen,PVscen)
,p_En_PL2G(t,Windscen,PVscen)
,p_Res_PL2G(t,Windscen,PVscen)
,p_Reg_PL2G(t,Windscen,PVscen)
,p_Reg_G2PL(t,Windscen,PVscen)
,q_En_PL2G(t,windscen,PVscen)
,q_En_G2PL(t,windscen,PVscen)
,Res_up(N,t,w)
,Res_Dn(N,t,w)

,Res_Dns(N,t,w,s)
,Res_up(N,t,w)
,Res_ups(N,t,w,s)
,Q_Res_up(N,t,w)
,Q_Res_Dn(N,t,w)
,TEst(N,t,s)
,Q_Res_ups(N,t,w,s)
,Q_Res_DNs(N,t,w,s)
 ,DR(t)
,DRs(t,s)
,P_DR1(N,t)
,P_DR(N,NP,t)
,P_DRs(N,t,s)



**----------------------------------------------------------------------------------------**


,V2(N,t,w)                          Voltage square at each bus                   [kV]
,I2(N,NP,t,w)                       Current square at each branch                [A]
,P_POS(N,NP,t,w)                    Positive active power flow                   [kW]
,P_NEG(N,NP,t,w)                    Negative active power flow                   [kW]
,Q_POS(N,NP,t,w)                    Positive reactive power flow                 [kW]
,Q_NEG(N,NP,t,w)                    Negative reactive power flow                 [kW]
,DELTA_P(N,NP,R,t,w)                Active power discretization block            [kW]
,DELTA_Q(N,NP,R,t,w)                Rective power discretization block           [kW]
,P_POS_S(NP,N,t,w,s)
,P_NEG_S(NP,N,t,w,s)
,I2_S(N,NP,t,w,s)
,V2_S(NP,t,w,s)
,WV_S(N,NP,t,w,s)
,Q_POS_S(N,NP,t,w,s)
,Q_NEG_S(N,NP,t,w,s)
,DELTA_P_S(N,NP,R,t,w,s)
,DELTA_Q_S(N,NP,R,t,w,s)
,Y_Dual1(t)
,Y_Dual2(t,N)
,Y_Dual3(t,N,NP)


;

binary variable
u1_PL(t,Windscen,PVscen),
u2_PL(t,Windscen,PVscen),
u_DG(N,t,w)
**------------------------------------------------------------------------------------------**
,var_switch(N,NP,t)                 Binary variable for switches (1-ON 0-OFF)
;
equations


**------------------------------------------------------------------------------------
Max_Profit    ,
E_BAL_POT_P(N,t,PDQ,w,PVscen,Windscen)             Active power balance                         [kW],
constES(t, N) ,
constES_s(t, N,s)     ,
E_BAL_POT_Q(N,t,PDQ,w,PVscen,Windscen)             Reactive power balance                       [kW],

E_BAL_POT_P_S(N,t,PDQ,w,PVscen,Windscen,s)    ,
E_BAL_POT_Q_S,
E_BAL_V_S(N,NP,t,w,s),
E_POT_S_SS(N,NP,t,w,s),
E_POT_P_S(N,NP,t,w,s),
E_POT_Q_S,
E_COTA_DELTA_P_S(N,NP,R,t,w,s),
E_COTA_DELTA_Q_S(N,NP,R,t,w,s),
E_INT_MAX_S(N,NP,t,w,s),
E_W_MAX_S(N,NP,t,w,s),
E_W_MIN_S(N,NP,t,w,s),
E_L_POT1_S(N,NP,t,w,s),
E_L_POT2_S(N,NP,t,w,s),
E_RENCOS_S(N,t,w,s),
E_RENCOS1_S,
E_RENCOS2_S,
E_RENCOS12_S,
*E_SUB_S,
*E_SUB1_S,
*reserve_limit_dns_L(N,t,w,s),

E_BAL_V(N,Np,t,w)                  Node voltage balance                         [V],
E_POT_S(N,NP,t,w)                  Liniearization balance                       [kW],
E_POT_P(N,NP,t,w)                  Expression for the linearization block of active power    [kW],
E_POT_Q(N,NP,t,w)                  Expression for the linearization block of reactive power  [kVar]
E_COTA_DELTA_P(N,NP,R,t,w)
E_COTA_DELTA_Q(N,NP,R,t,w)
E_INT_MAX(N,NP,t,w)
E_SWITCH(N,NP,t)
E_RADIALIDAD                       Radiality constraints
E_W_MAX(N,NP,t,w)                  Upper limits of voltage drop in switches.                     [kV]
E_W_MIN(N,NP,t,w)                  Lower limits of voltage drop in switches.                     [kV]
E_L_POT1(N,NP,t,w)                 Active power limits for the whole network                     [kW]
E_L_POT2(N,NP,t,w)                 Reactive power limits for the whole network                   [kVAr]
E_RENCOS(N,t,w)                    Active & reactive power relation for positive renewable generation (wind)
E_RENCOS1(N,t,w)                   Active & reactive power relation for negative renewable generation (wind)
E_RENCOS2(N,t,w)                    PV
E_RENCOS12(N,t,w)                     PV
E_MARCOS(N,t,w)                    Active & reactive power relation for positive conventional generation
E_MARCOS1(N,t,w)                   Active & reactive power relation for negative conventional generation
E_SUB(N,t,w)                       Active & reactive power relation for positive substation
E_SUB1(N,t,w)                      Active & reactive power relation for negative substation
P_DG_limit(N,t,w)                  Conventional active power generation limit                    [kW]
reserve_limit_ups
reserve_limit_dns
reserve_limit_dns_L
reserve_limit_ups_L

Qreserve_limit_ups
Qreserve_limit_dns
Qreserve_limit_dns_L
Qreserve_limit_ups_L

DR_limit2
active_load
DR_limit3
DR_limit4
ESS_con1
ESS_con2
ESS_con3

    ;

Max_Profit..                                           ProfitTest=e=sum((w,t,N), Price_Energy(t)*PSS(N,t,w)$(ord(N)=1)+Price_Reserve(t)*P_DG(N,t,w)$N_DG(N)+ Lambda_DR(t,N)*P_DR1(N,t)$N_DR(N)+ Price_Reg_up(t)*(Res_up(N,t,w)+Res_Dn(N,t,w))$(ord(N)=1)+sum(s, prob(s,w)*( Lambda_DR(t,N)*P_DRs(N,t,s)$N_DR(N)+Price_Reg_down(t)*(Res_ups(N,t,w,s)+Res_Dns(N,t,w,s))$(ord(N)=1))));
E_BAL_POT_P(N,t,'pt',w,PVscen,Windscen) ..             NumEolica *(P_RN(N,t,w)$ N_REN(N))+ (P_DR1(N,t)$N_DR(N))+NumEolica *(P_RN2(N,t,w)$ N_REN2(N)) +Pd(t, N)$N_ES(N)-Pc(t,N)$N_ES(N)+  PSS(N,t,w) $(ord(N)=1)+SUM[NP $ CONEX(NP,N), (P_POS(NP,N,t,w) - P_NEG(NP,N,t,w))] - SUM[NP $ CONEX(N,NP), (P_POS(N,NP,t,w) - P_NEG(N,NP,t,w)) + LineData(N,NP,'RL')* I2(N,NP,t,w)] + (P_DG(N,t,w) $ N_DG(N))  =E= 1.1*Pdem(N,'pt',t) ;
*ES constraint
constES(t, N)$N_ES(N)..                                SOC(t,N)=e= SOC0$(ord(t)=1)+SOC(t-1, N)$(ord(t)>1)+Pc(t,N)*eta_c-Pd(t,N)/eta_d;
constES_s(t, N,s)$N_ES(N)..                                SOC_s(t,N,s)=e= SOC0$(ord(t)=1)+SOC_s(t-1, N,s)$(ord(t)>1)+Pc_s(t,N,s)*eta_c-Pd_s(t,N,s)/eta_d;
*Reactive Power Balance [KW]
E_BAL_POT_Q(N,t,'qt',w,PVscen,Windscen) ..             NumEolica *Q_RN(N,t,w)$ N_REN(N)+NumEolica *Q_RN2(N,t,w)$ N_REN2(N) +  QSS(N,t,w) $(ord(N)=1) + SUM[NP $ CONEX(NP,N), (Q_POS(NP,N,t,w) - Q_NEG(NP,N,t,w))] - SUM[NP $ CONEX(N,NP), (Q_POS(N,NP,t,w) - Q_NEG(N,NP,t,w)) + LineData(N,NP,'XL')* I2(N,NP,t,w)] + (Q_DG(N,t,w) $ N_DG(N))  =E= 1.1*Qdem(N,'qt',t) ;
*second stage----------------------------------------
E_BAL_POT_P_S(N,t,'pt',w,PVscen,Windscen,s) ..            Res_ups(N,t,w,s)$(ord(N)=1)-Res_dns(N,t,w,s)$(ord(N)=1)+Pd_s(t,N,S)$N_ES(N)-Pc_s(t,N,S)$N_ES(N) +NumEolica *(P_RNs(N,t,w,s)$ N_REN(N)-P_RN(N,t,w)$ N_REN(N))+P_DRs(N,t,s)$N_DR(N)+ NumEolica *(P_RNs2(N,t,w,s)$ N_REN2(N)-P_RN2(N,t,w)$ N_REN2(N))+(SUM[NP $ CONEX(NP,N), (P_POS_S(NP,N,t,w,s) - P_NEG_S(NP,N,t,w,s))]-SUM[NP $ CONEX(NP,N), (P_POS(NP,N,t,w) - P_NEG(NP,N,t,w))])-(SUM[NP $ CONEX(N,NP), (P_POS_S(N,NP,t,w,s) - P_NEG_S(N,NP,t,w,s)) + LineData(N,NP,'RL')* I2_S(N,NP,t,w,s)] - SUM[NP $ CONEX(N,NP), (P_POS(N,NP,t,w) - P_NEG(N,NP,t,w)) + LineData(N,NP,'RL')* I2_S(N,NP,t,w,s)])=E= 0 ;
*Reactive Power Balance-second-stage [KW]
E_BAL_POT_Q_S(N,t,'qt',w,PVscen,Windscen,s) ..          Q_Res_ups(N,t,w,s)$(ord(N)=1)-Q_Res_dns(N,t,w,s)$(ord(N)=1)+ NumEolica *(Q_RNs(N,t,w,s)$ N_REN(N)-Q_RN(N,t,w)$ N_REN(N))+ NumEolica *(Q_RN2s(N,t,w,s)$ N_REN2(N)-Q_RN2(N,t,w)$ N_REN2(N))+(SUM[NP $ CONEX(NP,N), (Q_POS_S(NP,N,t,w,s) - Q_NEG_S(NP,N,t,w,s))]-SUM[NP $ CONEX(NP,N), (Q_POS(NP,N,t,w) - Q_NEG(NP,N,t,w))])-(SUM[NP $ CONEX(N,NP), (Q_POS_S(N,NP,t,w,s) - Q_NEG_S(N,NP,t,w,s)) + LineData(N,NP,'XL')* I2_S(N,NP,t,w,s)] - SUM[NP $ CONEX(N,NP), (Q_POS(N,NP,t,w) - Q_NEG(N,NP,t,w)) + LineData(N,NP,'XL')* I2_S(N,NP,t,w,s)])=E= 0 ;
*Node voltage Balance [V]                                                                                                                  hg
E_BAL_V(N,NP,t,w) $ CONEX(N,NP)..                      V2(N,t,w) + WV(N,NP,t,w) - 2 * [LineData(N,NP,'RL') * (P_POS(N,NP,t,w)-P_NEG(N,NP,t,w))+ LineData(N,NP,'XL') * (Q_POS(N,NP,t,w)-Q_NEG(N,NP,t,w))]- [POWER(LineData(N,NP,'RL'),2) + POWER(LineData(N,NP,'XL'), 2)] * I2(N,NP,t,w) - V2(NP,t,w)=E=0;
*Node voltage Balance_Second stage [V]
E_BAL_V_S(N,NP,t,w,s) $ CONEX(N,NP)..                  V2_S(N,t,w,s) + WV_S(N,NP,t,w,s) - 2 * [LineData(N,NP,'RL') * (P_POS_S(N,NP,t,w,s)-P_NEG_S(N,NP,t,w,s))+ LineData(N,NP,'XL') * (Q_POS_S(N,NP,t,w,s)-Q_NEG_S(N,NP,t,w,s))]- [POWER(LineData(N,NP,'RL'),2) + POWER(LineData(N,NP,'XL'), 2)] * I2_S(N,NP,t,w,s) - V2_S(NP,t,w,s)=E=0;
*Linearizaton Balance [KW]
E_POT_S(N,NP,t,w) $ CONEX(N,NP)..                      V0_2(N,t,w) * I2(N,NP,t,w) =E= SUM[R, (2 * ord(R) - 1)* DELTA_S(N,NP,t,w) * DELTA_P(N,NP,R,t,w)] + SUM[R, (2 * ord(R) - 1)* DELTA_S(N,NP,t,w) * DELTA_Q(N,NP,R,t,w)];
*Linearizaton Balance -second stage[KW]
E_POT_S_SS(N,NP,t,w,s) $ CONEX(N,NP)..                      V0_2(N,t,w) * I2_S(N,NP,t,w,s) =E= SUM[R, (2 * ord(R) - 1)* DELTA_S(N,NP,t,w) * DELTA_P_S(N,NP,R,t,w,s)] + SUM[R, (2 * ord(R) - 1)* DELTA_S(N,NP,t,w) * DELTA_Q_S(N,NP,R,t,w,s)];
*Expression for the linearization block of active power [KW]
E_POT_P(N,NP,t,w) $ CONEX(N,NP)..                      P_POS(N,NP,t,w) + P_NEG(N,NP,t,w) =E= SUM[R, DELTA_P(N,NP,R,t,w)];
*Expression for the linearization block of active power-second stage [KW]
E_POT_P_S(N,NP,t,w,s) $ CONEX(N,NP)..                      P_POS_S(N,NP,t,w,s) + P_NEG_S(N,NP,t,w,s) =E= SUM[R, DELTA_P_S(N,NP,R,t,w,s)];
*Expression for the linearization block of the reactive power [KVar]
E_POT_Q(N,NP,t,w) $ CONEX(N,NP)..                      Q_POS(N,NP,t,w) + Q_NEG(N,NP,t,w) =E= SUM[R, DELTA_Q(N,NP,R,t,w)];
*Expression for the linearization block of the reactive power-second stage [KVar]
E_POT_Q_S(N,NP,t,w,s) $ CONEX(N,NP)..                      Q_POS_S(N,NP,t,w,s) + Q_NEG_S(N,NP,t,w,s) =E= SUM[R, DELTA_Q_S(N,NP,R,t,w,s)];

E_COTA_DELTA_P(N,NP,R,t,w) $ CONEX(N,NP)..             DELTA_P(N,NP,R,t,w) =L=  DELTA_S(N,NP,t,w);

E_COTA_DELTA_Q(N,NP,R,t,w) $ CONEX(N,NP)..             DELTA_Q(N,NP,R,t,w) =L=  DELTA_S(N,NP,t,w);
*second stage-------------------
E_COTA_DELTA_P_S(N,NP,R,t,w,s) $ CONEX(N,NP)..             DELTA_P_S(N,NP,R,t,w,s) =L=  DELTA_S(N,NP,t,w);

E_COTA_DELTA_Q_S(N,NP,R,t,w,s) $ CONEX(N,NP)..             DELTA_Q_S(N,NP,R,t,w,s) =L=  DELTA_S(N,NP,t,w);

E_INT_MAX(N,NP,t,w) $ CONEX(N,NP)..                    I2(N,NP,t,w) =L= LineData(N,NP,'IMAX') * LineData(N,NP,'IMAX') * var_switch(N,NP,t);
*second-stage
E_INT_MAX_S(N,NP,t,w,s) $ CONEX(N,NP)..                    I2_S(N,NP,t,w,s) =L= LineData(N,NP,'IMAX') * LineData(N,NP,'IMAX') * var_switch(N,NP,t);
E_SWITCH(N,NP,t) $ CONEX_SIN(N,NP) ..                  var_switch(N,NP,t) =E= 1;
*Upper limits of voltage drop in switches
E_W_MAX(N,NP,t,w) $ CONEX(N,NP)..                      WV(N,NP,t,w) =L= W_MAX2 * (1- var_switch(N,NP,t));
*Lower limits of voltage drop in switches
E_W_MIN(N,NP,t,w) $ CONEX(N,NP)..                      WV(N,NP,t,w) =G= -W_MAX2 * (1- var_switch(N,NP,t));
*Upper limits of voltage drop in switches-second-stage
E_W_MAX_S(N,NP,t,w,s) $ CONEX(N,NP)..                      WV_S(N,NP,t,w,s) =L= W_MAX2 * (1- var_switch(N,NP,t));
*Lower limits of voltage drop in switches-second-stage
E_W_MIN_S(N,NP,t,w,s) $ CONEX(N,NP)..                      WV_S(N,NP,t,w,s) =G= -W_MAX2 * (1- var_switch(N,NP,t));
*Active power limits for the whole network [KW]
E_L_POT1(N,NP,t,w) $ CONEX(N,NP) ..                    P_POS(N,NP,t,w) + P_NEG(N,NP,t,w) =L= VNOM * LineData(N,NP,'IMAX');
*Reactive power limits for the whole network [KW]
E_L_POT2(N,NP,t,w) $ CONEX(N,NP) ..                    Q_POS(N,NP,t,w) + Q_NEG(N,NP,t,w) =L= VNOM * LineData(N,NP,'IMAX');
*Active power limits for the whole network [KW]-second-stage
E_L_POT1_S(N,NP,t,w,s) $ CONEX(N,NP) ..                    P_POS_S(N,NP,t,w,s) + P_NEG_S(N,NP,t,w,s) =L= VNOM * LineData(N,NP,'IMAX');
*Reactive power limits for the whole network [KW]-second stage
E_L_POT2_S(N,NP,t,w,s) $ CONEX(N,NP) ..                    Q_POS_S(N,NP,t,w,s) + Q_NEG_S(N,NP,t,w,s) =L= VNOM * LineData(N,NP,'IMAX');
*_____________________________________________________________________________________________________________________________________________________________________________________________________
*Active and reactive power relation for positive renewable generation
E_RENCOS(N,t,w) $ N_REN(N) ..                          Q_RN(N,t,w) =L= P_RN(N,t,w) * TAN ( ARCCOS(0.95));
*Active and reactive power relation for negative renewable generation
E_RENCOS1(N,t,w) $ N_REN(N) ..                         Q_RN(N,t,w) =g= P_RN(N,t,w) * TAN ( ARCCOS(-0.95));
*second stage ------------------------------------------
*Active and reactive power relation for positive renewable generation
E_RENCOS_S(N,t,w,s) $ N_REN(N) ..                          Q_RNs(N,t,w,s) =L= P_RNs(N,t,w,s) * TAN ( ARCCOS(0.95));
*Active and reactive power relation for negative renewable generation
E_RENCOS1_S(N,t,w,s) $ N_REN(N) ..                         Q_RNs(N,t,w,s) =g= P_RNs(N,t,w,s) * TAN ( ARCCOS(-0.95));
*Active & reactive power relation for positive convetional generation
*PV------------------------------------
E_RENCOS2(N,t,w) $ N_REN2(N) ..                          Q_RN2(N,t,w) =L= P_RN2(N,t,w) * TAN ( ARCCOS(0.95));
*Active and reactive power relation for negative renewable generation
E_RENCOS12(N,t,w) $ N_REN2(N) ..                         Q_RN2(N,t,w) =g= P_RN2(N,t,w) * TAN ( ARCCOS(-0.95));
*tow-stage------------------------------------------------------------------------------------------------------
E_RENCOS2_S(N,t,w,s) $ N_REN2(N) ..                          Q_RN2s(N,t,w,s) =L= P_RNs2(N,t,w,s) * TAN ( ARCCOS(0.95));
*Active and reactive power relation for negative renewable generation
E_RENCOS12_S(N,t,w,s) $ N_REN2(N) ..                         Q_RN2s(N,t,w,s) =g= P_RNs2(N,t,w,s) * TAN ( ARCCOS(-0.95));
*---------------------------------------------
E_MARCOS(N,t,w) $ N_DG(N) ..                           Q_DG(N,t,w) =L= P_DG(N,t,w) * TAN( ARCCOS(0.95));
*Active & reactive power relation for negative convetional generation
E_MARCOS1(N,t,w) $ N_DG(N) ..                          Q_DG(N,t,w) =g= P_DG(N,t,w) * TAN( ARCCOS(-0.95));
*Active & reactive power relation for positive subsation
E_SUB(N,t,w) $(ord(N)=1) ..                            QSS(N,t,w) =L= PSS(N,t,w) * TAN( ARCCOS(0.95));
*Active & reactive power relation for negative subsation
E_SUB1(N,t,w)$(ord(N)=1) ..                            QSS(N,t,w) =g= PSS(N,t,w) * TAN( ARCCOS(-0.95));
*two stage--------------------------------------
*E_SUB_S(N,t,w,s) $(ord(N)=1) ..                            Q_Res_ups(N,t,w,s) =L= Res_ups(N,t,w,s) * TAN( ARCCOS(0.95));
*Active & reactive power relation for negative subsation
*E_SUB1_S(N,t,w,s)$(ord(N)=1) ..                            Q_Res_ups(N,t,w,s) =g= Res_ups(N,t,w,s) * TAN( ARCCOS(-0.95));
*Convetional active power generation limit

P_DG_limit(N,t,w)$ N_DG(N)..                           P_DG(N,t,w) =L= .7;

reserve_limit_ups(N,t,w,s)$(ord(N)=1)..                           Res_ups(N,t,w,s)=l=Res_up(N,t,w);
reserve_limit_dns(N,t,w,s)$(ord(N)=1)..                           Res_dns(N,t,w,s)=l=Res_dn(N,t,w);
reserve_limit_ups_L(N,t,w,s)$(ord(N)=1)..                            Res_up(N,t,w)+PSS(N,t,w)=l=4;
reserve_limit_dns_L(N,t,w,s)$(ord(N)=1)..                            -Res_dn(N,t,w)+PSS(N,t,w)=g=-4;

Qreserve_limit_ups(N,t,w,s)$(ord(N)=1)..                           Q_Res_ups(N,t,w,s)=l=Q_Res_up(N,t,w);
Qreserve_limit_dns(N,t,w,s)$(ord(N)=1)..                           Q_Res_dns(N,t,w,s)=l=Q_Res_dn(N,t,w);
Qreserve_limit_ups_L(N,t,w,s)$(ord(N)=1)..                            Q_Res_up(N,t,w)+QSS(N,t,w)=l=Res_up(N,t,w) * TAN( ARCCOS(0.95));
Qreserve_limit_dns_L(N,t,w,s)$(ord(N)=1)..                            -Q_Res_dn(N,t,w)+QSS(N,t,w)=g=Res_dn(N,t,w) * TAN( ARCCOS(-0.95));

DR_limit2(N,t)$N_DR(N)..     P_DR1(N,t)=l= .1*Pdem(N,'pt',t);
DR_limit3(N,t,s)$N_DR(N)..   P_DRs(N,t,s)=l= .1*Pdem(N,'pt',t);
DR_limit4(N,t,s)$N_DR(N)..   P_DRs(N,t,s)=l= P_DR1(N,t);
ESS_con1(t,N,s)$N_ES(N)..    Pd_s(t, N,s)=l= Pd(t, N);
ESS_con2(t,N,s)$N_ES(N)..    Pc_s(t, N,s)=l= Pc(t, N);
ESS_con3(t,N,s)$N_ES(N)..    SOC_s(t, N,s)=l= SOC(t, N);

VMAX=1;
VNOM=1;
VMIN=0.950;
VMIN2=VMIN*VMIN;
VMAX2=VMAX*VMAX;
IMIN2=0;
VNOM2=VNOM*VNOM;
*The maximum renewable limits are controlled by the number of turbines
DELTA_S(N,NP,t,w) = VNOM * LineData(N,NP,'IMAX')/card(R);
*______________________________________________
*node voltage limits p.u.
V0_2(N,t,w)=VNOM2;
V2.lo(N,t,w)=VMIN2;
V2.up(N,t,w)=VMAX2;
**PUNTO INICIAL (FLATSTART POINT)
V2.l(N,t,w)=VNOM2;
**PUNTO INICIAL (BARRA SLACK)
V2.fx(SL,t,w)=VNOM2;

V2_S.lo(N,t,w,s)=VMIN2;
V2_S.up(N,t,w,s)=VMAX2;
**PUNTO INICIAL (FLATSTART POINT)
V2_S.l(N,t,w,s)=VNOM2;
**PUNTO INICIAL (BARRA SLACK)
V2_S.fx(SL,t,w,s)=VNOM2;
*_______________________________________________
W_MAX2= VMAX2 - VMIN2;

*

MODEL totalcost
*$ontext
/
*E_FOBJ
Max_Profit
constES
constES_s
E_BAL_POT_P
E_BAL_POT_Q
E_BAL_POT_P_S
E_BAL_V
E_POT_S
E_POT_P
E_POT_Q
E_COTA_DELTA_P
E_COTA_DELTA_Q
E_INT_MAX
E_SWITCH
E_W_MAX
E_W_MIN
E_L_POT1
E_L_POT2
E_RENCOS
E_RENCOS1
E_MARCOS
E_MARCOS1
P_DG_limit
reserve_limit_ups

E_BAL_POT_Q_S,
E_BAL_V_S,
E_POT_S_SS,
E_POT_P_S,
E_POT_Q_S,
E_COTA_DELTA_P_S,
E_COTA_DELTA_Q_S,
E_INT_MAX_S,
E_W_MAX_S,
E_W_MIN_S,
E_L_POT1_S,
E_L_POT2_S,
E_RENCOS_S,
E_RENCOS1_S,
E_RENCOS2_S,
E_RENCOS12_S,
*E_SUB_S,
*E_SUB1_S,
*reserve_limit_dns_L,
reserve_limit_dns,
*reserve_limit_ups_L,

*Qreserve_limit_dns_L,
Qreserve_limit_dns,
*Qreserve_limit_ups_L,
Qreserve_limit_ups,

DR_limit2
DR_limit3
DR_limit4

ESS_con1
ESS_con2
ESS_con3

/

OPTION LIMROW =100000000;
*Tolerancia para la convergencia de la optimizaci?n con variables binarias
OPTION OPTCR =0.15;
* N?mero m?ximo de iteraciones que se le permite haceral optimizador
OPTION ITERLIM =1e8;
* N?mero m?ximo de iteraciones que se le permite hacer al optimizador
OPTION RESLIM =1e10;
*OPTION threads=0;
OPTION DECIMALS =8;
OPTION mip=cplex;
SOLVE  totalcost USING MIP MINIMIZING ProfitTest;
$ontext
*_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________
*Execute "xlstalk -S Outlets.xls"
Execute_Unload "Model.gdx" ProfitTest.l,I2.l,PSS.l,QSS.l,Q_RN.l,V2.l,P_DG.l,Q_DG.l,P_POS.l,P_NEG.l,Q_POS.l,Q_NEG.l,WV.l,var_switch;
Execute 'GDXXRW.EXE Model.gdx O=Outlets.xls SQ=N Var=ProfitTest.l Rng=ProfitTest!B1';
Execute 'GDXXRW.EXE Model.gdx O=Outlets.xls SQ=N Var=I2.l Rng=CurSQR!B1 Var=PSS.l Rng=SubActive!B1 Var=QSS.l Rng=SubReactive!B1 Var=V2.l Rng=NodeVolt!B1';
Execute 'GDXXRW.EXE Model.gdx O=Outlets.xls SQ=N Var=WV.l  Rng=W!B1 Var=P_POS.l Rng=P_POS!B1 Var=P_NEG.l Rng=P_NEG!B1 Var=Q_POS.l Rng=Q_POS!B1 Var=Q_NEG.l Rng=Q_NEG!B1';
Execute 'GDXXRW.EXE Model.gdx O=Outlets.xls SQ=N Var=P_DG.l Rng=P_DG!B2 Var=Q_DG.l Rng=Q_DG!B2 Var=var_switch Rng=Binary-Switch!B1';
Execute 'GDXXRW.EXE Model.gdx O=Outlets.xls SQ=N Var=Q_RN.l Rng=Q_RN!A2';
Execute 'GDXXRW.EXE Model.gdx O=Outlets.xls SQ=N Par=Windexp Rng=P_RN!A2';
Execute "=shellexecute Outlets.xls" ;
*____________________________________________________________________________________________________________________________________________________________________________________________________
**-------------------------------------------------------------------------------------------------------------**
$offtext
$LIBINCLUDE XLDUMP  PSS.l ECM1-res.xls Sheet1
$LIBINCLUDE XLDUMP  QSS.l ECM1-res.xls Sheet2
$LIBINCLUDE XLDUMP  res_up.l ECM1-res.xls Sheet3
$LIBINCLUDE XLDUMP  res_dn.l ECM1-res.xls Sheet4
$LIBINCLUDE XLDUMP  Q_res_up.l ECM1-res.xls Sheet5
$LIBINCLUDE XLDUMP  Q_res_dn.l ECM1-res.xls Sheet6
*$LIBINCLUDE XLDUMP  load_final.l DRresult_secondwork3.xls Sheet2

*$LIBINCLUDE XLDUMP  P_DR1.l secondwork3-severalDR.xls Sheet3
*$LIBINCLUDE XLDUMP  P_DR.l secondwork3-severalDR.xls Sheet4
*$LIBINCLUDE XLDUMP  TOTALDR_COST_T.l secondwork3-severalDR.xls Sheet1
*$LIBINCLUDE XLDUMP  TOTALDR_COST_N.l secondwork3-severalDR.xls Sheet2
