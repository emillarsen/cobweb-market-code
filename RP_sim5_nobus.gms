option miqcp=cplex;
*avoid limit on iterations
option iterlim=999999999;
*timelimit for solver in sec.
option reslim=120;
*gap tolerance
option optcr=0.1;

sets
g        generating units
t        time steps
l        loads                   /l1*l2/
;
alias(t, tt)
;

parameters
Psch(g,t)        set point of generators
PupR(g,t)        price for up regulation
QupR(g,t)        quantity of up regulation available
PdoR(g,t)        price for down regulation
QdoR(g,t)        quantity of down regulation available
ramp             ramping of generators
alpha(l,t)       price elasticity of load
SpotPrice(t)     day ahead price
WS_cost          wind shedding cost
LS_cost          load shedding cost
FC_cost          frequency control cost
load_bid(l,t)    two loads - first fixed second flexible
LoadFlex(l)      bus information about loads /l1 0, l2 1/
balance(t)       balancing power
change_l_max(l,t)     what is max flexibility in MW
change_l_min(l,t)     what is min flexibility in MW
wind_bid(t)      wind bid
bin(t)
ChangeL_history(l) energy used by load in last three hours
Pur_intial(g)
Pdr_intial(g)
;

$gdxin InputDR
$load t, g
$load Psch, PupR, QupR, PdoR, QdoR, ramp
$load alpha, SpotPrice, WS_cost, LS_cost, FC_cost
$load load_bid, balance, wind_bid, bin, ChangeL_history
$load change_l_max, change_l_min, Pur_intial, Pdr_intial
$gdxin

display Psch, PupR, QupR, PdoR, QdoR, ramp, alpha, SpotPrice, bin
display WS_cost, LS_cost, FC_cost, load_bid, wind_bid, change_l_max, change_l_min, balance

positive variables
Pur(g,t)         up regulation provided per generator
Pdr(g,t)         down regulation provided  per generator
Shed(l,t)        load shedding
SpillW(t)        wind spill
PG(g,t)          production total (day ahead bid + regulation)
PL(l,t)          total load (inflexible + flexible)
PW(t)            total wind power injected
;

binary variables
x(g,t)           up regulation delivered
z(g,t)           down regulation delivered
k(t)             up regulation global for t
y(t)             down regulation global for t
m(g,t)           start up of up regulation
n(g,t)           start up of down regulation
;

free variables
TotalCost
GenCost(g,t)
ChangeL(l,t)
RT(l,t)            real-time price for price-elastic consumers
gradient_up(g,t)
gradient_dr(g,t)
balance_sanity(t)
new_check(t)
cost_over_time(t)
cost_over_time_with_wind(t)
sw
sw_with_wind_load_shed(t)
;

equations
objective        objective function
social_welfare_time
social_welfare_time2
sw_nowindorloadshed
sw_wind_and_load

con1             cost for each conventional power plant
con2             electricity produced by each unit
con3             load consumed
con4             wind injected
con5             balance constraint
con6             up regulation gradient change every 30 minutes
con7             up regulation gradient fixed inbetween
con8             plateu in up regulation for 15 minutes
con9             down regulation gradient change every 30 minutes
con10            down regulation gradient fixed inbetween
con11            plateu in down regulation for 15 minutes
con12            change in up regulation
con13            change in down regulation
con14            maximum up regulation quantity
con15            maximum down regulation quantity
con16            minimum up regulation quantity (with ramp rate built in)
con17            minimum down regulation quantity (with ramp rate built in)
con18            up regulation off when quantity available is zero
con19            down regulation off when quantity available is zero
*con20            bid activation (start-up) variable (up)
*con21            bid activation (start-up) variable (down)
*con22            minimum on time for up regulation
*con23            minimum on time for down regulation
con24            is up regulation being delivered?
con25            is down regulation being delivered?
con26            stop up and down regulation from being active at the same time
con27            real-time price equation
con28            max load consumption
con29            min load consumption
con30            load shifting
con31
con32
;

objective..                      TotalCost =e= sum((l,t),(SpotPrice(t)*ChangeL(l,t))-(0.5*alpha(l,t)*ChangeL(l,t)*ChangeL(l,t))) -
                                 sum((g,t),GenCost(g,t)) - sum(t,WS_cost*SpillW(t)) - sum((l,t),(LS_cost*Shed(l,t)));

social_welfare_time(t)..         cost_over_time(t) =e= sum(l,(SpotPrice(t)*ChangeL(l,t))-(0.5*alpha(l,t)*ChangeL(l,t)*ChangeL(l,t)))
                                 - sum(g,GenCost(g,t));

sw_nowindorloadshed..            sw =e= sum((l,t), (SpotPrice(t)*ChangeL(l,t))-(0.5*alpha(l,t)*ChangeL(l,t)*ChangeL(l,t)))
                                 - sum((g,t), GenCost(g,t));

social_welfare_time2(t)..        cost_over_time_with_wind(t) =e= sum(g,GenCost(g,t)) +
                                 WS_cost*SpillW(t) + sum(l,(LS_cost*Shed(l,t))) +
                                 sum(l,(-SpotPrice(t)*ChangeL(l,t))+(0.5*alpha(l,t)*ChangeL(l,t)*ChangeL(l,t)));

sw_wind_and_load(t)..            sw_with_wind_load_shed(t) =e= sum(l,(SpotPrice(t)*ChangeL(l,t))-(0.5*alpha(l,t)*ChangeL(l,t)*ChangeL(l,t))) -
                                 sum(g,GenCost(g,t)) - WS_cost*SpillW(t) - sum(l,(LS_cost*Shed(l,t)));

con1(g,t)..                      GenCost(g,t) =e= PupR(g,t)*Pur(g,t) - PdoR(g,t)*Pdr(g,t);
* + FC_cost*(Pur(g,t)+Pdr(g,t))
con2(g,t)..                      PG(g,t) =e= Psch(g,t) + Pur(g,t) - Pdr(g,t);
con3(l,t)..                      PL(l,t) =e= load_bid(l,t) - Shed(l,t)  + ChangeL(l,t);
con4(t)..                        PW(t) =e= wind_bid(t) - SpillW(t);
con5(t)..                        sum(g,PG(g,t)) - balance(t) + PW(t) =e=  sum(l,PL(l,t));

* generator constraints

con6(g,t)$(ord(t)=1)..          gradient_up(g,t) =e= Pur(g,t) - Pur_intial(g);
con7(g,t)$(ord(t)=1)..          gradient_dr(g,t) =e= Pdr(g,t) - Pdr_intial(g);

con8(g,t)$(bin(t) eq 0 and ord(t) > 1).. gradient_up(g,t) =e= gradient_up(g,t-1);
con9(g,t)$(bin(t) eq 0 and ord(t) > 1).. gradient_dr(g,t) =e= gradient_dr(g,t-1);

*con10(g,t)$(bin(t) eq 1 and ord(t) > 1).. gradient_up(g,t) =e= Pur(g,t)-Pur(g,t-1);
*con11(g,t)$(bin(t) eq 1 and ord(t) > 1).. gradient_dr(g,t) =e= Pdr(g,t)-Pdr(g,t-1);

con10(g,t)$(ord(t) > 1).. gradient_up(g,t) =e= Pur(g,t)-Pur(g,t-1);
con11(g,t)$(ord(t) > 1).. gradient_dr(g,t) =e= Pdr(g,t)-Pdr(g,t-1);

con12(g,t)$(bin(t) eq 2)..       gradient_up(g,t) =e= 0;
con13(g,t)$(bin(t) eq 2)..       gradient_dr(g,t) =e= 0;

con14(g,t)..                     Pur(g,t) =l= x(g,t)*QupR(g,t);
con15(g,t)..                     Pdr(g,t) =l= z(g,t)*QdoR(g,t);
con16(g,t)$(bin(t) eq 2)..       Pur(g,t) =g= x(g,t)*ramp*QupR(g,t);
con17(g,t)$(bin(t) eq 2)..       Pdr(g,t) =g= z(g,t)*ramp*QdoR(g,t);

con18(g,t)$(QdoR(g,t) eq 0)..    z(g,t) =e= 0;
con19(g,t)$(QupR(g,t) eq 0)..    x(g,t) =e= 0;

* minimum on time
*con20(g,t)$(ord(t)>1)..          m(g,t) =g= x(g,t)-x(g,t-1);
*con21(g,t)$(ord(t)>1)..          n(g,t) =g= z(g,t)-z(g,t-1);
*con22(g,t)..                     sum((tt)$(ord(tt) >= ord(t) and ord(tt) < ord(t) + 6), x(g,tt)) =g= 6 * m(g,t);
*con23(g,t)..                     sum((tt)$(ord(tt) >= ord(t) and ord(tt) < ord(t) + 6), z(g,tt)) =g= 6 * n(g,t);

* stop up and down bids being activated at the same time
con24(g,t)..                      Pur(g,t) =l= k(t)*QupR(g,t);
con25(g,t)..                      Pdr(g,t) =l= y(t)*QdoR(g,t);
con26(t)$(bin(t) eq 2)..          k(t) + y(t) =l= 1;

* load constraints
con27(l,t)$(LoadFlex(l) eq 1)..  -alpha(l,t)*ChangeL(l,t) =e= RT(l,t)-SpotPrice(t);
con28(l,t)..                     ChangeL(l,t) =l= change_l_max(l,t);
con29(l,t)..                     ChangeL(l,t) =g= change_l_min(l,t);
con30(l)..                       sum(t,ChangeL(l,t)) + ChangeL_history(l) =e= 0;

con31(g,t)..                     Pur(g,t) =g= x(g,t)*ramp*0.3*QupR(g,t);
con32(g,t)..                     Pdr(g,t) =g= z(g,t)*ramp*0.3*QdoR(g,t);


model NewModel /all/;

solve NewModel maximizing TotalCost using miqcp;

display con5.m;

Parameter
marginal_cost(t);
marginal_cost(t) = con5.m(t);

* to see minimum on times etc.
display Pur.l, x.l, Pdr.l, z.l;

display QupR;

* to see shedding
display Shed.l, SpillW.l;

* to see production and consumption totals
display PG.l, PL.l, PW.l, ChangeL.l;

* to see costs
display TotalCost.l, marginal_cost, SpotPrice, sw.l;

execute_unload 'OutputDR', gradient_up, gradient_dr, alpha, Pur.l, Pdr.l, Shed.l, SpillW.l, PG.l, PL.l, PW.l, RT.l, TotalCost.l, ChangeL.l, marginal_cost, Psch, SpotPrice, balance_sanity.l, new_check.l, load_bid, wind_bid, x.l, z.l, cost_over_time.l, cost_over_time_with_wind.l, sw_with_wind_load_shed.l, cost_over_time_with_wind.l;
