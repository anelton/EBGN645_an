set i "parcels" /A, B, C/ ;
set t "years"   /1*2/ ;

alias (i, j) ;
alias (t, tt) ;

parameter
wfp(i) "annual prob. a wildfire starts in parcel i"  /A 0.25,   B 0.18,   C 0.15/,
wfd(i) "damage if wildfire affects parcel i"         /A 400000, B 475000, C 600000/,
wfr(i) "reduction in prob. after treatment"          /A 0.90,   B 0.92,   C 0.88/,
c(i)   "cost of treating parcel i"                   /A 25000,  B 40000,  C 35000/,
trd(i) "damage from prescribed burning"              /A 2000,   B 4000,   C 2000/,
mb(t)  "max treatments per year"                     /1 1, 2 1/ ;

scalar tb "treatment budget across both years" /75000/ ;

* wfs(j,i): prob. a fire starting in j spreads directly to i (before treatment)
table wfs(j,i)
        A       B       C
A               0.45    0.30
B       0.40            0.35
C       0.35    0.35          ;

* total prob. parcel i is affected: own ignition + fire spreading in from j
parameter pa(i) "total prob. parcel i is affected by wildfire" ;
pa(i) = wfp(i) + sum(j$wfs(j,i), wfp(j) * wfs(j,i)) ;

binary variable X(i,t) "1 if parcel i is treated in year t" ;
positive variable Y(i,t) "1 if parcel i treated in current or prior year" ;
variable TC "total costs" ;

equation objfn, budget, maxtreat(t), once(i), track(i,t) ;

objfn.. TC =e=
    sum((i,t), pa(i) * (1 - wfr(i) * Y(i,t)) * wfd(i))
  + sum((i,t), (c(i) + trd(i)) * X(i,t)) ;

budget..      tb =g= sum((i,t), c(i) * X(i,t)) ;
maxtreat(t).. mb(t) =g= sum(i, X(i,t)) ;
once(i)..     sum(t, X(i,t)) =l= 1 ;
track(i,t)..  Y(i,t) =e= sum(tt$(ord(tt) <= ord(t)), X(i,tt)) ;

model wildfire /all/ ;
option optcr = 0 ;

solve wildfire using MIP minimizing TC ;
parameter rep ;
rep("Baseline","X",i,t) = X.l(i,t) ;
rep("Baseline","Y",i,t) = Y.l(i,t) ;
rep("Baseline","TC","total","-") = TC.l ;

execute_unload 'alldata_wildfire.gdx' ;



* example scenario: allow two treatments per year
mb(t) = 2 ;

solve wildfire using MIP minimizing TC ;
rep("Double Treatment","X",i,t) = X.l(i,t) ;
rep("Double Treatment","Y",i,t) = Y.l(i,t) ;
rep("Double Treatment","TC","total","-") = TC.l ;

display pa, rep ;

execute_unload 'alldata_wildfire.gdx' ;
