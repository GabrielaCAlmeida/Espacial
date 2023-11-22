%macro lisa(tab,var,matrizw,tipo=);
proc iml;
use &tab var{&var};
read all into y;
use &matrizw;
read all into w;
n=nrow(y);
yb=sum(y)/n;
yi=y-yb;
yi2=yi#yi;
wj=j(n,1,0);
do i=1 to n;
do j=1 to n;
wj[i]=wj[i]+w[i,j]#yi[j];
end;
end;
L=(yi#wj)/(sum(yi2)/n);
*use &matrizw;
*read all into w;
wi=j(n,1,0);                                                                                                                         
do i=1 to n;                                                                                                                           
do j=1 to n;                                                                                                                           
wi[i]=wi[i]+w[i,j];                                                                                                                 
end;                                                                                                                                    
end;
Ei=-wi/(n-1);
wi2=j(n,1,0);                                                                                                                         
do i=1 to n;                                                                                                                           
do j=1 to n;                                                                                                                           
wi2[i]=wi2[i]+w[i,j]#w[i,j];                                                                                                                 
end;                                                                                                                                    
end;
m2=sum(yi#yi)/n;
m4=sum(yi#yi#yi#yi)/n;
m22=m2*m2;
b2=m4/m22;
Vari=(wi2*(n-b2)/(n-1))+(wi*(2*b2-n)/((n-1)*(n-2)))-(wi#wi/((n-1)*(n-1)));
z=(l-ei)/sqrt(vari);
create L var{L z};
append;
quit;
data l;
merge l &tab;
prob=2*(1-probnorm(abs(z)));
run;
data l1;
set l;
if prob>0.05 then sig95='N�o Significativo (95%)';
if prob>0.01 then sig99='N�o Significativo (99%)';
if prob>0.001 then sig999='N�o Significativo (99.9%)';
run;
proc greplay igout=work.gseg nofs;
delete lis_fin lisa br_lisa;
run;
quit;
proc sql noprint;
create table l2 as
select distinct I from l1
where sig95='';
select n(i) into:nc from l2;
quit;
data l2;
set l2;
if i='Low-Low' then c='pink';
if i='High-High' then c='red';
if i='High-Low' then c='blue';
if i='Low-High' then c='vpab';
call symput('c'||trim(left(_n_)),c);
run;
goptions reset=all;
%do j=1 %to &nc;
pattern&j c=&&c&j;
%end;
goptions reset=footnote;
goptions ftitle='Verdana' ftext='Verdana';
title1 "Moran Map (95%)";
title2 "vari�vel &var";
%if %upcase(&tipo)=MICRO %then %do;
%let t=brasilmicro;
%let t1=microcod;
%end;
%else %if %upcase(&tipo)=MUN %then %do;
%let t=brasil;
%let t1=codigo;
%end;
%else %if %upcase(&tipo)=MESO %then %do;
%let t=brasilmeso;
%let t1=mesocod;
%end;
%else %do;
%let t=&tipo;
%let t1=&codigo;
%end;
proc gmap data=l1 map=&t all;
id &t1;
choro i / coutline=black name='lisa';
where sig95='';
run;
quit;
goptions ftitle='Verdana' ftext='Verdana';
pattern v=s c=white repeat=4;
title1 "Moran Map (95%)";
title2 "vari�vel &var";
proc gmap data=l1
%if %upcase(&tipo)=MESO or %upcase(&tipo)=MUN or %upcase(&tipo)=MICRO %then map=mapas3.brasilest;
%else map=&t; all;
id &t1;
choro i / cempty=black name='br_lisa';
where sig95='';
run;
proc greplay igout=work.gseg gout=work.gseg nofs tc=sashelp.templt template=whole;
treplay 1:br_lisa 1:lisa name='lis_fin';
run;
quit;
%mend lisa;