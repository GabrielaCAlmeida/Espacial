data triangulo;
input id x y;
cards;
1 0 0
1 0.5 1
1 1 0
;
data a;id=2;v=2;run;
proc gmap data=a map=triangulo all;
id id;
choro v /nolegend;
run;
quit;
data anno;length text $2.;retain hsys '3' xsys ysys '2' when 'a';
x=0.5;y=0.3;
function='label';color='black';size=5;style='special';text='J';output;
run;
proc gmap data=a map=triangulo all;
id id;
choro v / anno=anno nolegend;
run;
quit;

%macro ginside(map,data,saida);
proc iml;
use &map;
read all var{x} into _p1;
read all var{y} into _p2;
p=_p1||_p2;
free _p1 _p2;

start verificapontointercarestadireita(p1,p2,ponto);
estaaesquerda=0;
interceptasegmento=0;
numerodearestasinterceptadas=0;
ymin=min(p1[2],p2[2]);
ymax=max(p1[2],p2[2]);
y=ponto[2];
if (p2[2]-p1[2])^= 0 then x=p1[1]+(y-p1[2])*(p2[1]-p1[1])/(p2[2]-p1[2]);
else x=p1[1];
if ponto[1]<x then estaaesquerda=1;
if y>ymin & y<=ymax then interceptasegmento=1;
if estaaesquerda=1 & interceptasegmento=1 then numerodearestasinterceptadas=1;
return(numerodearestasinterceptadas);
finish verificapontointercarestadireita;

start aresta(p,ponto);
p1=p[1,];
p2=p[2,];
ponto=ponto;
cond=verificapontointercarestadireita(p1,p2,ponto);
return(cond);
finish aresta;

start isdentrodopoligono(ponto) global(p);
_tipo_=0;
ponto=ponto;
numerodearestasinterceptadas=0;
do i=1 to nrow(p)-1;
z=aresta(p[i:i+1,],ponto);
if z=1 then numerodearestasinterceptadas=numerodearestasinterceptadas+1;
*print numerodearestasinterceptadas;
end;
if mod(numerodearestasinterceptadas,2)=1 then _tipo_=1;
result=ponto||_tipo_;
append from result;
finish isdentrodopoligono;
use &data;
read all into ponto;
result=j(1,3,0);
create &saida from result[colname={"x" "y" "_tipo_"}];
do j=1 to nrow(ponto);
run isdentrodopoligono(ponto[j,]);
end;
close &saida;
quit;
%mend ginside;
data dados;
do x=0 to 1 by 0.005;
do y=0 to 1 by 0.005;
output;
end;
end;
run;
%ginside(triangulo,dados,saida);

data anno;length text $2.;set saida;retain hsys '3' xsys ysys '2' when 'a';
function='label';text='U';position='5';style='marker';size=2;
if _tipo_=1 then color="red";
run;
proc gmap data=a map=triangulo all;
id id;
choro v / anno=anno nolegend;
run;
quit;

proc mapimport out=sao_carlos datafile='/home/alansilva0/my_shared_file_links/alansilva0/Estatistica_Espacial/sao_carlos.shp';run;
data a;segment=2;v=2;run;
proc gmap data=a map=sao_carlos(where=(SPRPERIMET>10000)) all;
id segment;
choro v / nolegend;
run;
quit;
proc sql noprint;
select min(x) into:minx from sao_carlos (where=(SPRPERIMET>10000));
select min(y) into:miny from sao_carlos (where=(SPRPERIMET>10000));
select max(x) into:maxx from sao_carlos (where=(SPRPERIMET>10000));
select max(y) into:maxy from sao_carlos (where=(SPRPERIMET>10000));
quit;
%put minx=&minx maxx=&maxx miny=&miny maxy=&maxy;
data nums;
do x=&minx to &maxx by 100;
do y=&miny to &maxy by 100;
output;
end;
end;
run;
data Sao_carlos1;set sao_carlos (where=(SPRPERIMET>10000));run;
%ginside(Sao_carlos1,nums,saida);
data anno;length text $2.;set saida;
retain hsys '3' xsys ysys '2' when 'a';
function='label';text='U';position='5';style='marker';size=.6;
if _tipo_=1 then color="red";
run;
proc gmap data=a map=sao_carlos(where=(SPRPERIMET>10000)) all;
id segment;
choro v / anno=anno nolegend;
run;
quit;

/*********** calculando percentual de pontos dentro de um poligono ***/

proc iml;
use sao_carlos1;
read all var{x} into _p1;
read all var{y} into _p2;
p=_p1||_p2;
free _p1 _p2;

a1=min(p[,1])||min(p[,2]);
a2=min(p[,1])||max(p[,2]);
a3=max(p[,1])||max(p[,2]);
a4=max(p[,1])||min(p[,2]);

use saida;
read all into pontos;
totaldepontos=nrow(pontos);
pontosdentro=pontos[+,3];
ladoquadrado=a2[2]-a1[2];
relacaoentreareas=(pontosdentro/totaldepontos);
area=relacaoentreareas*(ladoquadrado*2)**2;
print relacaoentreareas area;
quit;




data quadrado;
input id x y;
cards;
1 0 0
1 0 1
1 1 1
1 1 0
;
data a;id=2;v=2;run;
proc gmap data=a map=quadrado all;
id id;
choro v / nolegend;
run;
quit;

proc iml;
use quadrado;
read all var{x} into _p1;
read all var{y} into _p2;
p=_p1||_p2;
free _p1 _p2;

start verificapontointercarestadireita(p1,p2,ponto);
estaaesquerda=0;
interceptasegmento=0;
numerodearestasinterceptadas=0;
ymin=min(p1[2],p2[2]);
ymax=max(p1[2],p2[2]);
y=ponto[2];
if (p2[2]-p1[2])^= 0 then x=p1[1]+(y-p1[2])*(p2[1]-p1[1])/(p2[2]-p1[2]);
else x=p1[1];
if ponto[1]<x then estaaesquerda=1;
if y>ymin & y<=ymax then interceptasegmento=1;
if estaaesquerda=1 & interceptasegmento=1 then numerodearestasinterceptadas=1;
return(numerodearestasinterceptadas);
finish verificapontointercarestadireita;

start aresta(p,ponto);
p1=p[1,];
p2=p[2,];
ponto=ponto;
cond=verificapontointercarestadireita(p1,p2,ponto);
return(cond);
finish aresta;

start isdentrodopoligono(ponto) global(p);
_tipo_=0;
numerodearestasinterceptadas=0;
do i=1 to nrow(p)-1;
z=aresta(p[i:i+1,],ponto);
if z=1 then numerodearestasinterceptadas=numerodearestasinterceptadas+1;
*print numerodearestasinterceptadas;
end;
if mod(numerodearestasinterceptadas,2)=1 then _tipo_=1;
return(_tipo_);
finish isdentrodopoligono;

pontosdentro=0;
totaldepontos=1000000;
ladoquadrado=2;

do j=1 to totaldepontos;
ponto=ranuni(0)*ladoquadrado||ranuni(0)*ladoquadrado;
tp=isdentrodopoligono(ponto);
if tp=1 then pontosdentro=pontosdentro+1;
end;
relacaoentreareas=1*pontosdentro/totaldepontos;
area=relacaoentreareas*(ladoquadrado*ladoquadrado);
print relacaoentreareas area;
quit;

data circulo;
id=1;
do x=-1 to 1 by 0.01;
if x>=1 then y=0;else y=sqrt(1-x**2);
output;
end;
run;
proc sort data=circulo out=circulo1;by descending x;run;
data circulo;set circulo circulo1(in=a);
if a then y=-y;
run;
data a;id=2;v=2;run;
proc gmap data=a map=circulo all;
id id;
choro v /nolegend;
run;
quit;

data dados;
do x=-1 to 1 by 0.05;
do y=-1 to 1 by 0.05;
output;
end;
end;
run;
%ginside(circulo,dados,saida);

data anno;length text $2.;set saida;retain hsys '3' xsys ysys '2' when 'a';
function='label';text='U';position='6';style='marker';size=2;
if _tipo_=1 then color="red";
run;
proc gmap data=a map=circulo all;
id id;
choro v / anno=anno nolegend;
run;
quit;

proc iml;
use circulo;
read all var{x} into _p1;
read all var{y} into _p2;
p=_p1||_p2;
free _p1 _p2;

start verificapontointercarestadireita(p1,p2,ponto);
estaaesquerda=0;
interceptasegmento=0;
numerodearestasinterceptadas=0;
ymin=min(p1[2],p2[2]);
ymax=max(p1[2],p2[2]);
y=ponto[2];
if (p2[2]-p1[2])^= 0 then x=p1[1]+(y-p1[2])*(p2[1]-p1[1])/(p2[2]-p1[2]);
else x=p1[1];
if ponto[1]<x then estaaesquerda=1;
if y>ymin & y<=ymax then interceptasegmento=1;
if estaaesquerda=1 & interceptasegmento=1 then numerodearestasinterceptadas=1;
return(numerodearestasinterceptadas);
finish verificapontointercarestadireita;

start aresta(p,ponto);
p1=p[1,];
p2=p[2,];
ponto=ponto;
cond=verificapontointercarestadireita(p1,p2,ponto);
return(cond);
finish aresta;

start isdentrodopoligono(ponto) global(p);
_tipo_=0;
numerodearestasinterceptadas=0;
do i=1 to nrow(p)-1;
z=aresta(p[i:i+1,],ponto);
if z=1 then numerodearestasinterceptadas=numerodearestasinterceptadas+1;
*print numerodearestasinterceptadas;
end;
if mod(numerodearestasinterceptadas,2)=1 then _tipo_=1;
return(_tipo_);
finish isdentrodopoligono;

pontosdentro=0;
totaldepontos=1000;
ladoquadrado=2;

do j=1 to totaldepontos;
ponto=ranuni(0)*ladoquadrado||ranuni(0)*ladoquadrado;
tp=isdentrodopoligono(ponto);
if tp=1 then pontosdentro=pontosdentro+1;
end;
relacaoentreareas=1*pontosdentro/totaldepontos;
area=relacaoentreareas*(arcos(-1)*ladoquadrado**2);
print relacaoentreareas area;
quit;
proc iml;
pi=arcos(-1);
print pi;
quit;


/************** gerando 2 circulos centrados em pontos diferentes ******/

data circulo;
id=1;
do x=-1 to 1 by 0.01;
if x>=1 then y=0;else y=sqrt(1-x**2);
output;
end;
run;
proc sort data=circulo out=circulo1;by descending x;run;
data circulo;set circulo circulo1(in=a);
if a then y=-y;
run;
data circulo2;
id=2;
do x=1 to 3 by 0.01;
if x>=3 then y=2;else y=sqrt(1**2-(x-2)**2)+2;
output;
end;
run;
data circulo3;
id=2;
do x=3 to 1 by -0.01;
if x>=3 then y=2;else y=-sqrt(1**2-(x-2)**2)+2;
output;
end;
run;
data circulo2;set circulo2 circulo3;run;
data circulo;set circulo circulo2;run;
data a;id=99;v=2;run;
proc gmap data=a map=circulo all;
id id;
choro v /nolegend;
run;
quit;

