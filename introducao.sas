data mapa;
input x y id;
cards;
0 0 1
0 2 1
2 2 1
2 0 1
2 2 2
4 4 2
5 2 2
;
data dados;
input id valor $;
cards;
2 triângulo
1 quadrado
;

proc gmap data=dados map=mapa all;
id id;
choro valor;
run;
quit;

proc mapimport datafile='/home/alansilva0/my_shared_file_links/alansilva0/Estatistica_Espacial/RJ_Municipios_2022.shp'
out=rj;
run;

data dados_rj;
cd_mun='3300100';
v=2;
run;

proc gmap data=dados_rj map=rj all;
id cd_mun;
choro v;
run;
quit;

proc gmap data=dados_rj map=rj;
id cd_mun;
choro v;
run;
quit;

proc sort data=rj out=rj2;by x y;run;
proc gmap data=dados_rj map=rj2;
id cd_mun;
choro v;
run;
quit;

data rj3;set rj;
u=ranuni(2);
run;
proc sort data=rj3;by u;run;
proc gmap data=dados_rj map=rj3 all;
id cd_mun;
choro v;
run;
quit;

proc sql noprint;
select min(x) into:minx from rj;
select min(y) into:miny from rj;
select max(x) into:maxx from rj;
select max(y) into:maxy from rj;
quit;
%put &minx &miny &maxx &maxy;

data grid;
do x=&minx to &maxx by 0.01;
	do y=&miny to &maxy by 0.01;
		output;
	end;
end;
run;
data grid;set grid;
length function style $10. color $8.;
retain line 1 xsys ysys '2' hsys '3' color 'red';
function='label';text='U';position='5';style='marker';size=0.1;
run;

data dados_rj2;cd_mun='330000';v=2;run;

proc gmap data=dados_rj2 map=rj all;
id cd_mun;
choro v /anno=grid nolegend;
run;
quit;

proc ginside data=grid out=grid2 insideonly map=rj;
id cd_mun;
run;
proc gmap data=dados_rj2 map=rj all;
id cd_mun;
choro v /anno=grid2 nolegend;
run;
quit;

