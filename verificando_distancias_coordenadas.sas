proc mapimport out=sao_carlos 
datafile='/home/alansilva0/my_shared_file_links/alansilva0/Estatistica_Espacial/sao_carlos.shp';run;
proc mapimport out=sao_carlos_pt 
datafile='/home/alansilva0/my_shared_file_links/alansilva0/Estatistica_Espacial/sao_carlos_pt.shp';run;

data a;
SEGMENT=0;
v=2;
run;
proc gmap data=a map=sao_carlos(where=(SPRPERIMET>1000)) all;
id SEGMENT;
choro v /nolegend;
run;
quit;


data anno;set Sao_carlos_pt;
retain hsys '3' xsys ysys "2";
function='label';text='J';color='orange';style='special';size=4;when='a';
run;
proc gmap data=a map=sao_carlos(where=(SPRPERIMET>1000)) all 
anno=anno;
id SEGMENT;
choro v /nolegend;
run;
quit;

proc iml;
use sao_carlos_pt;
read all var{x y} into COORD;
n=nrow(coord[,1]);
d=j(1,3,0);                                                                                                                             
nome={"idi" "idj" "d"};                                                                                                                 
create dist from d[colname=nome];                                                                                                       
do i=1 to n;
	do j=i+1 to n;
	if abs(coord[,1])<180 then do;
		dif=abs(COORD[i,1]-COORD[j,1]);
		raio=arcos(-1)/180;                                                                                                              
		if dif=0 then arco=0;
		else
		/* Law of Cosines */  
		arco=arcos(sin(COORD[i,2]*raio)*sin(COORD[j,2]*raio)+cos(COORD[i,2]*raio)*
cos(COORD[j,2]*raio)*cos(dif*raio));
		d[1]=i;  
		d[2]=j;
		d[3]=arco*6371 /*Earth's Radius = 6371 (aproximately)*/;
		append from d;
	end;
	else do;
		d[1]=i;  
		d[2]=j;
		d[3]=sqrt((COORD[i,1]-COORD[j,1])**2+(COORD[i,2]-COORD[j,2])**2);   
		append from d;
	end;
	end;
end;
close dist;
quit;

/*** calculando distancia pela funcao SAS ***/
proc iml;
use sao_carlos_pt;
read all var{x y} into COORD;
d2 = distance(COORD, "L2");
print (d2[1,2]);
create distance from d2;
append from d2;
n=nrow(coord)**2;
d=shape(d2,n,1);
create dist2 from d;
append from d;
quit;

proc sort data=dist out=dist2;by idi d;run;

data a;segment=2;v=2;run;
%macro m(dist);
data _null_;set Sao_carlos_pt;
call symput('x'||trim(left(_n_)),x);
call symput('y'||trim(left(_n_)),y);
run;
%put &x1 &y1;
data circulo;
id=1;
do x=&x1-&dist to &x1+&dist by 100;
if x>=&x1+&dist then y=&y1;else y=sqrt(&dist**2-(x-&x1)**2)+&y1;
output;
end;
run;
data circulo2;
id=1;
do x=&x1+&dist to &x1-&dist by -100;
if x>=&x1+&dist then y=&y1;else y=-sqrt(&dist**2-(x-&x1)**2)+&y1;
output;
end;
run;
data circulo;set circulo circulo2;run;
data circulo;set circulo;x1=lag1(x);y1=lag1(y);if x1=. then delete;run;
proc sql noprint;select count(*) into:max from circulo;quit;%put &max;
data circulo_1;set circulo;if _n_=&max;x1=x;y1=y;run;
data circulo_2;set circulo;if _n_=1;x1=x;y1=y;run;
data circulo_2;merge circulo_2 circulo_1(drop=x1 y1);run;
data circulo;set circulo circulo_2;run;
data anno1;length function style color $8 position $1 text $20;
retain xsys ysys '2' hsys '3' when 'a';
set circulo;
style='music';text='A';color='orange';size=0.2;position='5';line=1;
function='move';x=x;y=y;output;
function='draw';x=x1;y=y1;output;
run;

goptions ftitle='Verdana';
data Sao_carlos_pt;set Sao_carlos_pt;call symput('sec'||trim(left(FID)),trim(left(FID)));run;
%do i=1 %to 1;
data annoidi;set Sao_carlos_pt(where=(FID=&i));
retain hsys '3' xsys ysys "2";
function='label';text='J';color='orange';style='special';size=2;when='a';
run;
data annoidi2;set dist(where=(idi=&i and d<=&dist));
retain hsys '3' xsys ysys "2";
function='label';text='J';color='black';style='special';size=3;when='a';
run;
data annoidi2;merge annoidi2(in=a rename=idj=fid) Sao_carlos_pt;by fid;if a;run;
data annoidi;set annoidi annoidi2;run;
title "Vizinhos do Ponto &&sec&i - com raio &dist";
proc gmap data=a map=Sao_carlos(where=(SPRPERIMET>10000)) all anno=annoidi;
id segment;
choro v /discrete nolegend anno=anno1;
run;
quit;
%end;
%mend m;
%m(3000);

/*********** usando distancia lat-long **************/

proc iml;
use sao_carlos_pt;
read all var{x y} into COORD;
n=nrow(coord[,1]);
d=j(1,3,0);                                                                                                                             
nome={"idi" "idj" "d"};                                                                                                                 
create dist from d[colname=nome];                                                                                                       
do i=1 to n;
	do j=i+1 to n;
		dif=abs(COORD[i,1]-COORD[j,1]);
		raio=arcos(-1)/180;                                                                                                              
		if dif=0 then arco=0;
		else
		/* Law of Cosines */  
		arco=arcos(sin(COORD[i,2]*raio)*sin(COORD[j,2]*raio)+cos(COORD[i,2]*raio)*cos(COORD[j,2]*raio)*cos(dif*raio));
		d[1]=i;  
		d[2]=j;
		d[3]=arco*6371 /*Earth's Radius = 6371 (aproximately)*/;
		append from d;
	end;
end;
close dist;
quit;

data a;segment=2;v=2;run;
%macro m(dist);
goptions ftitle='Verdana';
data Sao_carlos_pt;set Sao_carlos_pt;call symput('sec'||trim(left(FID)),trim(left(FID)));run;
%do i=1 %to 1;
data annoidi;set Sao_carlos_pt(where=(FID=&i));
retain hsys '3' xsys ysys "2";
function='label';text='J';color='orange';style='special';size=2;when='a';
run;
data annoidi2;set dist(where=(idi=&i and d<=&dist));
retain hsys '3' xsys ysys "2";
function='label';text='J';color='black';style='special';size=3;when='a';
run;
data annoidi2;merge annoidi2(in=a rename=idj=fid) Sao_carlos_pt;by fid;if a;run;
data annoidi;set annoidi annoidi2;run;
title "Vizinhos do Ponto &&sec&i - com raio &dist";
proc gmap data=a map=Sao_carlos(where=(SPRPERIMET>10000)) all anno=annoidi;
id segment;
choro v /discrete nolegend;
run;
quit;
%end;
%mend m;
%m(1000);

/***************************** importando arquivo coordenada lat-long *******************/

%include '/home/alansilva0/my_shared_file_links/alansilva0/Estatistica_Espacial/lat_long_utm.sas';
%conversao(sao_carlos_pt,0,22);
data sao_carlos_pontos;merge sao_carlos_pt conversao(rename=(lat_y=y long_x=x));run;
%conversao(sao_carlos,0,22);
data sao_carlos2;merge sao_carlos conversao(rename=(lat_y=y long_x=x));
if _n_>=182 then delete;run;

data a;
SEGMENT=0;
v=2;
run;
proc gmap data=a map=sao_carlos2(where=(SPRPERIMET>1000)) all;
id SEGMENT;
choro v /nolegend;
run;
quit;

data anno;set Sao_carlos_pontos;
retain hsys '3' xsys ysys "2";
function='label';text='J';color='orange';style='special';size=4;when='a';
run;
proc gmap data=a map=sao_carlos2(where=(SPRPERIMET>1000)) all anno=anno;
id SEGMENT;
choro v /nolegend;
run;
quit;

proc sgmap plotdata=sao_carlos_pontos mapdata=sao_carlos2(where=(SPRPERIMET>1000));
  choromap / mapid=SEGMENT lineattrs=(color=gray77);
  bubble x=x y=y size=AVG_Z /
                             fillattrs=(color=red
                                        transparency=0.5)
                             datalabel=fid 
                             datalabelpos=topright 
                             datalabelattrs=(color=black size=14 style=normal weight=bold);
run;
proc sgmap plotdata=sao_carlos_pontos mapdata=sao_carlos2(where=(SPRPERIMET>1000));
  choromap / mapid=SEGMENT lineattrs=(color=gray77);
  bubble x=x y=y size=AVG_Z /
                             colorresponse=AVG_Z 
                             datalabel=fid 
                             datalabelpos=topright 
                             datalabelattrs=(color=black size=14 style=normal weight=bold)
 							 colormodel=(yellow orange red) name='bub';
   gradlegend 'bub' / title='' position=right;
run;
proc sgmap plotdata=sao_carlos_pontos mapdata=sao_carlos2(where=(SPRPERIMET>1000));
  choromap / mapid=SEGMENT lineattrs=(color=gray77);
  scatter x=x y=y /          datalabel=fid 
                             datalabelpos=topright 
                             datalabelattrs=(color=black size=14 style=normal weight=bold);
run;
proc sgmap plotdata=sao_carlos_pontos;                     /* Bubble locations and labels data    */
  openstreetmap;                                           /* Use OSM base map                    */
  bubble x=x y=y size=AVG_Z /                   /* Bubble size based on number of beds */
                             fillattrs=(color=red          /* Bubble color                        */
                                        transparency=0.5)  /* Let base map show though bubble     */
                             datalabel=fid            /* Label bubble with hospital name     */
                             datalabelpos=topright         /* Label location at bubble            */
                             datalabelattrs=(color=black   /* Hospital name label info            */
                                             size=14
                                             style=normal
                                             weight=bold);
run;

proc iml;
use sao_carlos_pontos;
read all var{x y} into COORD;
n=nrow(coord[,1]);
d=j(1,3,0);                                                                                                                             
nome={"idi" "idj" "d"};                                                                                                                 
create dist2 from d[colname=nome];                                                                                                       
do i=1 to n;
	do j=i+1 to n;
	if abs(coord[,1])<180 then do;
		dif=abs(COORD[i,1]-COORD[j,1]);
		raio=arcos(-1)/180;                                                                                                              
		if dif=0 then arco=0;
		else
		/* Law of Cosines */  
		arco=arcos(sin(COORD[i,2]*raio)*sin(COORD[j,2]*raio)+cos(COORD[i,2]*raio)*cos(COORD[j,2]*raio)*cos(dif*raio));
		d[1]=i;  
		d[2]=j;
		d[3]=arco*6371 /*Earth's Radius = 6371 (aproximately)*/;
		append from d;
	end;
	else do;
		d[1]=i;  
		d[2]=j;
		d[3]=sqrt((COORD[i,1]-COORD[j,1])**2+(COORD[i,2]-COORD[j,2])**2);   
		append from d;
	end;
	end;
end;
close dist2;
quit;

data a;segment=2;v=2;run;
%macro m(dist);
%let dist2=%sysevalf(&dist/111);/*cada grau corresponde a aprox. 111km*/
data _null_;set Sao_carlos_pontos;
call symput('x'||trim(left(_n_)),x);
call symput('y'||trim(left(_n_)),y);
run;
%put &x1 &y1;
data circulo;
id=1;
do x=&x1-&dist2 to &x1+&dist2 by 0.001;
if x>=&x1+&dist2 then y=&y1;else y=sqrt(&dist2**2-(x-&x1)**2)+&y1;
output;
end;
run;
data circulo2;
id=1;
do x=&x1+&dist2 to &x1-&dist2 by -0.001;
if x>=&x1+&dist2 then y=&y1;else y=-sqrt(&dist2**2-(x-&x1)**2)+&y1;
output;
end;
run;
data circulo;set circulo circulo2;run;
data circulo;set circulo;x1=lag1(x);y1=lag1(y);if x1=. then delete;run;
proc sql noprint;select count(*) into:max from circulo;quit;%put &max;
data circulo_1;set circulo;if _n_=&max;x1=x;y1=y;run;
data circulo_2;set circulo;if _n_=1;x1=x;y1=y;run;
data circulo_2;merge circulo_2 circulo_1(drop=x1 y1);run;
data circulo;set circulo circulo_2;run;
data anno1;length function style color $8 position $1 text $20;retain xsys ysys '2' hsys '3' when 'a';
set circulo;
style='music';text='A';color='orange';size=0.2;position='5';line=1;
function='move';x=x;y=y;output;
function='draw';x=x1;y=y1;output;
run;
goptions ftitle='Verdana';
data sao_carlos_pontos;set sao_carlos_pontos;call symput('sec'||trim(left(FID)),trim(left(FID)));run;
%do i=1 %to 1;
data annoidi;set sao_carlos_pontos(where=(FID=&i));
retain hsys '3' xsys ysys "2";
function='label';text='J';color='orange';style='special';size=2;when='a';
run;
data annoidi2;set dist2(where=(idi=&i and d<=&dist));
retain hsys '3' xsys ysys "2";
function='label';text='J';color='black';style='special';size=3;when='a';
run;
data annoidi2;merge annoidi2(in=a rename=idj=fid) sao_carlos_pontos;by fid;if a;run;
data annoidi;set annoidi annoidi2;run;
title "Vizinhos do Ponto &&sec&i - com raio &dist";
proc gmap data=a map=Sao_carlos2 all anno=annoidi;
id segment;
choro v /discrete nolegend anno=anno1;
run;
quit;
%end;
%mend m;
%m(3);

/*********** usando distancia UTM **************/

proc iml;
use sao_carlos_pontos;
read all var{x y} into COORD;
n=nrow(coord[,1]);
d=j(1,3,0);                                                                                                                             
nome={"idi" "idj" "d"};                                                                                                                 
create dist2 from d[colname=nome];                                                                                                       
do i=1 to n;
	do j=i+1 to n;
		d[1]=i;  
		d[2]=j;
		d[3]=sqrt((COORD[i,1]-COORD[j,1])**2+(COORD[i,2]-COORD[j,2])**2);   
		append from d;
	end;
end;
close dist2;
quit;

data a;segment=2;v=2;run;
%macro m(dist);
%put &dist;
%let dist=%sysevalf(&dist);%put &dist;

data _null_;set sao_carlos_pontos;
call symput('x'||trim(left(_n_)),x);
call symput('y'||trim(left(_n_)),y);
run;
%put &x1 &y1;
data circulo;
id=1;
do x=&x1-&dist to &x1+&dist by 0.001;
if x>=&x1+&dist then y=&y1;else y=sqrt(&dist**2-(x-&x1)**2)+&y1;
output;
end;
run;
data circulo2;
id=1;
do x=&x1+&dist to &x1-&dist by -0.001;
if x>=&x1+&dist then y=&y1;else y=-sqrt(&dist**2-(x-&x1)**2)+&y1;
output;
end;
run;
data circulo;set circulo circulo2;run;
data circulo;set circulo;x1=lag1(x);y1=lag1(y);if x1=. then delete;run;
proc sql noprint;select count(*) into:max from circulo;quit;%put &max;
data circulo_1;set circulo;if _n_=&max;x1=x;y1=y;run;
data circulo_2;set circulo;if _n_=1;x1=x;y1=y;run;
data circulo_2;merge circulo_2 circulo_1(drop=x1 y1);run;
data circulo;set circulo circulo_2;run;
data anno1;length function style color $8 position $1 text $20;retain xsys ysys '2' hsys '3' when 'a';
set circulo;
style='music';text='A';color='orange';size=0.2;position='5';line=1;
function='move';x=x;y=y;output;
function='draw';x=x1;y=y1;output;
run;

goptions ftitle='Verdana';
data sao_carlos_pontos;set sao_carlos_pontos;call symput('sec'||trim(left(FID)),trim(left(FID)));run;
%do i=1 %to 1;
data annoidi;set sao_carlos_pontos(where=(FID=&i));
retain hsys '3' xsys ysys "2";
function='label';text='J';color='orange';style='special';size=2;when='a';
run;
data annoidi2;set dist2(where=(idi=&i and d<=&dist));
retain hsys '3' xsys ysys "2";
function='label';text='J';color='black';style='special';size=3;when='a';
run;
data annoidi2;merge annoidi2(in=a rename=idj=fid) sao_carlos_pontos;by fid;if a;run;
data annoidi;set annoidi annoidi2;run;
title "Vizinhos do Ponto &&sec&i - com raio &dist";
proc gmap data=a map=Sao_carlos2 all anno=annoidi;
id segment;
choro v /discrete nolegend anno=anno1;
run;
quit;
%end;
%mend m;
%m(5/111);

data dist2;set dist2;dist_km=d*111;run;
