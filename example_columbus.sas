%let path=/home/alansilva0/my_shared_file_links/alansilva0/Estatistica_Espacial;
libname spatial "&path";
%include "&path/wmatrix_macro.sas";
%include "&path/moranscatterplot.sas";
%include "&path/lisa.sas";

data columbus;set spatial.columbus;run;
data crimeoh;set spatial.crimeoh;run;

/*
proc mapimport datafile='C:\spatial\COLUMBUS.SHP' out=columbus;run;
data columbus;set columbus(rename=code=code2);
if code2=1 then code=5;if code2=2 then code=1;if code2=3 then code=6;
if code2=4 then code=2;if code2=5 then code=7;if code2=6 then code=8;
if code2=7 then code=4;if code2=8 then code=3;if code2=9 then code=18;
if code2=10 then code=10;if code2=11 then code=38;if code2=12 then code=37;
if code2=13 then code=39;if code2=14 then code=40;if code2=15 then code=9;
if code2=16 then code=36;if code2=17 then code=11;if code2=18 then code=42;
if code2=19 then code=41;if code2=20 then code=17;if code2=21 then code=43;
if code2=22 then code=19;if code2=23 then code=12;if code2=24 then code=35;
if code2=25 then code=32;if code2=26 then code=20;if code2=27 then code=21;
if code2=28 then code=31;if code2=29 then code=33;if code2=30 then code=34;
if code2=31 then code=45;if code2=32 then code=13;if code2=33 then code=22;
if code2=34 then code=44;if code2=35 then code=23;if code2=36 then code=46;
if code2=37 then code=30;if code2=38 then code=24;if code2=39 then code=47;
if code2=40 then code=16;if code2=41 then code=14;if code2=42 then code=49;
if code2=43 then code=29;if code2=44 then code=25;if code2=45 then code=28;
if code2=46 then code=48;if code2=47 then code=15;if code2=48 then code=27;
if code2=49 then code=26;
drop code2;
run;*/ 


data a;
code=100;
v=2;
run;
%annomac;
proc sort data=columbus out=columbus2;by code;run;
%centroid(columbus2,centroid_columbus,code);
data anno;set centroid_columbus;
function='label';
style='simplex';
text=left(put(code,$2.));
size=0.85;
color='black';
xsys='2';
ysys='2';
when='a';
run;

title 'Shape of Columbus, Ohio';
proc gmap data=a map=columbus all;
id code;
choro v /nolegend anno=anno;
run;
quit;

title 'Shape of Columbus, Ohio';
proc gmap data=crimeoh map=columbus all;
id code;
choro crime;
run;
quit;

data anno2;set crimeoh(rename=(lat=y lon=x));
function='label';
style='simplex';
text=left(put(code,$2.));
size=0.65;
color='black';
xsys='2';
ysys='2';
when='a';
run;
proc ganno anno=anno2;run;

%WMATRIX(map=columbus, id=code, lat=y, long=x, type=contiguity, out=neighbor_columbus);

%let codigo=code;
%moran(Crimeoh,crime,Wmatrix_sdr,tipo=columbus);

proc spatialreg data=crimeoh Wmat=Wmatrix;
   model crime=income hvalue / type=SAR;
run;

proc spatialreg data=crimeoh Wmat=Wmatrix_sdr NONORMALIZE;
   model crime=income hvalue / type=SAR;
run;

proc spatialreg data=crimeoh Wmat=Wmatrix_compact;
   model crime=income hvalue / type=SAR;
   spatialid code;
run;

/*
data Wmatrix_compact2;set Wmatrix_compact;
if code=18 and ccode=32 then delete;
if code=32 and ccode=18 then delete;
if code=20 and ccode=33 then delete;
if code=33 and ccode=20 then delete;
if code=45 and ccode=47 then delete;
if code=47 and ccode=45 then delete;
run;
proc sql;
insert into Wmatrix_compact2
values (37,42,1) values (42,37,1);
quit;

proc spatialreg data=crimeoh Wmat=Wmatrix_compact2;
   model crime=income hvalue / type=SAR;
   spatialid code;
run;
proc spatialreg data=crimeoh Wmat=spatial2.matrizwpdr nonormalize;
   model crime=income hvalue / type=SAR;
run;
proc spatialreg data=spatial2.anselin Wmat=spatial2.matrizwpdr nonormalize;
   model crime=renda valor_casa / type=SAR;
run;
*/

data a1;
code=18;v='18';output;
code=32;v='32';output;
run;
title 'Regions 18 and 32 of Columbus, Ohio';
proc gmap data=a1 map=columbus;
id code;
choro v /nolegend anno=anno;
run;
quit;

data a2;
code=20;v='20';output;
code=33;v='33';output;
run;
title 'Regions 20 and 33 of Columbus, Ohio';
proc gmap data=a2 map=columbus;
id code;
choro v /nolegend anno=anno;
run;
quit;

data a3;
code=45;v='45';output;
code=47;v='47';output;
run;
title 'Regions 45 and 47 of Columbus, Ohio';
proc gmap data=a3 map=columbus;
id code;
choro v /nolegend anno=anno;
run;
quit;

data a4;
code=37;v='37';output;
code=42;v='42';output;
run;
title 'Regions 37 and 42 of Columbus, Ohio';
proc gmap data=a4 map=columbus;
id code;
choro v /nolegend anno=anno;
run;
quit;


data a5;
code=2;v='2';output;
code=3;v='3';output;
code=4;v='4';output;
run;
title 'Regions 2 and 4 of Columbus, Ohio';
proc gmap data=a5 map=columbus;
id code;
choro v /nolegend anno=anno;
run;
quit;


%WMATRIX(map=columbus, id=code, lat=y, long=x, type=distance, out=neighbor_columbus);

proc spatialreg data=crimeoh Wmat=Wmatrix;
   model crime=income hvalue / type=SAR;
run;

proc spatialreg data=crimeoh Wmat=Wmatrix_sdr NONORMALIZE;
   model crime=income hvalue / type=SAR;
run;

proc spatialreg data=crimeoh Wmat=Wmatrix_compact;
   model crime=income hvalue / type=SAR;
   spatialid code;
run;


%WMATRIX(map=columbus, id=code, lat=y, long=x, type=contiguity, neighbors=4, out=neighbor_columbus);

proc spatialreg data=crimeoh Wmat=Wmatrix;
   model crime=income hvalue / type=SAR;
run;

proc spatialreg data=crimeoh Wmat=Wmatrix_sdr NONORMALIZE;
   model crime=income hvalue / type=SAR;
run;

proc spatialreg data=crimeoh Wmat=Wmatrix_compact;
   model crime=income hvalue / type=SAR;
   spatialid code;
run;

%WMATRIX(map=columbus, id=code, lat=y, long=x, type=distance, distance=1, out=neighbor_columbus);

proc spatialreg data=crimeoh Wmat=Wmatrix;
   model crime=income hvalue / type=SAR;
run;

proc spatialreg data=crimeoh Wmat=Wmatrix_sdr NONORMALIZE;
   model crime=income hvalue / type=SAR;
run;

proc spatialreg data=crimeoh Wmat=Wmatrix_compact;
   model crime=income hvalue / type=SAR;
   spatialid code;
run;

/*** SIMULANDO DADOS ESPACIAIS ***/


proc iml;
use Wmatrix_sdr;
read all into wpdr;
n=nrow(wpdr);print n;
rho=0.3;
beta={-0.9,1.2,0.7,0.6,-1.0,0.8,-0.6};
x=j(n,7,0);
do i=1 to 7;
do j=1 to 49;
x[j,i]=rannor(2)*3+30*i;
end;
end;
z=1.8+x*beta+rannor(4);
print z x;
A=I(n)-rho*wpdr;
y=inv(A)*z;
print y x;
dados=y||x;
create datasim from dados[colname={"y" "x1" "x2" "x3" "x4" "x5" "x6" "x7"}];
append from dados;
quit;
proc spatialreg data=datasim Wmat=Wmatrix_sdr;
model y = x1 x2 x3 x4 x5 x6 x7 / type=SAR;
run;
