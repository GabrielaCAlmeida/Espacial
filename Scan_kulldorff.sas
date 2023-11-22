libname scan '/home/alansilva0/my_shared_file_links/alansilva0/Estatistica_Espacial';
data veiculos;set scan.veiculos;run;
data setores_censitariosDF;set scan.setores_censitariosDF;run;
data anno_via_rodo;set scan.anno_via_rodo;run;

proc mapimport datafile='/home/alansilva0/my_shared_file_links/alansilva0/Estatistica_Espacial/53SEE250GC_SIR.shp' 
out=setores_df;
run;
data a;id=0;v=2;run;
proc gmap data=a map=setores_df(where=(CD_GEOCODS='53001080506')) all;
id id;
choro v /nolegend;
run;
quit;

data Veiculos2;
length function style color $8 position $1 text $20;
retain function 'label' xsys ysys '2' hsys '3' when 'a';
set Veiculos;
style='special';text='J';color='blue';size=3;position='5';
if caso=1 then color='red';
run;

data b;points_=0;v=2;run;
proc gmap data=b map=setores_censitariosDF anno=anno_via_rodo all;
id points_;
choro v /annotate=Veiculos2 cempty=gray nolegend name='br';
run;
quit;

%include '/home/alansilva0/my_shared_file_links/alansilva0/Estatistica_Espacial/lat_long_utm.sas';
%conversao(veiculos,0,23);
data Veiculos3;
length function style color $8 position $1 text $20;
retain function 'label' xsys ysys '2' hsys '3' when 'a';
set conversao;
style='special';text='J';color='blue';size=3;position='5';
if caso=1 then color='red';
rename lat_y=y long_x=x;
run;

data a;id=0;v=2;run;
proc gmap data=a map=setores_df(where=(CD_GEOCODS='53001080506')) all;
id id;
choro v /nolegend anno=Veiculos3;
run;
quit;

proc sgmap plotdata=veiculos3;
openstreetmap;
bubble x=x y=y size=size;
run;
quit;

data conversao;set conversao;size=1;id=_n_;run;
proc sgmap plotdata=conversao;
openstreetmap;
bubble x=long_x y=lat_y size=size;
run;
quit;


%macro scan_kulldorff(data=,id=,x=,y=,caso=,controle=,sim=999);
proc iml;
use &data;
read all var{&x &y &id} into coord;
read all var{&id &caso &controle} into _freq2_;
n=_freq2_[+,2];
pz=_freq2_[,3]/_freq2_[+,3];print _freq2_;
uz=n*pz;
_freq2_=_freq2_||pz||uz;
print _freq2_[label="Dados"];
_freq_=(1:nrow(_freq2_))`||_freq2_[,2:ncol(_freq2_)];

logLz=j(1,4,0);
create kulldorff from logLz[colname={"i" "j" "&id" "lambda"}];
do i=1 to nrow(_freq_);
 dist=j(nrow(_freq_),3,0);
 do j=1 to nrow(_freq_);
  dist[,1]=i;
  dist[j,2]=j;
  dist[j,3]=sqrt((COORD[i,1]-COORD[j,1])**2+(COORD[i,2]-COORD[j,2])**2);
 end;
 call sort(dist,{1,3});
pos=dist[loc(dist[,1]=_freq_[i,1]),2]`;
*print pos;
do j=1 to round(nrow(_freq_)/2);
logLz[1]=i;
logLz[2]=j;
logLz[3]=_freq2_[pos[j],1];
if _freq_[pos[1:j],2][+]>_freq_[pos[1:j],5][+] then
logLz[4]=_freq_[pos[1:j],2][+]#log(_freq_[pos[1:j],2][+]/_freq_[pos[1:j],5][+])+
(_freq_[+,2]-_freq_[pos[1:j],2][+])#log((_freq_[+,2]-_freq_[pos[1:j],2][+])/
(_freq_[+,2]-_freq_[pos[1:j],5][+]));
else logLz[4]=0;
*print i j (pos[1:j]) logLz;
/*if i=1 & j=1 then lambda=logLz;
else do;
if logLz >0 then lambda=lambda//logLz;*/
append from logLz;
*end;
end;
end;
close kulldorff;

use kulldorff;
read all into lambda;

call sort(lambda,{1 2});
loc1=lambda[loc(lambda[,4]=lambda[<>,4]),1:4];
loc2=lambda[loc(lambda[,1]=loc1[1]),1:4];
cluster=loc2[loc(loc2[,2]<=loc1[2]),];
print (cluster[,3])[label="Cluster"] (loc1[4])[label="Lambda"];
create _cluster_ from cluster[colname={"i" "j" "&id" "lambda"}];
append from cluster;
close _cluster_;

/**** gerando simula��es **********/

call randseed(1); 
prob = _freq_[,3]/_freq_[+,3]; 
NumTrials = _freq_[+,2]; 
sim = RANDMULTINOMIAL(&sim,NumTrials,prob);
print sim;

do k=1 to nrow(sim);
_freq_[,2]=sim[k,]`;
do i=1 to nrow(_freq_);
dist=j(nrow(_freq_),3,0);
 do j=1 to nrow(_freq_);
  dist[,1]=i;
  dist[j,2]=j;
  dist[j,3]=sqrt((COORD[i,1]-COORD[j,1])**2+(COORD[i,2]-COORD[j,2])**2);
 end;
 call sort(dist,{1,3});
pos=dist[loc(dist[,1]=_freq_[i,1]),2]`;
*print pos;
do j=1 to round(nrow(_freq_)/2);
if _freq_[pos[1:j],2][+]>_freq_[pos[1:j],5][+] then
logLz=_freq_[pos[1:j],2][+]#log(_freq_[pos[1:j],2][+]/_freq_[pos[1:j],5][+])+
(_freq_[+,2]-_freq_[pos[1:j],2][+])#log((_freq_[+,2]-_freq_[pos[1:j],2][+])/
(_freq_[+,2]-_freq_[pos[1:j],5][+]));
else logLz=0;
if i=1 & j=1 then lambdas=logLz;
else lambdas=lambdas//logLz;
end;
end;
call sort(lambdas,{1},{1});
lambdasm=lambdas[1];
if k=1 then lambda_s=lambdasm;
else lambda_s=lambda_s//lambdasm;
end;
create dist from dist;append from dist;
lambda_s=lambda_s//lambda[1,4];
create lambda_s from lambda_s;
append from lambda_s;
call sort(lambda_s,{1});
critico=lambda_s[round(0.95*nrow(lambda_s))];
p_valor=ncol(loc(lambda_s>=loc1[4]))/(&sim+1);
print critico p_valor;
quit;
%mend scan_kulldorff;

%scan_kulldorff(data=veiculos,id=id,x=x,y=y,caso=caso,controle=pop,sim=29);


data Veiculos2;length function style color $8 position $1 text $20;retain function 'label' xsys ysys '2' hsys '3' when 'a';
set Veiculos;
style='special';text='J';color='black';size=1;position='5';
run;

data _cluster_;set Kulldorff;
where i=137 and j<=6;
run;
proc sort data=_cluster_;by id;run;
data _cluster_;merge veiculos2 _cluster_;by id;
if lambda ne . then color='red';
run;
data b;points_=0;v=2;run;
proc gmap data=b map=setores_censitariosDF anno=anno_via_rodo all;
id points_;
choro v /annotate=_cluster_ cempty=gray nolegend name='br2';
run;
quit;
proc gchart data=lambda_s;
vbar3d col1 /space=0;
run;
quit;


/************ acidentes area *********/

data descricaoDF4;set scan.descricaoDF4;run;
data loc1_dentroDF;set scan.loc1_dentroDF;run;
data End1_dentrodf;set scan.End1_dentrodf;run;
data anno;set scan.anno;run;
data centroides;set scan.centroides;run;

proc gmap data=descricaoDF4 map=descricaoDF4 anno=anno;
id setor1;
choro setor1 /cempty=gray annotate=loc1_dentroDF discrete;
footnote font='special' color='black' "  J"  font='Verdana' color='black' '   Acidentes';
run;
quit;

proc freq data=Loc1_dentrodf noprint;
tables setor1 / out=casos(drop=percent rename=count=caso);
run;
proc freq data=End1_dentrodf noprint;
tables setor1 / out=controles(drop=percent rename=count=controle);
weight veiculos;
run;
data dadosBSB;merge centroides(rename=id=setor1) casos controles;by setor1;run;

%scan_kulldorff(data=dadosbsb,id=setor1,x=x,y=y,caso=caso,controle=controle);

proc gchart data=lambda_s;
vbar3d col1 /space=0;
run;
quit;
proc means data=lambda_s mean;run;
data poisson;
do i=1 to 1000;
col1=ranpoi(2,1.99);
output;
end;
run;
proc gchart data=poisson;
vbar3d col1 /discrete space=0;
run;
quit;

proc gmap data=_cluster_ map=descricaoDF4 anno=anno all;
id setor1;
choro lambda /cempty=gray annotate=loc1_dentroDF discrete;
footnote font='special' color='black' "  J"  font='Verdana' color='black' '   Acidentes';
run;
quit;

proc sql;
create table _cluster_ as
select a.*,x,y from _cluster_ a,dadosbsb b
where a.setor1=b.setor1;
quit;
%conversao(_cluster_,0,23);
data _cluster_;merge _cluster_ conversao;run;

proc sgmap plotdata=_cluster_;
openstreetmap;
bubble x=long_x y=lat_y size=lambda;
run;
quit;

proc sgmap plotdata=_cluster_;
esrimap url='https://services.arcgisonline.com/arcgis/rest/services/World_Imagery/MapServer';
bubble x=long_x y=lat_y size=lambda;
run;
quit;
