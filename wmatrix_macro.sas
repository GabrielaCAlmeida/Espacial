%macro Wmatrix(map=,id=,lat=,long=,type=contiguity,centroid_lat=,centroid_long=,
distance=,neighbors=,out=neighbor);

%if %upcase(&type)=CONTIGUITY %then %do;
%if &neighbors = %then %do;
/********Defining the Neighborhood for contiguity ***********/
proc sort data=&map out=&map.2 nodupkey; by &long &lat &id;run;
data _hashmap_ (keep=&long &lat z &id);
retain z;
set &map.2(where=(&long NE . and &lat NE .));
by &long &lat;
if first.&long or first.&lat then z=1;
z=z+1;
run;
proc means data=_hashmap_ noprint;
output out=maxz max(z)=mz; run;
/*put the maximum value of z into macro variable maxnb*/
data _null_;
set maxz;
call symputx("maxnb",mz);
run;
%put &maxnb;
proc sort data=&map.2; by &id;run;
data nonb (keep=myid);
/*length &long &lat z 8; format &id 16.;*/
if _n_=1 then do;
/*this hash object holds a second copy of the entire map for comparison*/
declare hash fc(dataset: "_hashmap_");
fc.definekey("&long","&lat","z");
fc.definedata("&long","&lat","z","&id");
fc.definedone();
call missing (&long,&lat,z,&id);
/*this hash object will hold the rook neighbors for each area: they have two
points in common*/
declare hash rook();
rook.definekey("&id","myid");
rook.definedata("&id","myid");
rook.definedone();
end;
/*this hash object holds the bishop neighbors for each area: they have a point
in common*/
declare hash bishop();
bishop.definekey("&id","myid");
bishop.definedata("&id","myid");
bishop.definedone();
foundnb="N";
do until (last.myid);
set &map.2 (keep=&id &long &lat rename=(&id=myid &long=myx &lat=myy) where=(myx NE .
and myy NE .)) end=eof;
by myid;
do n=1 to &maxnb.; /*this is max number of points in common =max z*/
rc=fc.find(key:myx, key:myy, key:n);
if rc=0 and myid NE &id then do;
nbrc=rook.check(key:&id, key:myid);
if nbrc=0 then do;
rc2=rook.add(key:&id, key:myid, data:&id, data:myid);
foundnb="Y";
end;
else rc1=rook.add(key:&id, key:myid, data:&id, data:myid);
end;
end;*do &maxnb.;
end;*end DOW loop;
if foundnb="N" then output nonb;
if eof then rook.output(dataset:"&out");
run;
proc sort data=&out;by &id myid;run;

proc sort data=&out(keep=&id) nodupkey out=_tab_;by &id;run;
data _tab_;set _tab_;idi=_n_;run;
proc sort data=&out;by &id;run;
data _tab2_;merge &out _tab_;by &id;run;
proc sort data=_tab2_;by myid;run;
data _tab2_;merge _tab2_ _tab_(rename=(idi=idj &id=myid));by myid;run;
proc sort data=_tab2_;by &id;run;
proc sql;drop table _tab_,_hashmap_,Nonb,Maxz,&map.2;quit;

data Wmatrix_compact;set &out(rename=myid=c&id);
value=1;
run;
%end;
%else %do;
%if &centroid_lat= or &centroid_long= %then %do;
%annomac;
proc sort data=&map(keep=&id &lat &long rename=(&long=y &lat=x)) out=_&map;by &id;run;
%centroid(_&map,&out,&id);
proc sql;drop table _&map;quit;
%end;
%else %do;
data &out;set &map(keep=&id &centroid_lat &centroid_long rename=(&centroid_long=y &centroid_lat=x));run;
%end;

proc iml;
use &out;
read all var{x y} into COORD;
n=nrow(coord[,1]);
*print n;
d2 = distance(COORD, "L2");
n2=n**2;
d=shape(d2,n2,1);
create _dist_ from d;
append from d;
/*d=j(1,3,0);                                                                                                                             
nome={"idi" "idj" "d"};                                                                                                                 
create _dist_ from d[colname=nome];                                                                                                       
do i=1 to n;
	do j=1 to n;
		d[1]=i;  
		d[2]=j;
		d[3]=sqrt((COORD[i,1]-COORD[j,1])**2+(COORD[i,2]-COORD[j,2])**2);   
		append from d;
	end;
end;*/
close _dist_;
quit;
proc sql noprint;select count(*) into:nrw from &out;quit;
data _null;do idi=1 to &nrw;do idj=1 to &nrw;output;end;end;run;
data _dist_;merge _null _dist_;rename col1=d;if idi>=idj then delete;run;
data _dist_;set _dist_;if d=0 then delete;run;
proc sort data=_dist_;by idi d;run;
data _dist_;retain seq 0;set _dist_;by idi;
if first.idi then seq=1;
else seq+1;
run;
data _tab2_;set _dist_(where=(seq<=&neighbors));run;
data &out;set _tab2_;run;
data Wmatrix_compact;set &out;
rename idi=&id idj=c&id;value=1;
drop seq d;
run;
proc sql;drop table _dist_,_null;quit;
%end;

/********creating W Matrix ***********/

proc iml;
use _tab2_ var{idi idj};
read all;
n=max(idj);
w=j(n,n,0);
do h=1 to nrow(idi);                                                                                                                   
w[idi[h],idj[h]]=1;
end;
wpdr=j(n,n,0);                                                                                                                        
soma=j(n,1,0);                                                                                                                         
do i=1 to n;                                                                                                                           
do j=1 to n;                                                                                                                           
soma[i]=soma[i]+w[i,j];                                                                                                                 
end;                                                                                                                                    
end;                                                                                                                                    
do i=1 to n;                                                                                                                           
do j=1 to n;      
if soma[i]=0 then wpdr[i,j]=0;
else wpdr[i,j]=w[i,j]/soma[i];  
end;                                                                                                                                    
end; 
create Wmatrix from w;                                                                                                                  
append from w;                                                                                                                          
create Wmatrix_sdr from wpdr;                                                                                                            
append from wpdr;                                                                                                                       
quit;
proc sql;drop table _tab2_;quit;
%end;
%else %if %upcase(&type)=DISTANCE %then %do;

%if &centroid_lat= or &centroid_long= %then %do;
%annomac;
proc sort data=&map(keep=&id &lat &long rename=(&long=y &lat=x)) out=_&map;by &id;run;
%centroid(_&map,&out,&id);
proc sql;drop table _&map;quit;
%end;
%else %do;
data &out;set &map(keep=&id &centroid_lat &centroid_long rename=(&centroid_long=y &centroid_lat=x));run;
%end;

proc iml;
use &out;
read all var{x y} into COORD;
n=nrow(coord[,1]);
*print n;
d2 = distance(COORD, "L2");
n2=n**2;
d=shape(d2,n2,1);
create _dist_ from d;
append from d;
/*d=j(1,3,0); 
nome={"idi" "idj" "d"}; 
create _dist_ from d[colname=nome];                                                                                                       
do i=1 to n;
	do j=i+1 to n;
		d[1]=i;  
		d[2]=j;
		d[3]=sqrt((COORD[i,1]-COORD[j,1])**2+(COORD[i,2]-COORD[j,2])**2);   
		append from d;
	end;
end;*/
close _dist_;
quit;
proc sql noprint;select count(*) into:nrw from &out;quit;
data _null;do idi=1 to &nrw;do idj=1 to &nrw;output;end;end;run;
data _dist_;merge _null _dist_;rename col1=d;if idi>=idj then delete;run;
data &out;set _dist_;run;
proc sql;drop table _dist_,_null;quit;

data Wmatrix_compact;
set &out(rename=(idi=&id idj=c&id)) &out(rename=(idi=c&id idj=&id));
%if &distance= %then %do; 
value=1/d;
%end;
%else %do;
if d<=&distance then value=1;
if value=. then delete;
%end;
keep &id c&id value;
run;
proc sort data=Wmatrix_compact;by &id c&id;run;

/********creating W Matrix ***********/
proc iml;
use &out var{idi idj d};
read all;
n=max(idj);
w=j(n,n,0);
do h=1 to nrow(idi);
%if &distance= %then %do; 
w[idi[h],idj[h]]=1/d[h];
w[idj[h],idi[h]]=1/d[h];
%end;
%else %do;
if d[h]<=&distance then do;
w[idi[h],idj[h]]=1;
w[idj[h],idi[h]]=1;
end;
%end;
end;

wpdr=j(n,n,0);                                                                                                                        
soma=j(n,1,0);                                                                                                                         
do i=1 to n;                                                                                                                           
do j=1 to n;                                                                                                                           
soma[i]=soma[i]+w[i,j];                                                                                                                 
end;                                                                                                                                    
end;                                                                                                                                    
do i=1 to n;                                                                                                                           
do j=1 to n;      
if soma[i]=0 then wpdr[i,j]=0;
else wpdr[i,j]=w[i,j]/soma[i];  
end;                                                                                                                                    
end; 
create Wmatrix from w;                                                                                                                  
append from w;                                                                                                                          
create Wmatrix_sdr from wpdr;                                                                                                            
append from wpdr;                                                                                                                       
quit;
%end;

%mend Wmatrix;
