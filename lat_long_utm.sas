%macro llUTM(tab);
data &tab;set &tab;
/*y=-21.75223512;
x=-52.18464271;*/
_a_=6378137;
_1f_=298.2572235630;
e2=2*(1/_1f_)-(1/_1f_)**2;
e22=+(e2)/(1-e2);
k0=0.9996;

ry=(3.14159/180)*y;
rx=(3.14159/180)*x;
LongTemp=+(x+180)-int((x+180)/360)*360-180;
zone=+int((LongTemp+180)/6)+1;
LongOrigin=+(zone-1)*6-180+3;
rLongOrigin=(3.14159/180)*LongOrigin;
N=+_a_/(sqrt(1-e2*sin(ry)**2));
T=tan(ry)**2;
C=e22*cos(ry)**2;
A=+(rx-rLongOrigin)*cos(ry);
M=+_a_*((1-e22/4-3*e2**2/64-5*e2**3/256)*ry-(3*e2/8+3*e2**2/32+45*e2**3/1024)*sin(2*ry)+
(15*e2**2/256+45*e2**3/1024)*sin(4*ry)-(35*e2**3/3072)*sin(6*ry));
UTMy=+k0*(M+N*tan(ry)*(A**2/2+(5-T+9*C+4*C**2)*A**4)/24+((61-58*T+T**2+600*C-330*e22)*A**6)/720)+10000000;
UTMx=+K0*N*(A+((1-T+C)*A**3)/6+((5-18*T+T**2+72*C-58*e22)*A**5)/120)+500000;
drop _a_ _1f_ e2 e22 k0 ry rx LongTemp LongOrigin rLongOrigin N T C A M;
run;
%mend llUTM;

%macro UTMll(tab);
data &tab;set &tab;
/*y=7594597.009898804;
x=377432.8111434143;*/
x=x-500000;
y=y-10000000;
k0=0.9996;
_a_=6378137;
_1f_=298.2572235630;
e2=2*(1/_1f_)-(1/_1f_)**2;
e22=+(e2)/(1-e2);

M = y/k0;
mu = M/(_a_*(1 - e2/4 - 3*e2**2/64 - 5*e2**3/256 - 7*e2**4/1024));
e1 = (1 - sqrt(1 - e2))/(1 + sqrt(1 - e2));
J1 = (3*e1)/2 - (27*e1**3)/32;
J2 = (21*e1**2)/16 - (55*e1**4)/32;
J3 = (151*e1**3)/96;
J4 = (1097*e1**4)/512; 
fp = mu + J1*sin(2*mu) + J2*sin(4*mu) + J3*sin(6*mu) + J4*sin(8*mu);

C1 = e22*cos(fp)**2;
T1 = tan(fp)**2;
R1 = (_a_*(1-e2))/(1-e2*sin(fp)**2)**(3/2);
N1 = _a_/sqrt(1-e2*sin(fp)**2);
D = x/(N1*k0);

Q1 = N1*tan(fp)/R1;
Q2 = (D**2/2);
Q3 = (5 + 3*T1 + 10*C1 - 4*C1**2 -9*e22)*D**4/24;
Q4 = (61 + 90*T1 + 298*C1 +45*T1**2  - 3*C1**2 -252*e22)*D**6/720;
lat = fp - Q1*(Q2 - Q3 + Q4);

Q5 = D;
Q6 = (1 + 2*T1 + C1)*D**3/6;
Q7 = (5 - 2*C1 + 28*T1 - 3*C1**2 + 8*e22 + 24*T1**2)*D**5/120;
LongTemp=+(x+180)-int((x+180)/360)*360-180;
zone=abs(int((LongTemp+180)/6)+1);
long0=+(zone-1)*6-180+3;
long = long0 + (Q5 - Q6 + Q7)/cos(fp);
put lat long;
run;
%mend UTMll;


%macro conversao(tab,latlontoUTM,zona);
proc iml;
pi = arcos(-1);
sm_a = 6378160.0; 
sm_b = 6356752.314;
UTMScaleFactor = 0.9996;

start DegToRad(deg) global(pi);
deg=(deg / 180.0) * pi;
return (deg);
finish DegToRad;

start RadToDeg(rad) global(pi);
rad=(rad / pi) * 180.0;
return (rad);
finish RadToDeg;

start ArcLengthOfMeridian (phi) global(sm_a,sm_b);
n = (sm_a - sm_b) / (sm_a + sm_b);
alpha = ((sm_a + sm_b) / 2.0)* ((1.0 + n**2/ 4.0) + (n**4/ 64.0));
beta = (-3.0 * n / 2.0) + (9.0 * n**3 / 16.0)+ (-3.0 *n**5/ 32.0);
gamma = (15.0 * n**2 / 16.0)+ (-15.0 * n**4 / 32.0);
delta = (-35.0 * n**3 / 48.0)+ ((105.0 * n**5) / 256.0);
epsilon = (315.0 * n**4.0/ 512.0);
result = alpha* (phi + (beta * sin(2.0 * phi))
        + (gamma * sin(4.0 * phi))
        + (delta * sin(6.0 * phi))
        + (epsilon * sin(8.0 * phi)));
return (result);
finish ArcLengthOfMeridian;

start UTMCentralMeridian (zone);
cmeridian = DegToRad (-183.0 + (zone * 6.0));
return (cmeridian);
finish UTMCentralMeridian;

start FootpointLatitude (y) global(sm_a,sm_b);
n = (sm_a - sm_b) / (sm_a + sm_b);
alpha_ = ((sm_a + sm_b) / 2)* (1 + (n**2/ 4) + (n**4/64));
y_ = y / alpha_;
beta_ = (3 * n / 2) + (-27 * n**3 / 32)+ (269 * n**5/ 512);
gamma_ = (21 * n**2/ 16)+ (-55 * n**4/ 32);
delta_ = (151 * n**3/ 96)+ (-417 * n**5/ 128);
epsilon_ = (1097 * n**4/ 512);
result = y_ + (beta_ * sin (2 * y_))
            + (gamma_ * sin (4 * y_))
            + (delta_ * sin (6 * y_))
            + (epsilon_ * sin (8 * y_));
return (result);
finish FootpointLatitude;

start MapLatLonToXY (phi, lambda, lambda0) global(sm_a,sm_b,xy);
ep2 = (sm_a**2 - sm_b**2) / sm_b**2;
nu2 = ep2 * cos(phi)**2;
N = sm_a**2/(sm_b *sqrt (1 + nu2));
t = tan(phi);
t2 = t * t;
tmp = (t2 * t2 * t2) - t**6;
l = lambda - lambda0;
l3coef = 1.0 - t2 + nu2;
l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2 * nu2);
l5coef = 5.0 - 18.0 * t2 + (t2 * t2) + 14.0 * nu2 - 58.0 * t2 * nu2;
l6coef = 61.0 - 58.0 * t2 + (t2 * t2) + 270.0 * nu2 - 330.0 * t2 * nu2;
l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2 * t2) - (t2 * t2 * t2);
l8coef = 1385.0 - 3111.0 * t2 + 543.0 * (t2 * t2) - (t2 * t2 * t2);
*xy=j(1,2,0);
xy[1] = N * cos (phi) * l
       + (N / 6.0 * cos(phi)**3 * l3coef * l**3)
       + (N / 120.0 * cos(phi)**5.0* l5coef * l**5.0)
       + (N / 5040.0 * cos(phi)**7.0 * l7coef * l**7.0);
xy[2] = ArcLengthOfMeridian (phi)
       + (t / 2.0 * N * cos(phi)**2 * l**2)
       + (t / 24.0 * N * cos(phi)**4 * l4coef * l**4)
       + (t / 720.0 * N * cos(phi)**6 * l6coef * l**6)
       + (t / 40320.0 * N * cos(phi)**8 * l8coef * l**8);
return (xy);
finish MapLatLonToXY;

start MapXYToLatLon (x, y, lambda0) global(sm_a,sm_b,philambda);
phif = FootpointLatitude (y);
ep2 = (sm_a**2 - sm_b**2)/sm_b**2;
cf = cos(phif);
nuf2 = ep2 * cf**2;
Nf = sm_a**2/(sm_b * sqrt(1 + nuf2));
Nfpow = Nf;
tf = tan(phif);
tf2 = tf * tf;
tf4 = tf2 * tf2;
x1frac = 1.0 / (Nfpow * cf);
Nfpow = Nfpow*Nf;   /* now equals Nf**2) */
x2frac = tf / (2.0 * Nfpow);
Nfpow = Nfpow*Nf;   /* now equals Nf**3) */
x3frac = 1.0 / (6.0 * Nfpow * cf);
Nfpow = Nfpow*Nf;   /* now equals Nf**4) */
x4frac = tf / (24.0 * Nfpow);
Nfpow = Nfpow*Nf;   /* now equals Nf**5) */
x5frac = 1.0 / (120.0 * Nfpow * cf);
Nfpow = Nfpow*Nf;   /* now equals Nf**6) */
x6frac = tf / (720.0 * Nfpow);
Nfpow = Nfpow*Nf;   /* now equals Nf**7) */
x7frac = 1.0 / (5040.0 * Nfpow * cf);
Nfpow = Nfpow*Nf;   /* now equals Nf**8) */
x8frac = tf / (40320.0 * Nfpow);
x2poly = -1.0 - nuf2;
x3poly = -1.0 - 2 * tf2 - nuf2;
x4poly = 5.0 + 3.0 * tf2 + 6.0 * nuf2 - 6.0 * tf2 * nuf2
	- 3.0 * (nuf2 *nuf2) - 9.0 * tf2 * (nuf2 * nuf2);
x5poly = 5.0 + 28.0 * tf2 + 24.0 * tf4 + 6.0 * nuf2 + 8.0 * tf2 * nuf2;
x6poly = -61.0 - 90.0 * tf2 - 45.0 * tf4 - 107.0 * nuf2
	+ 162.0 * tf2 * nuf2;
x7poly = -61.0 - 662.0 * tf2 - 1320.0 * tf4 - 720.0 * (tf4 * tf2);
x8poly = 1385.0 + 3633.0 * tf2 + 4095.0 * tf4 + 1575 * (tf4 * tf2);
*philambda=j(1,2,0);
philambda[1] = phif + x2frac * x2poly * (x * x)
	+ x4frac * x4poly * x**4
	+ x6frac * x6poly * x**6
	+ x8frac * x8poly * x**8;
philambda[2] = lambda0 + x1frac * x
	+ x3frac * x3poly * x**3
	+ x5frac * x5poly * x**5
	+ x7frac * x7poly * x**7;
return (philambda);
finish MapXYToLatLon;

start LatLonToUTMXY (lat, lon, zone) global(UTMScaleFactor);
h=UTMCentralMeridian (zone);
xy=MapLatLonToXY (lat, lon, h);
/* Adjust easting and northing for UTM system. */
xy[1] = xy[1] * UTMScaleFactor + 500000.0;
xy[2] = xy[2] * UTMScaleFactor;
if (xy[2] < 0) then xy[2] = xy[2] + 10000000.0;
return (xy);
finish LatLonToUTMXY;

start UTMXYToLatLon (x, y, zone, southhemi) global(UTMScaleFactor);
x = x-500000;
x = x/UTMScaleFactor;
/* If in southern hemisphere, adjust y accordingly. */
if (southhemi=1) then do;
y = y-10000000;
y = y/UTMScaleFactor;
end;
cmeridian = UTMCentralMeridian (zone);
xy=MapXYToLatLon (x, y, cmeridian);
return (xy);
finish UTMXYToLatLon;

use &tab var{x y};
read all;
if &latlontoUTM=1 then do;
xy=j(1,2,0);
create conversao from xy[colname={"UTMx" "UTMy"}];   
do i=1 to nrow(y);
lat=y[i];
lon=x[i];
zone = int((lon + 180.0) / 6 + 1);*print zone;
utm=LatLonToUTMXY (DegToRad (lat), DegToRad (lon), zone);
xy=utm;
append from xy;
end;
end;
else do;
philambda=j(1,2,0);
create conversao from philambda[colname={"Lat_y" "Long_x"}];   
do i=1 to nrow(y);
utmx = x[i] - 75;
utmy = y[i] - 25;
latlong=UTMXYToLatLon(utmx,utmy,&zona,1);
latlong[1]=RadToDeg(latlong[1]);
latlong[2]=RadToDeg(latlong[2]);
philambda=latlong;
append from philambda;
end;
end;

/*x=-52.18532908148965;
y=-21.748195847064192;
lat=y;
lon=x;
xy=j(nrow(y),2,0);
zone = int((lon + 180.0) / 6 + 1);print zone;
utm=LatLonToUTMXY (DegToRad (lat), DegToRad (lon), zone);
print utm[colname={"UTMx" "UTMy"}];

y=7594597.009898804;
x=377432.8111434143;
x = x - 75;
y = y - 25;
philambda=j(nrow(y),2,0);
latlong=UTMXYToLatLon(x,y,22,1);
latlong[1]=RadToDeg(latlong[1]);
latlong[2]=RadToDeg(latlong[2]);
print latlong;*/
quit;
%mend conversao;
