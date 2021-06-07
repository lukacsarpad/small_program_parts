%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Friedmann-Robertson-Walker spacetime                                  %
% check metric                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% packages needed
load_package "taylor";
load_package "trigsimp";
load_package "compact";

in "../grfunc.red";
in "../utils.red";
in "../trn.red";

% spatial dimensions
spdim := 3;


% coordinates
operator x;

x(0) := t;
x(1) := r;
x(2) := theta;
x(3) := phi;

% scale factor
depend a,t;

% curvature +,-, or zero depends on kk=0,+-1

% metric
array gdd(3,3), guu(3,3);
gdd(0,0) := 1;
gdd(1,1) := -a^2/(1-kk*r^2);
gdd(2,2) := -a^2*r^2;
gdd(3,3) := -a^2*r^2*sin(theta)^2;

operator dxu;
dxu(0) := dt;
dxu(1) := dr;
dxu(2) := dth;
dxu(3) := dph;

ds2 := for j:=0:3 sum for k:=0:3 sum gdd(j,k) * dxu(j) * dxu(k);

% calculate inverse metric
comguu();
rdetg:=sqrt(-detg);
rdetg:=sub({abs(-sigm)=sigm,abs(sin(theta))=sin(theta)},rdetg);


% volume form
filleta();
filleps3();
filleps4();

comvolform4();

%%%%%%%%%%%%%%%%%%%%%%%%% Christoffel symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%
comchristoffel();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Riemann tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%
comriemann();

%for j:=0:3 do for k:=0:3 do for l:=0:3 do for m:=0:3 do <<
%    Riedddd(j,k,l,m) := trigsimp(Riedddd(j,k,l,m));
%    Rieuddd(j,k,l,m) := trigsimp(Rieuddd(j,k,l,m));
%    Rieuudd(j,k,l,m) := trigsimp(Rieuudd(j,k,l,m));
%    Rieuuuu(j,k,l,m) := trigsimp(Rieuuuu(j,k,l,m));
%>>;

riemsymm();

bianchialg();

%%%%%%%%%%%%%%%%%%%%%%%% Verify Bianchi identity %%%%%%%%%%%%%%%%%%%%%%%%
bianchidiff();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ricci tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%
comricci();

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Einstein tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%
comeinstein();

%for j:=0:3 do for k:=0:3 do eidd(j,k) := compact(trnr(eidd(j,k),theta),{sin(theta)^2+cos(theta)^2-1});

%%%%%%%%%%%%%%%%%%%%%%%%%% Friedmann equations %%%%%%%%%%%%%%%%%%%%%%%%%%
% fluid velocity
array vd(3), vu(3);
vu(0) := 1/sqrt(gdd(0,0));

for j:=0:3 do vd(j) := for k:=0:3 sum gdd(j,k)*vu(k);

% fluid energy momentum tensor
array Tdd(3,3), Tud(3,3), Tuu(3,3);

for j:=0:3 do for k:=0:3 do Tdd(j,k) := (pressure + edens ) * vd(j) * vd(k) - pressure * gdd(j,k);
for j:=0:3 do for k:=0:3 do Tud(j,k) := (pressure + edens ) * vu(j) * vd(k) - pressure * if j=k then 1 else 0;
for j:=0:3 do for k:=0:3 do Tuu(j,k) := (pressure + edens ) * vu(j) * vu(k) - pressure * guu(j,k);

% 1st Friedmann equation: 0,0 compt. of Einstein equations
FrEq1 := Eiud(0,0) - Tud(0,0) - La $

begin scalar qq;
    qq := coeffn(freq1,edens,1);
    freq1 := -freq1/qq;
end;

freq1;

freq1p := 3*( df(a,t)^2 + kk ) / a^2 - ( edens + la ) $

freq1 - freq1p;

clear(freq1p);

% 2nd Friedmann equation: trace of Einstein equations
FrEq2 := ( for j:=0:3 sum  Eiud(j,j) - Tud(j,j) ) - 4 * La;

FrEq2 := df(a,t,2)/a - solv(FrEq2,df(a,t,2))/a;

freq2:=freq2-coeffn(freq2,df(a,t),2)/coeffn(freq1,df(a,t),2)*freq1;

freq2p := df(a,t,2) / a - ( -1/6*(edens+3*pressure) + La / 3) $

freq2 - freq2p;

clear(freq2p);

%%%%%%%%%%%%%%% Verify other Einstein equation components %%%%%%%%%%%%%%%
array eieqdd(3,3);

for i1:=0:3 do for i2:=0:3 do eieqdd(i1,i2):=Eiud(i1,i2)-Tud(i1,i2)-La * (if i1=i2 then 1 else 0);

% substitute Friedmann equations into the other components, shall be zero

for i1:=0:3 do for i2:=i1:3 do begin scalar qq;
    qq:=eieqdd(i1,i2)-coeffn(eieqdd(i1,i2),df(a,t),2)/coeffn(freq1,df(a,t),2)*freq1
	- coeffn(eieqdd(i1,i2),df(a,t,2),1)/coeffn(freq2,df(a,t,2),1)*freq2;

    if qq neq 0 then write "Error, Einstein equation component ", i1, " ", i2, " not expressed with Friedmann equations: ", qql
end;

write "Expressing Einstein equations with Friedmann equations verified.";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Weyl tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comweyl();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Invariants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comdualriemann();

cominvars();

invar1;

invar2;

invar1a := sub({df(a,t,2)=q2,df(a,t)=sqrt(q1)},invar1) $

invar2a := sub({df(a,t,2)=q2,df(a,t)=sqrt(q1)},invar2) $

on factor; invar2a; off factor;

end;
