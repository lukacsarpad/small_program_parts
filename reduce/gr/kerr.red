%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kerr spacetime                                                        %
% check metric                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% packages needed
load_package "taylor";
load_package "trigsimp";
load_package "compact";

in "../grfunc.red";
in "../trn.red";

% spatial dimensions
spdim := 3;


% coordinates
operator x;

x(0) := t;
x(1) := r;
x(2) := theta;
x(3) := phi;

% metric
%sigm := r^2 + a^2*cos(theta)^2 $
delt := r^2 - rg*r+a^2 $
sigm := r^2 + a^2 - a^2 * sin(theta)^2;
array gdd(3,3);
gdd(0,0) := (1-rg*r/sigm);
gdd(1,1) := - sigm/delt;
gdd(2,2) := - sigm;
gdd(3,3) := - ( r^2 + a^2 + rg*r*a^2/sigm * sin(theta)^2 ) * sin(theta)^2;
gdd(0,3) := gdd(3,0) :=  rg*r*a / sigm * sin(theta)^2;

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

for j:=0:3 do for k:=0:3 do eidd(j,k) := compact(trnr(eidd(j,k),theta),{sin(theta)^2+cos(theta)^2-1});

verifyvac();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Weyl tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comweyl();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Asymptotics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newtonian potential
npot := k * mass;

% asymptotic form of gdd(0,0)
sub(rp=r,sub(r=rp/epsilon,gdd(0,0))) $

taylor(ws,epsilon,0,2);

% note: gdd(0,0) ~ 1 + 2 phi, where phi is the Newtonian gravitational potential
% solve for Schwarzschild radius

solrg:=solve(coeffn(-r*taylortostandard(ws),epsilon,1)=2*npot,rg);

% postnewtonian rotation
pnrot := 2 * k * Angmom;

sub(rp=r,sub(r=rp/epsilon,gdd(0,3))) $
taylor(ws,epsilon,0,2);

% note: gdd(0,3) ~ 2 k J / r sin(theta)^2, where J is angular momentum
sola:=solve(coeffn(r/sin(theta)^2*taylortostandard(ws),epsilon,1)=pnrot,a);

sub(solrg,sola);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Invariants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comdualriemann();

cominvars();

Invar1 := trigsimp(Invar1) $ length(ws);

Invar1 := sub(abs(sin(theta)^2*a^2-a^2-r^2)=r^2+a^2-a^2*sin(theta)^2, Invar1) $ length(ws);

Invar1 := trnr(Invar1, theta) $ length(ws);

Invar1 := compact(Invar1, {sin(theta)^2+cos(theta)^2-1}) $ length(ws);

Invar1 := trigsimp(Invar1,cos);

Invar2 := trigsimp(Invar2) $ length(ws);

Invar2 := sub(abs(sin(theta)^2*a^2-a^2-r^2)=r^2+a^2-a^2*sin(theta)^2, Invar2) $ length(ws);

Invar2 := trnr(Invar2, theta) $ length(ws);

Invar2 := compact(Invar2, {sin(theta)^2+cos(theta)^2-1}) $ length(ws);

Invar2 := trigsimp(Invar2,cos);

% nonsingular at r=rg
sub(r=rg,Invar1);

sub(r=rg,Invar2);

% singularity at r=0:
limit(Invar1,r,0);

limit(Invar2,r,0);

% solve({(lambda1^2+lambda2^2+lambda1*lambda2)/3=invar1,lambda1*lambda2*(lambda1+lambda2)/2=invar2},{lambda1,lambda2});


end;