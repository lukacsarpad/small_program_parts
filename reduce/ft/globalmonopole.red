%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global monopole radial equation                                      %
%                                                                      %
% References:                                                          %
%   * Polyakov, JETP Lett. 20, 194-195, (1974)                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

in "../grfunc.red";
in "../ftfunc.red";

% 3 d
spdim := 3;

% flat metric, cartesian coordinates
operator x;
x(0) := t;
x(1) := r;
x(2) := theta;
x(3) := phi;

setgpolar();
comguu3();

rdetg:=sqrt(-detg);
rdetg:=sub({ abs(sin(theta)) = sin(theta), abs(r) = r }, rdetg);

% Euclidean coordinates for Hedgehog Ansatz
operator xc;
xc(1) := r*sin(theta)*cos(phi);
xc(2) := r*sin(theta)*sin(phi);
xc(3) := r*cos(theta);

% Hedgehog Ansatz
array phi1(3);
depend u,r;
for a:=1:3 do phi1(a) := xc(a)*u/r;

% Lagrangian
lag := 1/2*( for j:=0:spdim sum for k:=0:spdim sum for a:=1:3 sum df( phi1(a), x(j) ) * df( phi1(a), x(k) ) * guu(j,k) )
	+ mu^2/2 * ( for a:=1:3 sum phi1(a)^2 ) - lambda1 / 4 * ( for a:=1:3 sum phi1(a)^2 )^2 $ length(ws);

lag := trigsimp(lag);

% derive equations of motion
begin scalar lagp, dldpiu;
    lagp := rdetg * sub( df(u, r) = piu, lag);
    dldpiu := df(lagp, piu);
    dldpiu := sub( piu = df(u, r), dldpiu);
    eomu := (1/rdetg) * ( df(dldpiu, r) - df(lagp, u) );

    eomu := eomu / coeffn(eomu, df(u,r,2),1);
end;

% verify
eomup := (1/r^2)* df( r^2*df(u,r), r) + (mu^2 -2/r^2)*u - lambda1*u^3;

eomu - eomup;

clear eomup;

% use full EoM, formula verified in scalareom.red
comchristoffel();

array dphid(3,spdim);
for a:=1:3 do for j:=0:spdim do dphid(a,j) := df(phi1(a), x(j));

array dphiu(3,spdim);
for a:=1:3 do for j:=0:spdim do dphiu(a,j) := for k:=0:spdim sum guu(j,k)*dphid(a,k);

array ddphi(3);
for a:=1:3 do ddphi(a) := for j:=0:spdim sum df( dphiu(a,j), x(j)) + for l:=0:spdim sum dphiu(a,j)*gaudd(l,j,l);

array eomphi(3);
for a:=1:3 do eomphi(a) := ddphi(a) - mu^2*phi1(a) + lambda1 * ( for b:=1:3 sum phi1(b)*phi1(b) ) * phi1(a);

for a:=1:3 do eomphi(a) := trigsimp(eomphi(a));

for a:=1:3 do begin scalar qq;
    qq := eomphi(a) - coeffn(eomphi(a), df(u,r,2),1) / coeffn(eomu, df(u,r,2),1) * eomu;
    if qq neq 0 then write "*** Error in eomphi(",a,") = ", qq;
end;

write "Consistency of the Ansatz verified: calculated all components of the field equations.";

end;

% alternate ending: read output of scalareom.red

off echo;
in "scalareom2.out";
%in "scalareom.out";
on echo;

for a:=1:3 do eomphi(a) := trigsimp(eomphi(a));

for a:=1:3 do begin scalar qq;
    qq := eomphi(a) - coeffn(eomphi(a), df(u,r,2),1) / coeffn(eomu, df(u,r,2),1) * eomu;
    if qq neq 0 then write "*** Error in eomphi(",a,") = ", qq;
end;

end;
