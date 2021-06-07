%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic monopole radial equations                                   %
%                                                                      %
% References:                                                          %
%   * Polyakov, JETP Lett. 20, 194-195, (1974)                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

in "../grfunc.red";
in "../ftfunc.red";

% 3 d
spdim := 3;

% flat metric, polar coordinates
operator x;
x(0) := t;
x(1) := r;
x(2) := theta;
x(3) := phi;

setgpolar();

comguu3();

rdetg:=sqrt(-detg);
rdetg:=sub({ abs(sin(theta)) = sin(theta), abs(r) = r }, rdetg);

% Euclidean coordinates for hedgehog Anatz
operator xc;
xc(1) := r*sin(theta)*cos(phi);
xc(2) := r*sin(theta)*sin(phi);
xc(3) := r*cos(theta);

% Hedgehog scalar field Ansatz
array phi1(3);
depend u,r;
for a:=1:3 do phi1(a) := xc(a)*u/r;

% Hedgehog gauge field Ansatz
array A2d(3,spdim), A2cd(3,spdim);

% need Levi-Civita symbol
filleps3();

depend a,r;

% Ansatz in cartesian coordinates
for a1:=1:3 do for j:=1:spdim do A2cd(a1,j) := for b:=1:3 sum eps3(j,a1,b) * xc(b) * ( a - 1/gg/r^2 );

% transform to polar
array Tpc(3,3);
for j:=1:3 do for k:=1:3 do Tpc(j,k) := df( xc(j), x(k));

for j:=1:3 do for a1:=1:3 do A2d(a1,j) := for k:=1:3 sum A2cd(a1,k) * Tpc(k,j);

settau();
comcomm2();

comfieldstrength2();
%setcovd2();



lagscalar := 1/2*( for j:=0:spdim sum for k:=0:spdim sum for a:=1:3 sum Dca2(phi1,a,j) * Dca2(phi1,a,k) * guu(j,k) )
    + mu^2/2 * ( for a:=1:3 sum phi1(a)^2 ) - lambda1 / 4 * ( for a:=1:3 sum phi1(a)^2 )^2 $ length(ws);


lagscalar := trigsimp(lagscalar);

laggauge := -1/4 * ( for a1:=1:3 sum for k:=0:spdim sum for l:=0:spdim sum F2uu(a1,k,l) * F2dd(a1,k,l) ) $ length(ws);

laggauge := trigsimp(laggauge);

lag := lagscalar + laggauge;

for all f,x let df(f,x,f) = df(f,f,x);

begin scalar lagp, dldpiu, dldpia;
    lagp := rdetg * sub( { df(u, r) = piu, df(a, r) = pia }, lag);
    dldpiu := df(lagp, piu);
    dldpia := df(lagp, pia);
    dldpiu := sub({ piu = df(u, r), pia = df(a,r) }, dldpiu);
    dldpia := sub({ piu = df(u, r), pia = df(a,r) }, dldpia);

    lagp := sub({ piu = df(u, r), pia = df(a,r) }, lagp);

    eomu := (1/rdetg) * ( df(dldpiu, r) - df(lagp, u) );
    eoma := (1/rdetg) * ( df(dldpia, r) - df(lagp, a) );

    eomu := eomu / coeffn(eomu, df(u,r,2),1);
    eoma := eoma / coeffn(eoma, df(a,r,2),1);
end;

eomup := (1/r^2)* df( r^2*df(u,r), r) + (mu^2 -2 * gg^2 * a^2 * r^2 )*u - lambda1*u^3;
%
% Note: r^2 in the term 2 * gg^2 * a^2 * r^2 was not there in the original paper !


eomu - eomup;

eomap := df(a,r,2) + 4/r * df(a,r) + 3/r^2 * a - gg^2 * r^2 * a^3 - gg^2 u^2 * a;
% Note: coeffn of a^2/r^2 was - 3 in the original paper!

eoma - eomap;

clear eomup, eomap;

% Verify consistency of the Ansatz
% scalar
array dphi1d(3,spdim);
for a:=1:3 do for j:=0:spdim do dphi1d(a,j) := Dca2(phi1, a, j);

array dphi1u(3, spdim);
for a:=1:3 do for j:=0:spdim do dphi1u(a,j) := for k:=0:spdim sum guu(j,k)*dphi1d(a,k);

comchristoffel();

array ddphi1(3);
for a:=1:3 do ddphi1(a) := for j:=0:spdim sum df( dphi1u(a,j), x(j)) + ( for k:=0:spdim sum Gaudd(k,j,k)*dphi1u(a,j) )
    + gg * for b:=1:3 sum for c:=1:3 sum c2(a,b,c)*A2d(b,j)*dphi1u(c,j);

array eomphi1(3);
for a:=1:3 do eomphi1(a) := ddphi1(a) - mu^2 * phi1(a) + lambda1 * ( for b:=1:3 sum phi1(b)^2 ) * phi1(a);
for a:=1:3 do eomphi1(a) := trigsimp(eomphi1(a));

for a:=1:3 do begin scalar qq;
    qq := eomphi1(a) - coeffn( eomphi1(a), df(u,r,2), 1)/coeffn( eomu, df(u,r,2), 1)*eomu;
    if qq neq 0 then write "Error in scalar field equation compt j = ",j," : ", qq;
end;
write "All scalar field equation components verified.";

% gauge
array dF2u(3,spdim);
for a:=1:3 do for j:=0:spdim do dF2u(a,j) := for k:=0:spdim sum df( F2uu(a,j,k), x(k) ) + ( for l:=0:spdim sum Gaudd(k,l,k)*F2uu(a,j,l) + Gaudd(j,k,l)*F2uu(a,k,l) )
    + gg * for b:=1:3 sum for c:=1:3 sum c2(a,b,c) * A2d(b,k) * F2uu(c,j,k);

for a:=1:3 do for j:=0:spdim do dF2u(a,j) := trigsimp(df2u(a,j));

array j2u(3,spdim);
for a:=1:3 do for j:=0:spdim do j2u(a,j) := for b:=1:3 sum for c:=1:3 sum gg * dphi1u(b,j)*c2(b,a,c)*phi1(c);
for a:=1:3 do for j:=0:spdim do j2u(a,j) := trigsimp(j2u(a,j));

array eoma2u(3,spdim);
for a:=1:3 do for j:=0:spdim do eoma2u(a,j) := df2u(a,j) - j2u(a,j);

for a1 :=1:3 do for j:=0:spdim do begin scalar qq;
    qq := eoma2u(a1,j) - coeffn(eoma2u(a1,j), df(a,r,2),1)/coeffn(eoma, df(a,r,2),1)*eoma;
    if qq neq 0 then write "Error in gauge field equation component a = ", a1, ", j = ", j,": ",qq;
end;
write "All gauge field equation components verified.";

end;
