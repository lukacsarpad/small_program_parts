%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic monopole radial equations                                   %
%                                                                      %
% References:                                                          %
%   * Breitenlohner, Forgács, Maison, Nucl Phys B383, 357-376 (1992)   %
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

% Euclidean coordinates for hedgehog Anatz
operator xc;
xc(1) := r*sin(theta)*cos(phi);
xc(2) := r*sin(theta)*sin(phi);
xc(3) := r*cos(theta);

% Hedgehog scalar field Ansatz
array phi1(3);
depend h,r;
for a:=1:3 do phi1(a) := xc(a)*h/r;

% Hedgehog gauge field Ansatz
array A2d(3,spdim), A2cd(3,spdim);

% need Levi-Civita symbol
filleps3();

depend w,r;

% Ansatz in cartesian coordinates
for a1:=1:3 do for j:=1:spdim do A2cd(a1,j) := for b:=1:3 sum eps3(j,a1,b) * xc(b) * (w - 1)/r^2/gg;

% transform to polar
array Tpc(3,3);
for j:=1:3 do for k:=1:3 do Tpc(j,k) := df( xc(j), x(k));

for j:=1:3 do for a1:=1:3 do A2d(a1,j) := for k:=1:3 sum A2cd(a1,k) * Tpc(k,j);

settau();
comcomm2();

comfieldstrength2();
%setcovd2();



lagscalar := 1/2*( for j:=0:spdim sum for k:=0:spdim sum for a:=1:3 sum Dca2(phi1,a,j) * Dca2(phi1,a,k) * guu(j,k) )
    - lambda1/8 * ( (for a:=1:3 sum phi1(a)^2 ) -H0^2 )^2 $ length(ws);

lagscalar := trigsimp(lagscalar);

laggauge := -1/4 * ( for a1:=1:3 sum for k:=0:spdim sum for l:=0:spdim sum F2uu(a1,k,l) * F2dd(a1,k,l) ) $ length(ws);

laggauge := trigsimp(laggauge);

lag := lagscalar + laggauge;

for all f,x let df(f,x,f) = df(f,f,x);

begin scalar lagp, dldpih, dldpiw;
    lagp := rdetg * sub( { df(h, r) = pih, df(w, r) = piw }, lag);
    dldpih := df(lagp, pih);
    dldpiw := df(lagp, piw);
    dldpih := sub({ pih = df(h, r), piw = df(w,r) }, dldpih);
    dldpiw := sub({ pih = df(h, r), piw = df(w,r) }, dldpiw);

    lagp := sub({ pih = df(h, r), piw = df(w,r) }, lagp);

    eomh := (1/rdetg) * ( df(dldpih, r) - df(lagp, h) );
    eomw := (1/rdetg) * ( df(dldpiw, r) - df(lagp, w) );

    eomh := eomh / coeffn(eomh, df(h,r,2),1);
    eomw := eomw / coeffn(eomw, df(w,r,2),1);
end;

eomhp := 1/r^2*df( r^2*df(h,r),r) - 2*w^2*h/r^2 - lambda1/2*(h^2-h0^2)*h;

eomh - eomhp;

eomwp := df(w,r,2) - w*( (w^2-1)/r^2 + gg^2 * h^2);

eomw-eomwp;


clear eomhp, eomwp;

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
for a:=1:3 do eomphi1(a) := ddphi1(a) + lambda1/2 * ( -h0^2 + for b:=1:3 sum phi1(b)^2 ) * phi1(a);
for a:=1:3 do eomphi1(a) := trigsimp(eomphi1(a));

for a:=1:3 do begin scalar qq;
    qq := eomphi1(a) - coeffn( eomphi1(a), df(h,r,2), 1)/coeffn( eomh, df(h,r,2), 1)*eomh;
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
    qq := eoma2u(a1,j) - coeffn(eoma2u(a1,j), df(w,r,2),1)/coeffn(eomw, df(w,r,2),1)*eomw;
    if qq neq 0 then write "Error in gauge field equation component a = ", a1, ", j = ", j,": ",qq;
end;
write "All gauge field equation components verified.";

end;
