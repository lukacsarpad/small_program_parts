%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravitating magnetic monopole radial equations                       %
%                                                                      %
% References:                                                          %
%   * Breitenlohner, Forgács, Maison, Nucl Phys B383, 357-376 (1992)   %
%   * Breitenlohner, Maison, Lavrelashvili, CQG 21, 1667-1684 (2004)   %
%          letter for the constants in the gravitational action        %
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

array gdd(3,3);
depend a, r;
depend mu,r;
gdd(0,0) := a^2*mu;
gdd(1,1) := -1/mu;
gdd(2,2) := -r^2;
gdd(3,3) := -r^2*sin(theta)^2;

comguu3();

rdetg:=sqrt(-detg);
rdetg:=sub({ abs(sin(theta)) = sin(theta), abs(r) = r, abs(a)=a }, rdetg);

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

% scalar lagrangian
lagscalar := 1/2*( for j:=0:spdim sum for k:=0:spdim sum for a:=1:3 sum Dca2(phi1,a,j) * Dca2(phi1,a,k) * guu(j,k) )
    - lambda1/8 * ( (for a:=1:3 sum phi1(a)^2 ) -H0^2 )^2 $ length(ws);

lagscalar := trigsimp(lagscalar);

% gauge lagrangian
comfieldstrength2();

laggauge := -1/4 * ( for a1:=1:3 sum for k:=0:spdim sum for l:=0:spdim sum F2uu(a1,k,l) * F2dd(a1,k,l) ) $ length(ws);

laggauge := trigsimp(laggauge);

% calculate gravitational quantities
comchristoffel();
comriemann();

riemsymm();
bianchialg();
bianchidiff();

comricci();

laggrav := -Rs/Gr/4;
% Note sign and constants!

lag := lagscalar + laggauge + laggrav;

for all f,x let df(f,x,f) = df(f,f,x);

for all f,x,y let df(f,x,y,f) = df(f,f,y,x);

begin scalar lagp, dldpih, dldpiw, dldpia, dldpimu, eommu0;
    lagp := rdetg * sub( { df(h, r) = pih, df(w, r) = piw , df(a,r) = pia, df(mu,r) = pimu, df(a,r,2) = pi2a, df(mu,r,2) = pi2mu }, lag);
    dldpih  := df(lagp, pih);
    dldpiw  := df(lagp, piw);
    dldpia  := df(lagp, pia);
    dldpi2a  := df(lagp, pi2a);
    dldpimu := df(lagp, pimu);
    dldpi2mu := df(lagp, pi2mu);
    dldpih  := sub({ pih = df(h, r), piw = df(w,r), pia = df(a,r), pimu = df(mu, r), pi2a = df(a,r,2), pi2mu = df(mu,r,2) }, dldpih);
    dldpiw  := sub({ pih = df(h, r), piw = df(w,r), pia = df(a,r), pimu = df(mu, r), pi2a = df(a,r,2), pi2mu = df(mu,r,2) }, dldpiw);
    dldpia  := sub({ pih = df(h, r), piw = df(w,r), pia = df(a,r), pimu = df(mu, r), pi2a = df(a,r,2), pi2mu = df(mu,r,2) }, dldpia);
    dldpimu := sub({ pih = df(h, r), piw = df(w,r), pia = df(a,r), pimu = df(mu, r), pi2a = df(a,r,2), pi2mu = df(mu,r,2) }, dldpimu);

    lagp := sub({ pih = df(h, r), piw = df(w,r), pia = df(a,r), pimu = df(mu, r), pi2a = df(a,r,2), pi2mu = df(mu,r,2) }, lagp);

    eomh  := (1/rdetg) * ( df(dldpih, r) - df(lagp, h) );
    eomw  := (1/rdetg) * ( df(dldpiw, r) - df(lagp, w) );
    eoma  := (1/rdetg) * ( -df(dldpi2a,r,2) + df(dldpia, r) - df(lagp, a) );
    eommu := (1/rdetg) * ( -df(dldpi2mu,r,2) + df(dldpimu, r) - df(lagp, mu) );


    eomh  := eomh / coeffn(eomh, df(h,r,2),1);
    eomw  := eomw / coeffn(eomw, df(w,r,2),1);
    eommu0 := eoma;
    eoma  := eommu / coeffn(eommu, df(a,r),1);
    eommu := eommu0 / coeffn(eommu0, df(mu,r),1);
end;


eomhp := 1/r^2*df( r^2*a*mu*df(h,r),r) - 2*a*w^2*h/r^2 - a*lambda1/2*(h^2-h0^2)*h $

a*mu*eomh - eomhp;

eomwp := df(a*mu*df(w,r),r) - a*w*( (w^2-1)/r^2 + gg^2 * h^2) $

a*mu*eomw-eomwp;

V1 := df(w,r)^2/gg^2 + r^2 * df(h,r)^2/2 $
eomap := df(a,r) - 2* Gr * V1 * a / r $

eoma - eomap;

V2 := (1-w^2)^2/2/gg^2/r^2 + w^2*h^2 + lambda1/8 * r^2 * (h^2-h0^2)^2 $

eommup := df(mu,r) - 1/r*( 1 - mu - 2*Gr*( mu*V1 + V2)) $

eommu - eommup;

clear eomhp, eomwp, eomap, eommup;

% Verify consistency of the Ansatz
% scalar
array dphi1d(3,spdim);
for a:=1:3 do for j:=0:spdim do dphi1d(a,j) := Dca2(phi1, a, j);

array dphi1u(3, spdim);
for a:=1:3 do for j:=0:spdim do dphi1u(a,j) := for k:=0:spdim sum guu(j,k)*dphi1d(a,k);

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

% gravity
comeinstein();
for j:=0:spdim do for k:=0:spdim do <<
    Eidd(j,k) := trigsimp(Eidd(j,k));
    Eiud(j,k) := trigsimp(Eiud(j,k));
>>;
comT2();
for j:=0:spdim do for k:=0:spdim do <<
    T2dd(j,k) := trigsimp(T2dd(j,k));
    T2ud(j,k) := trigsimp(T2ud(j,k));
    T2uu(j,k) := trigsimp(T2uu(j,k));
>>;

array Tsdd(spdim,spdim), Tsud(spdim,spdim);
% scalar stress-energy tensor, Hilbert form, LL2 94.4
for j:=0:spdim do for k:=0:spdim do Tsdd(j,k) := ( for a:=1:3 sum Dca2(phi1,a,j) * Dca2(phi1,a,k) )
    - gdd(j,k) * lagscalar;
for j:=0:spdim do for k:=0:spdim do Tsud(j,k) := for l:=0:spdim sum guu(j,l)*Tsdd(l,k);

for j:=0:spdim do for k:=0:spdim do <<
    Tsdd(j,k) := trigsimp(Tsdd(j,k));
    Tsud(j,k) := trigsimp(Tsud(j,k));
>>;

array Eiequd(spdim,spdim);
for j:=0:spdim do for k:=0:spdim do Eiequd(j,k) := Eiud(j,k) - 2* Gr * ( T2ud(j,k) + Tsud(j,k) );

for j:=0:spdim do for k:=0:spdim do Eiequd(j,k) := trigsimp(Eiequd(j,k));

% verify
for j:=0:spdim do for k:=0:spdim do begin scalar qq;
    qq := Eiequd(j,k);

    qq := qq - coeffn(qq, df(mu,r,2),1)/coeffn(eommu,df(mu,r),1) * df( eommu, r);
    qq := qq - coeffn(qq, df(a,r,2),1)/coeffn(eoma,df(a,r),1) * df( eoma, r);
    qq := qq - coeffn(qq, df(h,r,2),1)/coeffn(eomh,df(h,r,2),1) * eomh;
    qq := qq - coeffn(qq, df(w,r,2),1)/coeffn(eomw,df(w,r,2),1) * eomw;

    qq := qq - coeffn(qq, df(mu,r), 1)/coeffn(eommu, df(mu,r), 1) * eommu;
    qq := qq - coeffn(qq, df(a,r), 1)/coeffn(eoma, df(a,r), 1) * eoma;
    if qq neq 0 then write "Error in Einstein equations, components ",j," ",k,": ", qq;
end;
write "All components of the Einstein equations verified.";

end;
