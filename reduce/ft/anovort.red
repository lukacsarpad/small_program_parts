%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abrikosov-Nielsen-Olesen vortex equations                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

in "../grfunc.red";
in "../ftfunc.red";

% 2 d
spdim := 2;

% flat metric, polar coordinates
operator x;
x(0) := t;
x(1) := r;
x(2) := phi;

setgpolar();

% inverse metric
comguu();

rdetg := sqrt(detg);

rdetg := sub(abs(r)=r,rdetg);

% Electromagnetic interaction only
ggp := 2*e;

% Ansatz
% gauge field, radial Ansatz
array A1d(2);
depend a,r;
a1d(2) := n*a;
% scalar field
depend f,r;
h := f*exp(i*n*phi);

realvalued f, phi, n, a, r;
% calculate field strength tensor
comfieldstrength1();

% lagrangian
lagem := for k:=0:2 sum for l:= 0:2 sum -F1dd(k,l)*F1uu(k,l)/4;

lagscalar := ( for k:=0:2 sum for l:=0:2 sum guu(k,l)*conj(Dc1(h,k))*Dc1(h,l) ) - la*( conj(h) * h - eta^2)^2 $

lagscalar := trigsimp(lagscalar,expon);

lag := lagem + lagscalar;

for all f,x let df(f,x,f) = df(f,f,x);

realvalued pia, pif;
begin scalar lagp, dldpia, dldpih;
    lagp := rdetg * sub({ df(a,r) = pia, df(f,r) = pif }, lag);
    dldpia := df(lagp, pia);
    dldpif := df(lagp, pif);
    dldpia := sub({ pia = df(a,r), pif = df(f,r) }, dldpia);
    dldpif := sub({ pia = df(a,r), pif = df(f,r) }, dldpif);
    lagp := rdetg * lag;
    eoma := (1/rdetg) * ( df(dldpia, r) - df(lagp, a) );
    eomf := (1/rdetg) * ( df(dldpif, r) - df(lagp, f) );
end;

eoma := - eoma / coeffn(eoma, df(a,r,2),1);

eomf := - eomf / coeffn(eomf, df(f,r,2),1);

eomfp := -1/r * df( r*df(f,r),r) + f*( n^2*(1-e*a)^2/r^2 + 2*la*(f^2-eta^2)) $

eomf-eomfp;

eomap := -r * df( df(a,r)/r,r) + 2*e*f^2*(e*a-1) $

eoma - eomap;

clear(eomfp, eomap);

% verify consistency of the Ansatz
array dhd(spdim);
for j:=0:spdim do dhd(j) := df( h, x(j) ) - i*e*A1d(j)*h;

array dhu(spdim);
for j:=0:spdim do dhu(j) := for k:=0:spdim sum guu(j,k)*dhd(k);

comchristoffel();

array A1u(spdim);
for j:=0:spdim do A1u(j) := for k:=0:spdim sum guu(j,k)*A1d(k);

ddh := for j:=0:spdim sum df( dhu(j), x(j)) + ( for k:=0:spdim sum Gaudd(j,k,j) * dhu(k) )
	- i*e*A1u(j)*dhd(j) ;

eoms := ddh + 2*la*( conj(h) * h - eta^2)*h $

eoms := trigsimp( eoms, expon) ;

% verify scalar eom: shall be zero
eoms - coeffn(eoms,df(f,r,2),1)/coeffn(eomf,df(f,r,2),1)*eomf;

array dF1u(spdim);
for j:=0:spdim do dF1u(j) := for k:=0:spdim sum df( F1uu(j,k), x(k) ) + for l:=0:spdim sum F1uu(j,k)*Gaudd(l,k,l)
    + F1uu(l,k)*Gaudd(j,l,k);

array j1d(spdim);
for j:=0:spdim do j1d(j) := i*e*( conj(h)*dhd(j) - h*conj(dhd(j)) );
for j:=0:spdim do j1d(j) := trigsimp( j1d(j), expon);

array j1u(spdim);
for j:=0:spdim do j1u(j) := for k:=0:spdim sum guu(j,k)*j1d(k);

array eomau(spdim);
for j:=0:spdim do eomau(j) := df1u(j) - j1u(j);

for j:=0:spdim do begin scalar qq;
    qq := eomau(j) - coeffn(eomau(j), df(a,r,2), 1)/coeffn(eoma, df(a,r,2), 1) * eoma;
    if qq neq 0 then write "Error in gauge field equation, index j = ", j, " : ", qq;
end;
write "All gauge field equation components verified.";

end;
