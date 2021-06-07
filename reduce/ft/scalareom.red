%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scalar field equation of motion                                      %
% O(3) symmetric real scalar                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

in "../grfunc.red";

operator x;
operator gdd;

% spatial dimensions
spdim := 3;

% symmetric metric
for j:=0:spdim do for k:=j+1:spdim do gdd(k,j) := gdd(j,k);

for i1 := 0:spdim do for i2:=i1:spdim do for i3:=0:spdim do depend gdd(i1,i2), x(i3);

comguu();

rdetg := sqrt(-detg) $ length(ws) ;

% scalar field
operator phi1;

for i1 :=0:spdim do for a1:=1:3 do depend phi1(a1),x(i1);

phi2 := for a1:=1:3 sum phi1(a1)*phi1(a1);

for all xi2 let V(xi2) = - mu^2/2 * xi2 + lambda1/4 * xi2^2;

% SO(3) symmetric scalar Lagrangian
lag := 1/2 * ( for j:=0:spdim sum for k:=0:spdim sum for a:=1:3 sum guu(j,k)*df(phi1(a),x(j))*df(phi1(a),x(k)) )
	- V(phi2) $ length(ws);

% calculate Euler-Lagrange equations
array eomhpi(3);
operator dphi;

lagp := rdetg*lag $ length(ws);
lagps := lagp $
for a:=1:3 do for j:=0:spdim do lagps := sub( df(phi1(a),x(j)) = dphi(a,j) , lagps);

array piphi(3,spdim);

for all f,x let df(f,x,f) = df(f,f,x);

for a:=1:3 do for j:=0:spdim do piphi(a,j) := df( lagps, dphi(a,j));

for a1:=1:3 do for j1:=0:spdim do for a2:=1:3 do for j2:=0:spdim do
	    piphi(a1,j1) := sub( dphi(a2,j2) = df(phi1(a2),x(j2) ), piphi(a1,j1) );

array eomphi(3);

for a:=1:3 do eomphi(a) := ( ( for j := 0:spdim sum df( piphi(a,j), x(j) ) ) - df(lagp, phi1(a)) ) / rdetg;

%% save result to file
%off nat;
%off echo;
%out "scalareom.out";
%for a:=1:3 do write "eomphi(",a,") := ", eomphi(a), "$";
%write ";end;";
%shut "scalareom.out";
%on nat;
%on echo;

% known formula for scalar EoM
array dphid(3,spdim);
for a:=1:3 do for j:=0:spdim do dphid(a,j) := df(phi1(a), x(j));

array dphiu(3,spdim);
for a:=1:3 do for j:=0:spdim do dphiu(a,j) := for k:=0:spdim sum guu(j,k)*dphid(a,k);

comchristoffel();

array ddphi(3);
for a:=1:3 do ddphi(a) := for j:=0:spdim sum df( dphiu(a,j), x(j)) + for l:=0:spdim sum dphiu(a,j)*gaudd(l,j,l);

array eomphi2(3);
for a:=1:3 do eomphi2(a) := ddphi(a) - mu^2*phi1(a) + lambda1 * ( for b:=1:3 sum phi1(b)*phi1(b) ) * phi1(a);

%% save results to file
%off nat;
%off echo;
%out "scalareom2.out";
%%write "rdg := ", rdetg, "$";
%for a:=1:3 do write "eomphi(",a,") := ", eomphi2(a), "$";
%write ";end;";
%shut "scalareom2.out";
%on nat;
%on echo;

% compare Euler-Lagrange equations and the hand-typed version
for a:=1:3 do begin scalar qq;
    qq := eomphi(a) - eomphi2(a);
    if qq neq 0 then write "Error in scalar EoM, a = ", a;
end;

write "Scalar equation of motion: variation of Lagrangian and hand-typed formula compared.";

end;
