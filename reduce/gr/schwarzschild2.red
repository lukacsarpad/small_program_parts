%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schwarzschild spacetime                                               %
% check metric                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% packages needed
load_package "taylor";


% coordinates
operator x;

x(0) := t;
x(1) := r;
x(2) := theta;
x(3) := phi;

% metric
array gdd(3,3), guu(3,3);
gdd(0,0) := (1-rg/r);
gdd(1,1) := -1/(1-rg/r);
gdd(2,2) := -r^2;
gdd(3,3) := -r^2*sin(theta)^2;

% diag
for k:=0:3 do guu(k,k) := 1/gdd(k,k);
rdetg:=sqrt(- for k:=0:3 product gdd(k,k) );

array eta(3,3);
eta(0,0):=1;
eta(1,1):=eta(2,2):=eta(3,3):=-1;

%%%%%%%%%%%%%%%%%%%%%%%%% Levi-Civita symbol %%%%%%%%%%%%%%%%%%%%%%%%%%%%
array eps4uuuu(3,3,3,3), eps4dddd(3,3,3,3);
eps4uuuu(0,1,2,3):=eps4uuuu(0,2,3,1):=eps4uuuu(0,3,1,2):=1 $
eps4uuuu(0,3,2,1):=eps4uuuu(0,2,1,3):=eps4uuuu(0,1,3,2):=-1 $

eps4uuuu(1,0,2,3):=eps4uuuu(2,1,0,3):=eps4uuuu(3,1,2,0):=-1 $
eps4uuuu(2,0,3,1):=eps4uuuu(3,2,0,1):=eps4uuuu(1,2,3,0):=-1 $
eps4uuuu(3,0,1,2):=eps4uuuu(1,3,0,2):=eps4uuuu(2,3,1,0):=-1 $

eps4uuuu(3,0,2,1):=eps4uuuu(2,3,0,1):=eps4uuuu(1,3,2,0):=1 $
eps4uuuu(2,0,1,3):=eps4uuuu(1,2,0,3):=eps4uuuu(3,2,1,0):=1 $
eps4uuuu(1,0,3,2):=eps4uuuu(3,1,0,2):=eps4uuuu(2,1,3,0):=1 $

for i1:=0:3 do for i2:=0:3 do for i3:=0:3 do for i4:=0:3 do
    eps4dddd(i1,i2,i3,i4) := for j1:=0:3 sum for j2:=0:3 sum for j3:=0:3 sum for j4:=0:3 sum
	eta(i1,j1)*eta(i2,j2)*eta(i3,j3)*eta(i4,j4)*eps4uuuu(j1,j2,j3,j4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Volume form %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array E4dddd(3,3,3,3), Euuuu(3,3,3,3);

for i1:=0:3 do for i2:=0:3 do for i3:=0:3 do for i4:=0:3 do E4dddd(i1,i2,i3,i4) := eps4dddd(i1,i2,i3,i4)*rdetg;

for i1:=0:3 do for i2:=0:3 do for i3:=0:3 do for i4:=0:3 do Euuuu(i1,i2,i3,i4) := eps4uuuu(i1,i2,i3,i4)/rdetg;

%%%%%%%%%%%%%%%%%%%%%%%%% Christoffel symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%
array gaddd(3,3,3), gaudd(3,3,3);

for k:=0:3 do for l:=0:3 do for m:=0:3 do gaddd(k,l,m) := 1/2 * ( df(gdd(k,l),x(m)) + df(gdd(k,m),x(l)) - df(gdd(l,m),x(k)) );

for k:=0:3 do for l:=0:3 do for m:=0:3 do gaudd(k,l,m) := for n:=0:3 sum guu(k,n)*gaddd(n,l,m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Riemann tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%
array Rieuddd(3,3,3,3), Riedddd(3,3,3,3);

for k:=0:3 do for l:=0:3 do for m:=0:3 do for n:=0:3 do Rieuddd(k,l,m,n) := df(gaudd(k,l,n),x(m)) - df(gaudd(k,l,m),x(n))
		+ for j:=0:3 sum ( gaudd(k,j,m)*gaudd(j,l,n) - gaudd(k,j,n)*gaudd(j,l,m) );

for k:=0:3 do for l:=0:3 do for m:=0:3 do for n:=0:3 do Riedddd(k,l,m,n):=for j:=0:3 sum gdd(k,j)*Rieuddd(j,l,m,n);

% verify symmetry

for k:=0:3 do for l:=0:3 do for m:=0:3 do for n:=0:3 do begin scalar qq;
    qq := Riedddd(k,l,m,n)+Riedddd(l,k,m,n);
    if qq neq 0 then write "Error in Riemann tensor antisymmetry, first pair, indices: ", k,l,m,n;
end;

for k:=0:3 do for l:=0:3 do for m:=0:3 do for n:=0:3 do begin scalar qq;
    qq := Riedddd(k,l,m,n)+Riedddd(k,l,n,m);
    if qq neq 0 then write "Error in Riemann tensor antisymmetry, first pair, indices: ", k,l,m,n;
end;

for k:=0:3 do for l:=0:3 do for m:=0:3 do for n:=0:3 do begin scalar qq;
    qq := Riedddd(k,l,m,n)-Riedddd(m,n,k,l);
    if qq neq 0 then write "Error in Riemann tensor antisymmetry, first pair, indices: ", k,l,m,n;
end;

% verify cyclic

for k:=0:3 do for l:=0:3 do for m:=0:3 do for n:=0:3 do begin scalar qq;
    qq := Riedddd(k,l,m,n)+Riedddd(k,n,l,m)+Riedddd(k,m,n,l);
    if qq neq 0 then write "Error in Riemann tensor cyclic property, first pair, indices: ", k,l,m,n;
end;

%%%%%%%%%%%%%%%%%%%%%%%% Verify Bianchi identity %%%%%%%%%%%%%%%%%%%%%%%%
array DRieudddd(3,3,3,3,3);

for n:=0:3 do for j:=0:3 do for k:=0:3 do for l:=0:3 do for m:=0:3 do DRieudddd(n,j,k,l,m) := df( Rieuddd(n,j,k,l), x(m))
    + for p:=0:3 sum ( gaudd(n,p,m)*Rieuddd(p,j,k,l) - gaudd(p,j,m)*Rieuddd(n,p,k,l) - gaudd(p,k,m)*Rieuddd(n,j,p,l)
						     - gaudd(p,l,m)*Rieuddd(n,j,k,p) );

for n:=0:3 do for j:=0:3 do for k:=0:3 do for l:=0:3 do for m:=0:3 do begin scalar qq;
    qq := DRieudddd(n,j,k,l,m) + DRieudddd(n,j,m,k,l) + DRieudddd(n,j,l,m,k);
    if qq neq 0 then write "Error in Bianchi identity, indices: ", n,j,k,l,m;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ricci tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%
array Ricdd(3,3);

for k:=0:3 do for l:=0:3 do Ricdd(k,l) := for m:=0:3 sum Rieuddd(m,k,m,l);

% verify formula, 92.7
array Ricddp(3,3);

for j:=0:3 do for k:=0:3 do Ricddp(j,k) := for l:=0:3 sum ( df(gaudd(l,j,k),x(l)) - df(gaudd(l,j,l),x(k))
		+ for m:=0:3 sum ( gaudd(l,j,k)*gaudd(m,l,m) - gaudd(m,j,l)*gaudd(l,k,m) ) );

for k:=0:3 do for l:=0:3 do begin scalar qq;
    qq := Ricdd(k,l) - Ricddp(k,l);
    if qq neq 0 then write "Error in Ricci tensor, indices: ", k, l;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ricci scalar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rs := for k:=0:3 sum for l:=0:3 sum guu(k,l)*Ricdd(l,k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Einstein tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%
array Eidd(3,3), Eiud(3,3);

for k:=0:3 do for l:=0:3 do Eidd(k,l) := Ricdd(k,l) - 1/2 * gdd(k,l)*Rs;

for j:=0:3 do for k:=0:3 do Eiud(j,k) := for l:=0:3 sum guu(j,l)*Eidd(l,k);

% verify vacuum solution
for k:=0:3 do for l:=0:3 do if Eidd(k,l) neq 0 then write "Error in Einstein equations, indices: ", k, l;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Weyl tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array Wedddd(3,3,3,3);

for j:=0:3 do for k:=0:3 do for l:=0:3 do for m:=0:3 do Wedddd(j,k,l,m) := Riedddd(j,k,l,m)
    - 1/2 * Ricdd(j,l)*gdd(k,m) + 1/2 * Ricdd(j,m)*gdd(k,l) + 1/2*Ricdd(k,l)*gdd(j,l)
    + 1/6 * Rs * ( gdd(j,l)*gdd(k,m) - gdd(j,m)*gdd(k,l) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Asymptotics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newtonian potential
npot := k * mass;

% asymptotic form of gdd(0,0)
sub(rp=r,sub(r=rp/epsilon,gdd(0,0))) $

taylor(ws,epsilon,0,2);

% note: gdd(0,0) ~ 1 + 2 phi / c^2, where phi is the Newtonian gravitational potential
% solve for Schwarzschild radius

solve(coeffn(-r*taylortostandard(ws),epsilon,1)=2*npot/c^2,rg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Invariants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array Rieuudd(3,3,3,3), Rieuuuu(3,3,3,3), rsdddd(3,3,3,3), rsuuuu(3,3,3,3);

for i1:=0:3 do for i2:=0:3 do for i3:=0:3 do for i4:=0:3 do Rieuudd(i1,i2,i3,i4) :=
	for j2:=0:3 sum Rieuddd(i1,j2,i3,i4)*guu(i2,j2);

for i1:=0:3 do for i2:=0:3 do for i3:=0:3 do for i4:=0:3 do Rieuuuu(i1,i2,i3,i4) :=
	for j3:=0:3 sum for j4:=0:3 sum Rieuudd(i1,i2,j3,j4)*guu(i3,j3)*guu(i4,j4);

% dual of the Riemann tensor
for i1:=0:3 do for i2:=0:3 do for i3:=0:3 do for i4:=0:3 do rsdddd(i1,i2,i3,i4) :=
	for j1:=0:3 sum for j2:=0:3 sum 1/2*E4dddd(i1,i2,j1,j2)*Rieuudd(j1,j2,i3,i4);

for i1:=0:3 do for i2:=0:3 do for i3:=0:3 do for i4:=0:3 do rsuuuu(i1,i2,i3,i4) :=
	for j1:=0:3 sum for j2:=0:3 sum for j3:=0:3 sum for j4:=0:3 sum
	    rsdddd(j1,j2,j3,j4)*guu(i1,j1)*guu(i2,j2)*guu(i3,j3)*guu(i4,j4);

for i1:=0:3 do for i2:=0:3 do for i3:=0:3 do for i4:=0:3 do rsdduu(i1,i2,i3,i4) :=
	for j3:=0:3 sum for j4:=0:3 sum rsdddd(i1,i2,j3,j4)*guu(i3,j3)*guu(i4,j4);

Invar1 := for j:=0:3 sum for k:=0:3 sum for l:=0:3 sum for m:=0:3 sum
    1/48*( Riedddd(j,k,l,m) * Rieuuuu(j,k,l,m) - i * Riedddd(j,k,l,m)*rsuuuu(j,k,l,m));

Invar2 := for j:=0:3 sum for k:=0:3 sum for l:=0:3 sum for m:=0:3 sum for p:=0:3 sum for s:=0:3 sum
    1/96*( Riedddd(j,k,l,m)*Rieuuuu(l,m,p,s)*Rieuudd(j,k,p,s) + i * Riedddd(j,k,l,m)*Rieuuuu(l,m,p,s)*rsdduu(p,s,j,k) );

% nonsingular at r=rg
sub(r=rg,Invar1);

sub(r=rg,Invar2);

% singularity at r=0:
limit(Invar1,r,0);

limit(Invar2,r,0);

solve({(lambda1^2+lambda2^2+lambda1*lambda2)/3=invar1,lambda1*lambda2*(lambda1+lambda2)/2=invar2},{lambda1,lambda2});


end;