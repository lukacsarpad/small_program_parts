%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schwarzschild spacetime                                               %
% check metric                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% packages needed
load_package "taylor";


% coordinates
operator x;

x(1) := r;
x(2) := theta;
x(3) := phi;
x(4) := t;

% metric
array gdd(4,4), guu(4,4);
gdd(1,1) := -1/(1-rg/r);
gdd(2,2) := -r^2;
gdd(3,3) := -r^2*sin(theta)^2;
gdd(4,4) := (1-rg/r);

% diag
for k:=1:4 do guu(k,k) := 1/gdd(k,k);
rdetg:=sqrt(- for k:=1:4 product gdd(k,k) );

array eta(4,4);
eta(1,1):=eta(2,2):=eta(3,3):=-1;
eta(4,4):=1;

%%%%%%%%%%%%%%%%%%%%%%%%% Levi-Civita symbol %%%%%%%%%%%%%%%%%%%%%%%%%%%%
array eps4uuuu(4,4,4,4), eps4dddd(4,4,4,4);
eps4uuuu(4,1,2,3):=eps4uuuu(4,2,3,1):=eps4uuuu(4,3,1,2):=1 $
eps4uuuu(4,3,2,1):=eps4uuuu(4,2,1,3):=eps4uuuu(4,1,3,2):=-1 $

eps4uuuu(1,4,2,3):=eps4uuuu(2,1,4,3):=eps4uuuu(3,1,2,4):=-1 $
eps4uuuu(2,4,3,1):=eps4uuuu(3,2,4,1):=eps4uuuu(1,2,3,4):=-1 $
eps4uuuu(3,4,1,2):=eps4uuuu(1,3,4,2):=eps4uuuu(2,3,1,4):=-1 $

eps4uuuu(3,4,2,1):=eps4uuuu(2,3,4,1):=eps4uuuu(1,3,2,4):=1 $
eps4uuuu(2,4,1,3):=eps4uuuu(1,2,4,3):=eps4uuuu(3,2,1,4):=1 $
eps4uuuu(1,4,3,2):=eps4uuuu(3,1,4,2):=eps4uuuu(2,1,3,4):=1 $

for i1:=1:4 do for i2:=1:4 do for i3:=1:4 do for i4:=1:4 do
    eps4dddd(i1,i2,i3,i4) := for j1:=1:4 sum for j2:=1:4 sum for j3:=1:4 sum for j4:=1:4 sum
	eta(i1,j1)*eta(i2,j2)*eta(i3,j3)*eta(i4,j4)*eps4uuuu(j1,j2,j3,j4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Volume form %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array E4dddd(4,4,4,4), Euuuu(4,4,4,4);

for i1:=1:4 do for i2:=1:4 do for i3:=1:4 do for i4:=1:4 do E4dddd(i1,i2,i3,i4) := eps4dddd(i1,i2,i3,i4)*rdetg;

for i1:=1:4 do for i2:=1:4 do for i3:=1:4 do for i4:=1:4 do Euuuu(i1,i2,i3,i4) := eps4uuuu(i1,i2,i3,i4)/rdetg;

%%%%%%%%%%%%%%%%%%%%%%%%% Christoffel symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%
array gaddd(4,4,4), gaudd(4,4,4);

for k:=1:4 do for l:=1:4 do for m:=1:4 do gaddd(k,l,m) := 1/2 * ( df(gdd(k,l),x(m)) + df(gdd(k,m),x(l)) - df(gdd(l,m),x(k)) );

for k:=1:4 do for l:=1:4 do for m:=1:4 do gaudd(k,l,m) := for n:=1:4 sum guu(k,n)*gaddd(n,l,m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Riemann tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%
array Rieuddd(4,4,4,4), Riedddd(4,4,4,4);

for k:=1:4 do for l:=1:4 do for m:=1:4 do for n:=1:4 do Rieuddd(k,l,m,n) := df(gaudd(k,l,n),x(m)) - df(gaudd(k,l,m),x(n))
		+ for j:=1:4 sum ( gaudd(k,j,m)*gaudd(j,l,n) - gaudd(k,j,n)*gaudd(j,l,m) );

for k:=1:4 do for l:=1:4 do for m:=1:4 do for n:=1:4 do Riedddd(k,l,m,n):=for j:=1:4 sum gdd(k,j)*Rieuddd(j,l,m,n);

% verify symmetry

for k:=1:4 do for l:=1:4 do for m:=1:4 do for n:=1:4 do begin scalar qq;
    qq := Riedddd(k,l,m,n)+Riedddd(l,k,m,n);
    if qq neq 0 then write "Error in Riemann tensor antisymmetry, first pair, indices: ", k,l,m,n;
end;

for k:=1:4 do for l:=1:4 do for m:=1:4 do for n:=1:4 do begin scalar qq;
    qq := Riedddd(k,l,m,n)+Riedddd(k,l,n,m);
    if qq neq 0 then write "Error in Riemann tensor antisymmetry, first pair, indices: ", k,l,m,n;
end;

for k:=1:4 do for l:=1:4 do for m:=1:4 do for n:=1:4 do begin scalar qq;
    qq := Riedddd(k,l,m,n)-Riedddd(m,n,k,l);
    if qq neq 0 then write "Error in Riemann tensor antisymmetry, first pair, indices: ", k,l,m,n;
end;

% verify cyclic

for k:=1:4 do for l:=1:4 do for m:=1:4 do for n:=1:4 do begin scalar qq;
    qq := Riedddd(k,l,m,n)+Riedddd(k,n,l,m)+Riedddd(k,m,n,l);
    if qq neq 0 then write "Error in Riemann tensor cyclic property, first pair, indices: ", k,l,m,n;
end;

%%%%%%%%%%%%%%%%%%%%%%%% Verify Bianchi identity %%%%%%%%%%%%%%%%%%%%%%%%
array DRieudddd(4,4,4,4,4);

for n:=1:4 do for j:=1:4 do for k:=1:4 do for l:=1:4 do for m:=1:4 do DRieudddd(n,j,k,l,m) := df( Rieuddd(n,j,k,l), x(m))
    + for p:=1:4 sum ( gaudd(n,p,m)*Rieuddd(p,j,k,l) - gaudd(p,j,m)*Rieuddd(n,p,k,l) - gaudd(p,k,m)*Rieuddd(n,j,p,l)
						     - gaudd(p,l,m)*Rieuddd(n,j,k,p) );

for n:=1:4 do for j:=1:4 do for k:=1:4 do for l:=1:4 do for m:=1:4 do begin scalar qq;
    qq := DRieudddd(n,j,k,l,m) + DRieudddd(n,j,m,k,l) + DRieudddd(n,j,l,m,k);
    if qq neq 0 then write "Error in Bianchi identity, indices: ", n,j,k,l,m;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ricci tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%
array Ricdd(4,4);

for k:=1:4 do for l:=1:4 do Ricdd(k,l) := for m:=1:4 sum Rieuddd(m,k,m,l);

% verify formula, 92.7
array Ricddp(4,4);

for j:=1:4 do for k:=1:4 do Ricddp(j,k) := for l:=1:4 sum ( df(gaudd(l,j,k),x(l)) - df(gaudd(l,j,l),x(k))
		+ for m:=1:4 sum ( gaudd(l,j,k)*gaudd(m,l,m) - gaudd(m,j,l)*gaudd(l,k,m) ) );

for k:=1:4 do for l:=1:4 do begin scalar qq;
    qq := Ricdd(k,l) - Ricddp(k,l);
    if qq neq 0 then write "Error in Ricci tensor, indices: ", k, l;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ricci scalar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rs := for k:=1:4 sum for l:=1:4 sum guu(k,l)*Ricdd(l,k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Einstein tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%
array Eidd(4,4), Eiud(4,4);

for k:=1:4 do for l:=1:4 do Eidd(k,l) := Ricdd(k,l) - 1/2 * gdd(k,l)*Rs;

for j:=1:4 do for k:=1:4 do Eiud(j,k) := for l:=1:4 sum guu(j,l)*Eidd(l,k);

% verify vacuum solution
for k:=1:4 do for l:=1:4 do if Eidd(k,l) neq 0 then write "Error in Einstein equations, indices: ", k, l;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Weyl tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array Wedddd(4,4,4,4);

for j:=1:4 do for k:=1:4 do for l:=1:4 do for m:=1:4 do Wedddd(j,k,l,m) := Riedddd(j,k,l,m)
    - 1/2 * Ricdd(j,l)*gdd(k,m) + 1/2 * Ricdd(j,m)*gdd(k,l) + 1/2*Ricdd(k,l)*gdd(j,l)
    + 1/6 * Rs * ( gdd(j,l)*gdd(k,m) - gdd(j,m)*gdd(k,l) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Asymptotics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newtonian potential
npot := k * mass;

% asymptotic form of gdd(4,4)
sub(rp=r,sub(r=rp/epsilon,gdd(4,4))) $

taylor(ws,epsilon,0,2);

% note: gdd(4,4) ~ 1 + 2 phi / c^2, where phi is the Newtonian gravitational potential
% solve for Schwarzschild radius

solve(coeffn(-r*taylortostandard(ws),epsilon,1)=2*npot/c^2,rg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Invariants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array Rieuudd(4,4,4,4), Rieuuuu(4,4,4,4), rsdddd(4,4,4,4), rsuuuu(4,4,4,4);

for i1:=1:4 do for i2:=1:4 do for i3:=1:4 do for i4:=1:4 do Rieuudd(i1,i2,i3,i4) :=
	for j2:=1:4 sum Rieuddd(i1,j2,i3,i4)*guu(i2,j2);

for i1:=1:4 do for i2:=1:4 do for i3:=1:4 do for i4:=1:4 do Rieuuuu(i1,i2,i3,i4) :=
	for j3:=1:4 sum for j4:=1:4 sum Rieuudd(i1,i2,j3,j4)*guu(i3,j3)*guu(i4,j4);

% dual of the Riemann tensor
for i1:=1:4 do for i2:=1:4 do for i3:=1:4 do for i4:=1:4 do rsdddd(i1,i2,i3,i4) :=
	for j1:=1:4 sum for j2:=1:4 sum 1/2*E4dddd(i1,i2,j1,j2)*Rieuudd(j1,j2,i3,i4);

for i1:=1:4 do for i2:=1:4 do for i3:=1:4 do for i4:=1:4 do rsuuuu(i1,i2,i3,i4) :=
	for j1:=1:4 sum for j2:=1:4 sum for j3:=1:4 sum for j4:=1:4 sum
	    rsdddd(j1,j2,j3,j4)*guu(i1,j1)*guu(i2,j2)*guu(i3,j3)*guu(i4,j4);

for i1:=1:4 do for i2:=1:4 do for i3:=1:4 do for i4:=1:4 do rsdduu(i1,i2,i3,i4) :=
	for j3:=1:4 sum for j4:=1:4 sum rsdddd(i1,i2,j3,j4)*guu(i3,j3)*guu(i4,j4);

Invar1 := for j:=1:4 sum for k:=1:4 sum for l:=1:4 sum for m:=1:4 sum
    1/48*( Riedddd(j,k,l,m) * Rieuuuu(j,k,l,m) - i * Riedddd(j,k,l,m)*rsuuuu(j,k,l,m));

Invar2 := for j:=1:4 sum for k:=1:4 sum for l:=1:4 sum for m:=1:4 sum for p:=1:4 sum for s:=1:4 sum
    1/96*( Riedddd(j,k,l,m)*Rieuuuu(l,m,p,s)*Rieuudd(j,k,p,s) + i * Riedddd(j,k,l,m)*Rieuuuu(l,m,p,s)*rsdduu(p,s,j,k) );

% nonsingular at r=rg
sub(r=rg,Invar1);

sub(r=rg,Invar2);

% singularity at rg:
limit(Invar1,r,0);

limit(Invar2,r,0);

solve({(lambda1^2+lambda2^2+lambda1*lambda2)/3=invar1,lambda1*lambda2*(lambda1+lambda2)/2=invar2},{lambda1,lambda2});


end;