%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General relativity specific calculations                              %
%                                                                       %
% Preliminaries:                                                        %
%   * Set spdim to number of spatial dimensions                         %
%   * Set up operator x, x(i)=coordinates, i=0:spdim                    %
%   * Unless using spherical or cartesian, set gdd() to metric          %
%                                                                       %
% References:                                                           %
%       Landau & Lifshitz 2: coordinate basis formulae                  %
%       MacCallum & McCrea: grlib                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global '(oldecho!*) $
symbolic << oldecho!* := !*echo; off1 'echo >>;  % or use on1 to turn echo on.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inverse metric %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation of guu, a modified version of McCrea and MacCallum's comguu.red
% Coordinate basis
%                                 ij
%This file computes GUU(I,J)<--->g  , i.e. inverse of GDD(I,J)<--->g
%                                                                   ij
%must be preceded by file inputting gdd(i,j).

procedure comguu3; begin;
    array ag(2,2), bg(3,3),guu(3,3);
    if spdim neq 3 then write "COMGUU3 error: spdim not 3!";
    for i:=0:3 do for j:=i+1:3 do if gdd(j,i) neq gdd(i,j) then write "COMGUU3 error: gdd(i,j) not symmetric for indices ", i, ", ", j;
    for k:=0:3 do for l:=k:3 do begin
	integer ii;
	for i:=0:2 do begin
	    if i>=k then ii:=i+1 else ii:=i;
	    for j:=0:2 do if j>=l then ag(i,j):=gdd(ii,j+1) else ag(i,j):=gdd(ii,j) end;
	bg(k,l):=bg(l,k):=ag(0,0)*(ag(1,1)*ag(2,2)-ag(1,2)*ag(2,1))
	    -ag(0,1)*(ag(1,0)*ag(2,2)-ag(1,2)*ag(2,0))
	    +ag(0,2)*(ag(1,0)*ag(2,1)-ag(1,1)*ag(2,0));
    end;
    detg:=gdd(0,0)*bg(0,0)-gdd(0,1)*bg(0,1)+gdd(0,2)*bg(0,2)-gdd(0,3)*bg(0,3);

    for k:=0:3 do for l:=k:3 do guu(k,l):=guu(l,k):=(-1)**(k+l)*bg(l,k)/detg;
    clear ag,bg;
    write "GUU and detg computed.";
end;

procedure comguu; begin integer i1, i2;
    matrix gddm(spdim+1, spdim+1), guum(spdim+1,spdim+1);
    for i1:=0:spdim do for i2:=0:spdim do gddm(i1+1,i2+1) := gdd(i1,i2);
    detg := det(gddm);
    guum := gddm^(-1);
    array guu(spdim,spdim);
    for i1:=0:spdim do for i2:=0:spdim do guu(i1,i2) := guum(i1+1,i2+1);
    clear(gddm, guum);
    write "GUU and detg computed."
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flat metric %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure filleta; begin integer k;
    array eta(spdim,spdim);
    eta(0,0):=1;
    for k:=0:spdim do eta(k,k) := -1;
end;

procedure setgflat; begin integer k;
    array gdd(spdim,spdim), guu(spdim,spdim);
    guu(0,0) := gdd(0,0) := 1;
    for k:=1:spdim do guu(k,k) := gdd(k,k) := -1;
end;

% Spherical coordinates in spdim dimensions
% Assumptions: x(0) = t
%              x(1) = r
%              x(2) = theta1, 0..pi
%              ...
%              x(spdim-1) = theta(spdim-2). 0..pi
%              x(spdim)   = phi
procedure setgpolar; begin integer k,l;
    array gdd(spdim,spdim);
    gdd(0,0) := 1;
    gdd(1,1) := -1;
    for k:= 2:spdim do gdd(k,k) := -r^2*for l:=2:k-1 product sin(x(l))^2;
end;

%%%%%%%%%%%%%% Levi-Civita symbol in 3 and four dimensions %%%%%%%%%%%%%%
procedure filleps3; begin
    array eps3(3,3,3);
    eps3(1,2,3) := eps3(2,3,1) := eps3(3,1,2) := 1;
    eps3(3,2,1) := eps3(2,1,3) := eps3(1,3,2) := -1;
end;

procedure filleps4; begin
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
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Volume form %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rdetg must be calculated first; not here due to sign
procedure comvolform3; begin integer i1,i2,i3, j1, j2, j3;
    array E3ddd(2,2,2), E3uuu(2,2,2);

    E3uuu(0,1,2) := E3uuu(1,2,0) := E3uuu(2,0,1) := 1/rdetg;
    E3uuu(2,1,0) := E3uuu(1,0,2) := E3uuu(0,2,1) := -1/rdetg;

    E3ddd(0,1,2) := E3ddd(1,2,0) := E3ddd(2,0,1) := rdetg;
    E3ddd(2,1,0) := E3ddd(1,0,2) := E3ddd(0,2,1) := -rdetg;
end;

procedure comvolform4; begin integer i1,i2,i3,i4;
    array E4dddd(3,3,3,3), E4uuuu(3,3,3,3);

    for i1:=0:3 do for i2:=0:3 do for i3:=0:3 do for i4:=0:3 do E4dddd(i1,i2,i3,i4) := eps4dddd(i1,i2,i3,i4)*rdetg;

    for i1:=0:3 do for i2:=0:3 do for i3:=0:3 do for i4:=0:3 do E4uuuu(i1,i2,i3,i4) := eps4uuuu(i1,i2,i3,i4)/rdetg;
end;

%%%%%%%%%%%%%%%%%%%%%%%%% Christoffel symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates Levi-Civita connection (zero torsion);
procedure comchristoffel; begin integer k,l,m,n;
    array gaddd(spdim,spdim,spdim), gaudd(spdim,spdim,spdim);

    for k:=0:spdim do for l:=0:spdim do for m:=0:spdim do gaddd(k,l,m) := 1/2 * ( df(gdd(k,l),x(m)) + df(gdd(k,m),x(l)) - df(gdd(l,m),x(k)) );

    for k:=0:spdim do for l:=0:spdim do for m:=0:spdim do gaudd(k,l,m) := for n:=0:spdim sum guu(k,n)*gaddd(n,l,m) ;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Riemann tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure comriemann; begin integer k,l,m,n,s;
    array Rieuddd(spdim,spdim,spdim,spdim), Riedddd(spdim,spdim,spdim,spdim), Rieuudd(spdim,spdim,spdim,spdim), Rieuuuu(spdim,spdim,spdim,spdim);

    % LL2 92.1, faster according to MacCallum (no guu derivatives)
    for j:=0:spdim do for k:=0:spdim do for l:=0:spdim do for m:=0:spdim do Riedddd(j,k,l,m) := 
	1/2*( df(gdd(j,m),x(k),x(l)) +df(gdd(k,l),x(j),x(m)) -df(gdd(j,l),x(k),x(m)) -df(gdd(k,m),x(j),x(l)) )
	+ for n:=0:spdim sum for s:=0:spdim sum gdd(n,s)*( gaudd(n,k,l)*gaudd(s,j,m) - gaudd(n,k,m)*gaudd(s,j,l) );

    for i1:=0:spdim do for i2:=0:spdim do for i3:=0:spdim do for i4:=0:spdim do Rieuddd(i1,i2,i3,i4) := 
	for j1:=0:spdim sum guu(i1,j1)*Riedddd(j1,i2,i3,i4);

    for i1:=0:spdim do for i2:=0:spdim do for i3:=0:spdim do for i4:=0:spdim do Rieuudd(i1,i2,i3,i4) :=
	for j2:=0:spdim sum Rieuddd(i1,j2,i3,i4)*guu(i2,j2);

    for i1:=0:spdim do for i2:=0:spdim do for i3:=0:spdim do for i4:=0:spdim do Rieuuuu(i1,i2,i3,i4) :=
	for j3:=0:spdim sum for j4:=0:spdim sum Rieuudd(i1,i2,j3,j4)*guu(i3,j3)*guu(i4,j4);

    write "Riemann tensor computed.";
end;

% verify symmetry
procedure riemsymm; begin integer k,l,m,n;
    for k:=0:spdim do for l:=0:spdim do for m:=0:spdim do for n:=0:spdim do begin scalar qq;
	qq := Riedddd(k,l,m,n)+Riedddd(l,k,m,n);
	if qq neq 0 then write "Error in Riemann tensor antisymmetry, first pair, indices: ", k,l,m,n;
    end;

    for k:=0:spdim do for l:=0:spdim do for m:=0:spdim do for n:=0:spdim do begin scalar qq;
	qq := Riedddd(k,l,m,n)+Riedddd(k,l,n,m);
	if qq neq 0 then write "Error in Riemann tensor antisymmetry, first pair, indices: ", k,l,m,n;
    end;

    for k:=0:spdim do for l:=0:spdim do for m:=0:spdim do for n:=0:spdim do begin scalar qq;
	qq := Riedddd(k,l,m,n)-Riedddd(m,n,k,l);
	if qq neq 0 then write "Error in Riemann tensor antisymmetry, first pair, indices: ", k,l,m,n;
    end;

    write "Symmetry properties of the Riemann tensor verified.";
end;

% verify cyclic
procedure bianchialg; begin scalar k,l,m,n;
    for k:=0:spdim do for l:=0:spdim do for m:=0:spdim do for n:=0:spdim do begin scalar qq;
	qq := Riedddd(k,l,m,n)+Riedddd(k,n,l,m)+Riedddd(k,m,n,l);
	if qq neq 0 then write "Error in Riemann tensor cyclic property, first pair, indices: ", k,l,m,n;
    end;

    write "Algebraic Bianchi (Ricci) identity verified.";
end;

%%%%%%%%%%%%%%%%%%%%%%%% Verify Bianchi identity %%%%%%%%%%%%%%%%%%%%%%%%
procedure bianchidiff; begin integer n,j,k,l,m,n,p;
    array DRieudddd(spdim,spdim,spdim,spdim,spdim);

    for n:=0:spdim do for j:=0:spdim do for k:=0:spdim do for l:=0:spdim do for m:=0:spdim do DRieudddd(n,j,k,l,m) := 
	df( Rieuddd(n,j,k,l), x(m))
	+ for p:=0:spdim sum ( gaudd(n,p,m)*Rieuddd(p,j,k,l) - gaudd(p,j,m)*Rieuddd(n,p,k,l) - gaudd(p,k,m)*Rieuddd(n,j,p,l)
						    - gaudd(p,l,m)*Rieuddd(n,j,k,p) ) ;

    for n:=0:spdim do for j:=0:spdim do for k:=0:spdim do for l:=0:spdim do for m:=0:spdim do begin scalar qq;
	qq := DRieudddd(n,j,k,l,m) + DRieudddd(n,j,m,k,l) + DRieudddd(n,j,l,m,k);
	if qq neq 0 then write "Error in Bianchi identity, indices: ", n,j,k,l,m;
    end;

    clear(DRieudddd);
    write "Differential Bianchi identity verified.";
end;

%%%%%%%%%%%%%%%%%%%%%%%%% Ricci tensor and scalar %%%%%%%%%%%%%%%%%%%%%%%
procedure comricci; begin integer k,l,m;
    array Ricdd(spdim,spdim);

    for k:=0:spdim do for l:=0:spdim do Ricdd(k,l) := for m:=0:spdim sum Rieuddd(m,k,m,l) ;

    Rs := for k:=0:spdim sum for l:=0:spdim sum guu(k,l)*Ricdd(l,k) ;
    write "Ricci tensor and scalar computed.";
end;

% verify formula, 92.7
    procedure verifyricci; begin integer j,k,l,m;
    array Ricddp(spdim,spdim);

    for j:=0:spdim do for k:=0:spdim do Ricddp(j,k) := for l:=0:spdim sum ( df(gaudd(l,j,k),x(l)) - df(gaudd(l,j,l),x(k))
		+ for m:=0:spdim sum ( gaudd(l,j,k)*gaudd(m,l,m) - gaudd(m,j,l)*gaudd(l,k,m) ) ) ;

    for k:=0:spdim do for l:=0:spdim do begin scalar qq;
	qq := Ricdd(k,l) - Ricddp(k,l);
	if qq neq 0 then write "Error in Ricci tensor, indices: ", k, l;
    end;
    write "Formula LL2 92.7 for the Ricci tensor verified."
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Einstein tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure comeinstein; begin scalar j,k,l;
    array Eidd(spdim,spdim), Eiud(spdim,spdim);

    for k:=0:spdim do for l:=0:spdim do Eidd(k,l) := Ricdd(k,l) - 1/2 * gdd(k,l)*Rs;

    for j:=0:spdim do for k:=0:spdim do Eiud(j,k) := for l:=0:spdim sum guu(j,l)*Eidd(l,k);

    write "Einstein tensor computed.";
end;

% verify vacuum solution
procedure verifyvac; begin integer k,l;
    for k:=0:spdim do for l:=0:spdim do if Eidd(k,l) neq 0 then write "Error in Einstein equations, indices: ", k, l;
    write "Vacuum solution (vanishing of the Einstein tensor) verified.";
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Weyl tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
procedure comweyl; begin integer j,k,l,m;
    array Wedddd(spdim,spdim,spdim,spdim);

    for j:=0:spdim do for k:=0:spdim do for l:=0:spdim do for m:=0:spdim do Wedddd(j,k,l,m) := 
	Riedddd(j,k,l,m)
	- 1/2 * Ricdd(j,l)*gdd(k,m) + 1/2 * Ricdd(j,m)*gdd(k,l) + 1/2*Ricdd(k,l)*gdd(j,l)
	+ 1/6 * Rs * ( gdd(j,l)*gdd(k,m) - gdd(j,m)*gdd(k,l) );

    write "Weyl tensor computed.";
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Invariants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dual of the Riemann tensor
procedure comdualriemann; begin integer i1,i2,i3,i4, j1,j2,j3,j4;
    array rsdddd(spdim,spdim,spdim,spdim), rsuuuu(spdim,spdim,spdim,spdim), rsdduu(spdim,spdim,spdim,spdim);

    for i1:=0:spdim do for i2:=0:spdim do for i3:=0:spdim do for i4:=0:spdim do rsdddd(i1,i2,i3,i4) :=
	for j1:=0:spdim sum for j2:=0:spdim sum 1/2*E4dddd(i1,i2,j1,j2)*Rieuudd(j1,j2,i3,i4);

    for i1:=0:spdim do for i2:=0:spdim do for i3:=0:spdim do for i4:=0:spdim do rsuuuu(i1,i2,i3,i4) :=
	for j1:=0:spdim sum for j2:=0:spdim sum for j3:=0:spdim sum for j4:=0:spdim sum
	    rsdddd(j1,j2,j3,j4)*guu(i1,j1)*guu(i2,j2)*guu(i3,j3)*guu(i4,j4);

    for i1:=0:spdim do for i2:=0:spdim do for i3:=0:spdim do for i4:=0:spdim do rsdduu(i1,i2,i3,i4) :=
	for j3:=0:spdim sum for j4:=0:spdim sum rsdddd(i1,i2,j3,j4)*guu(i3,j3)*guu(i4,j4);

    write "Hodge dual of the Riemann tensor computed."
end;

procedure cominvars; begin integer j,k,l,m,p,s;

    Invar1 := for j:=0:spdim sum for k:=0:spdim sum for l:=0:spdim sum for m:=0:spdim sum
	1/48*( Riedddd(j,k,l,m) * Rieuuuu(j,k,l,m) - i * Riedddd(j,k,l,m)*rsuuuu(j,k,l,m));

    Invar2 := for j:=0:spdim sum for k:=0:spdim sum for l:=0:spdim sum for m:=0:spdim sum for p:=0:spdim sum for s:=0:spdim sum
	1/96*( Riedddd(j,k,l,m)*Rieuuuu(l,m,p,s)*Rieuudd(j,k,p,s) + i * Riedddd(j,k,l,m)*Rieuuuu(l,m,p,s)*rsdduu(p,s,j,k) );

    write "Invariants I1 and I2 computed."
end;

write "GR specific subroutines loaded. Warning! do not forget to set the variable spdim to the number of spatial dimensions!";

symbolic onoff('echo, oldecho!*)$

end;
