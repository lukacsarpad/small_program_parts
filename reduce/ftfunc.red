%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field theory related subroutines                                     %
%                                                                      %
%                                                                      %
% Preliminaries:                                                       %
%   * Set spdim to number of spatial dimensions                        %
%   * Set up operator x, x(i)=coordinates, i=0:spdim                   %
%   * Set array guu(spdim,spdim), gdd(spdim,spdim) to metric           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global '(oldecho!*) $
symbolic << oldecho!* := !*echo; off1 'echo >>;  % or use on1 to turn echo on.


%%%%%%%%%%%%%%%%%%%%%%%%%%% Electrodynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% index 1: U(1)
% must have A1d set, and gdd, guu.
% coupling constant: ggp (g=coupling constant, gp = gprime = g')
% if only electrodynamics, set ggp = 2*e.
procedure comfieldstrength1; begin integer j,k, l, m;
    array F1dd(spdim,spdim), F1uu(spdim,spdim);

    for j:=0:spdim do for k:=0:spdim do F1dd(j,k) := df( A1d(k), x(j) ) - df( A1d(j), x(k) );
    for j:=0:spdim do for k:=0:spdim do F1uu(j,k) := for l:=0:spdim sum for m:=0:spdim sum guu(j,l)*guu(k,m)*F1dd(l,m);

    write "Maxwell tensor F1dd, Fuu computed.";
end;

procedure Dc1(scalr,indk); df( scalr, x(indk) ) - i*ggp*A1d(indk)*scalr/2;

procedure Lag1; begin integer i1,i2;
    return -1/4*( for i1:=0:spdim sum for i2:=0:spdim sum F1uu(i1,i2)*F1dd(i1,i2) );
end;

% stress-energy tensor
procedure comT1; begin scalar lagr1; integer i1,i2,i3,i4;
    array T1dd(spdim,spdim), T1uu(spdim,spdim), T1ud(spdim,spdim);

    lagr1 := Lag1();

    for i1:=0:spdim do for i2:=0:spdim do T1dd(i1,i2) := ( for i3:=0:spdim sum for i4:=0:spdim sum
	guu(i3,i4) * F1dd(i1,i3)*F1dd(i4,i2) ) - gdd(i1,i2)*lagr1;

    for i1:=0:spdim do for i2:=0:spdim do T1ud(i1,i2) := ( for i3:=0:spdim sum
	F1uu(i1,i3)*F1dd(i3,i2) ) - (if i1 = i2 then 1 else 0)*lagr1;

    for i1:=0:spdim do for i2:=0:spdim do T1uu(i1,i2) := ( for i3:=0:spdim sum for i4:=0:spdim sum
	gdd(i3,i4) * F1uu(i1,i3)*F1uu(i4,i2) ) - guu(i1,i2)*lagr1;

    write "Maxwell stress energy tensor T1dd, T1ud, T1uu computed.";

end;

%%%%%%%%%%%%%%%%%%%%%%%%%% SU(2) field theory %%%%%%%%%%%%%%%%%%%%%%%%%%
% index 2: SU(2)
% must have A2d set, and gdd, guu.
% coupling constant: gg (g=coupling constant, g)
% generators: Pauli matrices
procedure settau; begin
    array tau(3,2,2);
    tau(1,1,2) := tau(1,2,1) := 1;
    tau(2,1,2) := -i;
    tau(2,2,1) := i;
    tau(3,1,1) := 1;
    tau(3,2,2) := -1;
    tau(0,1,1) := tau(0,2,2) := 1;
end;

% commutators
procedure comcomm2; begin integer a1,a2,a3,b1,b2,b3;
    array c2(3,3,3);
    for a1:=1:3 do for a2:=1:3 do for a3:=1:3 do c2(a1,a2,a3) := 1/4/i * ( for b1:=1:2 sum for b2:=1:2 sum for b3:=1:2 sum
	tau(a1,b1,b2) * ( tau(a2,b2,b3) * tau(a3,b3,b1) - tau(a3,b2,b3) * tau(a2,b3,b1) ) );
end;

procedure comfieldstrength2; begin integer a1,a2,a3,j,k,l,m;
    array F2dd(3,spdim,spdim), F2uu(3,spdim, spdim);

    for a1:=1:3 do for j:=0:spdim do for k:=0:spdim do F2dd(a1,j,k) := df( A2d(a1,k), x(j) ) - df( A2d(a1,j), x(k) ) +
	for a2:=1:3 sum for a3:=1:3 sum gg * c2(a1,a2,a3) * A2d(a2,j) * A2d(a3,k);

    for a1:=1:3 do for j:=0:spdim do for k:=0:spdim do F2uu(a1,j,k) := for l:=0:spdim sum for m:=0:spdim sum guu(j,l)*guu(k,m)*F2dd(a1,l,m);

    write "SU(2) field strength tensor F2dd, F2uu computed.";
end;


procedure Dc2(scalr,inda,indk); begin integer a1,a2;
    return df( scalr(inda), x(indk) ) - i * gg / 2 * ( for a1:=1:3 sum for a2:=1:2 sum
	tau(a1,inda,a2) * A2d(a1,indk) * scalr(a2) );
end;

procedure Dc12(scalr,inda,indk); begin integer a1,a2;
    return df( scalr(inda), x(indk) ) - i * gg / 2 * ( for a1:=1:3 sum for a2:=1:2 sum
	tau(a1,inda,a2) * A2d(a1,indk) * scalr(a2) ) -i * ggp / 2 * A1d(indk) * scalr(inda);
end;

% adjoint representation
procedure Dca2(scalr,inda,indk); begin integer a1,a2;
    return df( scalr(inda), x(indk) ) + gg * ( for a1:=1:3 sum for a2:=1:3 sum
	c2(inda,a1,a2) * A2d(a1,indk) * scalr(a2) );
end;

procedure Lag2; begin integer i1,i2;
    return -1/4 * ( for i1:=0:spdim sum for i2:=0:spdim sum for a1:=1:3 sum F2uu(a1,i1,i2)*F2dd(a1,i1,i2) );
end;

procedure comT2; begin scalar lagr2; integer i1,i2,i3,i4, a1;
    array T2dd(spdim,spdim), T2uu(spdim,spdim), T2ud(spdim,spdim);

    lagr2 := Lag2();

    for i1:=0:spdim do for i2:=0:spdim do T2dd(i1,i2) := ( for i3:=0:spdim sum for i4:=0:spdim sum for a1:=1:3 sum
	guu(i3,i4) * F2dd(a1,i1,i3)*F2dd(a1,i4,i2) ) - gdd(i1,i2)*lagr2;

    for i1:=0:spdim do for i2:=0:spdim do T2ud(i1,i2) := ( for i3:=0:spdim sum for a1:=1:3 sum
	F2uu(a1,i1,i3)*F2dd(a1,i3,i2) ) - (if i1 = i2 then 1 else 0)*lagr2;

    for i1:=0:spdim do for i2:=0:spdim do T2uu(i1,i2) := ( for i3:=0:spdim sum for i4:=0:spdim sum for a1:=1:3 sum
	gdd(i3,i4) * F2uu(a1,i1,i3)*F2uu(a1,i4,i2) ) - guu(i1,i2)*lagr2;

    write "Stress energy tensor T2dd, T2ud, T2uu of SU(2) gauge field computed.";

end;


%%%%%%%%%%%%%%%%%%%%%%%%%% SU(3) field theory %%%%%%%%%%%%%%%%%%%%%%%%%%
% index 3: SU(3)
% must have A3d set, and gdd, guu.
% coupling constant: g3 (g=coupling constant, g)
% generators: Gell-Mann matrices
procedure setgellmann; begin
    array gml(8,3,3);
    gml(1,1,2) := gml(1,2,1) := 1;
    gml(2,1,2) := -i;
    gml(2,2,1) := i;
    gml(3,1,1) := 1;
    gml(3,2,2) := -1;
    gml(4,1,3) := gml(4,3,1) := 1;
    gml(5,1,3) :=-i;
    gml(5,3,1) := i;
    gml(6,2,3) := gml(6,3,2) := 1;
    gml(7,2,3) := -i;
    gml(7,3,2) := i;
    gml(8,1,1) := 1/sqrt(3);
    gml(8,2,2) := 1/sqrt(3);
    gml(8,3,3) := -2/sqrt(3);
end;

% commutators
procedure comcomm3; begin integer a1,a2,a3,b1,b2,b3;
    array c3(8,8,8);
    for a1:=1:8 do for a2:=1:8 do for a3:=1:8 do c3(a1,a2,a3) := 1/4/i * ( for b1:=1:3 sum for b2:=1:3 sum for b3:=1:3 sum
	gml(a1,b1,b2) * ( gml(a2,b2,b3) * gml(a3,b3,b1) - gml(a3,b2,b3) * gml(a2,b3,b1) ) );
end;

procedure comacomm3; begin integer a1,a2,a3,b1,b2,b3;
    array ac3(8,8,8);
    for a1:=1:8 do for a2:=1:8 do for a3:=1:8 do ac3(a1,a2,a3) := 1/4 * ( for b1:=1:3 sum for b2:=1:3 sum for b3:=1:3 sum
	gml(a1,b1,b2) * ( gml(a2,b2,b3) * gml(a3,b3,b1) + gml(a3,b2,b3) * gml(a2,b3,b1) ) ) - 2/3 * ( if a2 = a3 then 1 else 0 );
end;

procedure comfieldstrenth3; begin integer a1,a2,a3,j,k,l,m;
    array F2dd(8,spdim,spdim), F2uu(8,spdim, spdim);

    for a1:=1:8 do for j:=0:spdim do for k:=0:spdim do F3dd(a1,j,k) := df( A3d(a1,k), x(j) ) - df( A3d(a1,j), x(k) ) +
	for a2:=1:8 sum for a3:=1:8 sum g3 * c3(a1,a2,a3) * A3d(a2,j) * A3d(a3,k);

    for a1:=1:8 do for j:=0:spdim do for k:=0:spdim do F3uu(a1,j,k) := for l:=0:spdim sum for m:=0:spdim sum guu(j,l)*guu(k,m)*F3dd(a1,l,m);

    write "SU(3) field strength tensor F3dd, F3uu computed.";
end;

procedure Dc3(scalr,inda,indk); begin integer a1,a2;
    return df( scalr(inda), x(indk) ) - i * g3 / 2 * ( for a1:=1:8 sum for a2:=1:3 sum
	glm(a1,indaa,a2) * A3d(a1,indk) * scalr(a2) );
end;

procedure Dc123(scalr,inda,indb,indk); begin integer a1,a2;
    return df( scalr(inda,indbb), x(indk) ) - i * g3 / 2 * ( for a1:=1:8 sum for a2:=1:3 sum
	gml(a1,inda,a2) * A3d(a1,indk) * scalr(a2,indb) ) - i * gg / 2 * ( for a1:=1:3 sum for a2:=1:2 sum
	tau(a1,indb,a2) * A2d(a1,indk) * scalr(inda,a2) ) - i * ggp / 2 * A1d(indk) * scalr(inda,indb);
end;


procedure Lag3; begin integer i1,i2;
    return -1/4 * ( for i1:=0:spdim sum for i2:=0:spdim sum for a1:=1:8 sum F3uu(a1,i1,i2)*F3dd(a1,i1,i2) );
end;

procedure comT3; begin scalar lagr3; integer i1,i2,i3,i4, a1;
    array T3dd(spdim,spdim), T3uu(spdim,spdim), T3ud(spdim,spdim);

    lagr3 := Lag3();

    for i1:=0:spdim do for i2:=0:spdim do T3dd(i1,i2) := ( for i3:=0:spdim sum for i4:=0:spdim sum for a1:=1:8 sum
	guu(i3,i4) * F3dd(a1,i1,i3)*F3dd(a1,i4,i2) ) - gdd(i1,i2)*lagr3;

    for i1:=0:spdim do for i2:=0:spdim do T3ud(i1,i2) := ( for i3:=0:spdim sum for a1:=1:8 sum
	F3uu(a1,i1,i3)*F3dd(a1,i3,i2) ) - (if i1 = i2 then 1 else 0)*lagr3;

    for i1:=0:spdim do for i2:=0:spdim do T3uu(i1,i2) := ( for i3:=0:spdim sum for i4:=0:spdim sum for a1:=1:8 sum
	gdd(i3,i4) * F3uu(a1,i1,i3)*F3uu(a1,i4,i2) ) - guu(i1,i2)*lagr3;

    write "Stress energy tensor T3dd, T3ud, T3uu of SU(3) gauge field computed.";

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dirac matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%

procedure setdirmatchir; begin;
    array dgamm(3);
    dgamm(0) := mat((0,0,1,0),(0,0,0,1),(1,0,0,0),(0,1,0,0));
    dgamm(1) := mat((0,0,0,1),(0,0,1,0),(0,-1,0,0),(-1,0,0,0));
    dgamm(2) := mat((0,0,0,-i),(0,0,i,0),(0,i,0,0),(-i,0,0,0));
    dgamm(3) := mat((0,0,1,0),(0,0,0,-1),(-1,0,0,0),(0,1,0,0));
    dgamm5 := i * dgamm(0) * dgamm(1) * dgamm(2) * dgamm(3);
end;

procedure setdirarrchir; begin;
    array dgama(3,4,4);
    dgama(0,3,1) := dgama(0,2,3) := 1;
    dgama(0,3,1) := dgama(0,4,2) := 1;

    dgama(1,1,4) := dgama(1,2,3) := 1;
    dgama(1,3,2) := dgama(1,4,1) := -1;

    dgama(2,1,4) := dgama(2,4,1) := -i;
    dgama(2,2,3) := dgama(2,3,2) := i;

    dgama(3,1,3) := dgama(3,4,2) := 1;
    dgama(3,2,4) := dgama(3,3,1) := -1;
end;

write "Field theory related subroutines loaded. Do not forget to have spdim set to number of spatial dimensions!";

symbolic onoff('echo, oldecho!*)$

end;
