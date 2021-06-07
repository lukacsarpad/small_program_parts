% Code by Z. Perjés
% moved into blocks by Á. Lukács
off echo;

%%% makes the polynomial f homogeneous in sin(x) and cos(x)   %%%
%%% clear substitutions like sin(x)^2:=1-cos(x)^2 before using it
 
procedure trn(f,x); begin scalar ftrn, dftrn, ftrn1, ftrn2;
    if den(f)=1 then <<
	cos(x):=atrn*ctrn;
	sin(x):=atrn*strn;
	ftrn:=f;
	dftrn:=deg(ftrn,atrn);
% write "degree=",dftrn;
	ftrn1:=for itrn:=dftrn step -2 until 0 sum
		    coeffn(ftrn,atrn,itrn)*atrn^itrn*(atrn^2*(strn^2+ctrn^2))^((dftrn-itrn)/2);
	ftrn2:=for itrn:=dftrn-1 step -2 until 0 sum
		    coeffn(ftrn,atrn,itrn)*atrn^itrn*(atrn^2*(strn^2+ctrn^2))^((dftrn-1-itrn)/2);
	ftrn:=ftrn1+ftrn2;
% if ftrn2=0 then
%              if evenp(dftrn) then write "only even degree terms"
%                              else write "only odd degree terms"
%            else write "both even and odd order terms";
	clear sin(x),cos(x);
	ctrn:=cos(x)/atrn;
	strn:=sin(x)/atrn;
	ftrn:=ftrn;
	clear ctrn,strn;
	return ftrn;
    >>
    else write "not a polynomial";
end;
 
%%% simplifies a rational expression,
%%%  making the numerator and the denominator homogeneous
%%%  and dividing off with common factors
 
procedure trnr(f,x); begin scalar nftrnr, dftrnr, gftrnr, ftrnr;
    nftrnr:=trn(num(f),x);
    dftrnr:=trn(den(f),x);
    gftrnr:=gcd(nftrnr,dftrnr);
    nftrnr:=nftrnr/gftrnr;
    dftrnr:=dftrnr/gftrnr;
    ftrnr:=nftrnr/dftrnr;
    if trn(num(f-ftrnr),x)=0 then return ftrnr
                        else write "ERROR";
end;

on echo;
end;
