off echo;

%%% Linear Eq. solver.
%  Z. Perjés

procedure solv(f,h);-coeffn(f,h,0)/coeffn(f,h,1);

%%% Sort the output of factorize
% Á. Lukács, 190208

procedure sortfac(l);
begin scalar i,j,t1,t1a,t2,t2a;
    for i:=1:length(l) do for j:=i+1:length(l) do <<
	t1  := part(l,i);
	t1a := first(t1);
	t2  := part(l,j);
	t2a := first(t2);
	if length(t2a) < length(t1a) then <<
	    l := ( part(l, i) := t2 );
	    l := ( part(l, j) := t1 );
	>>;
    >>;
    return l
end;


on echo;
end;
