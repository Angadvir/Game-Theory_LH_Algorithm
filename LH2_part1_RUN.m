function [y,z]=LH2_part1_RUN()  % returns Nash equilibrium found by the Lemke-Howson algorithm starting with (0,0) and label k, 1 <= k <= m+n
global epsilon M b m n flag en;  % global variables
clc; %counter =0;
    en =int16(8*rand(1)); 
    L = rand(en,en)-0.85*rand(en,en);
    [m n] = size(L);
    %Let's first write all the linear constraints {Ay <=1; x^TB <= 1; x>=0; y>=0 as Mz <= b  where z=(x y)
    Mp = [L -1.*ones(m,1)];
    M = [Mp; -eye(m,n) zeros(m,1); zeros(1,n) -1]
    bb = rand(m,1) - 0.85*rand(m,1);
    b = [bb; zeros(m+1,1)]
    epsilon = 10^(-6); % For floating pt. computation
    eps = 10^(-6); 
    z_ini = 0; 
    [x,y] = size(b);
    for i =1:1:x
        maxi = -b(i);
        if maxi >= z_ini
            z_ini = maxi;
            label = i;
        end
    end
    iimin = label;
    v = [zeros(1,m) z_ini]; % starting vertex
    %iimin%it gives the next vertex and the new label acquired when moving from v and leaving hyperplane corresponding to k
    count=0;  % it will store the number of pivotings
    if(b)>=0
        fprintf(1, 'Trivial Solution Found!\n');
        y=v(1:m);
        return;
    end
    while abs(v(m+1))> epsilon
        count=count+1;
        if iimin<=m
           [v,iimin]=computenextvertex(v,iimin+m);
        else
            [v,iimin]=computenextvertex(v,iimin-m);
        end;
        if flag==1
           fprintf(1, 'Secondary Ray Encountered\n');
         return;
        end
    end;
    y=v(1:m);
    z=v(m+1);
    return;	

% computenextvertex takes a vertex v and a tight inequality ineq at v,
% and outputs another vertex v after leaving ineq at v and iimin - which
% captures the newly tight inequality at the new v.
function [v,iimin,flag]=computenextvertex(v,ineq);
global M epsilon b m n flag totC counter;
flag = 0;
[I,s1]=computeI(v,ineq);
[k,l]=size(I);
iimin=0;
if l > m+1
	fprintf(1, 'Degeneracy\n');	
	return; 
elseif l < m+1
    fprintf(1, 'Not a vertex\n');
	return;
end;
dirs = -inv(M(I,:));
dd=dirs(:,s1);
I1=I;
betamin=10^10;
[k,l]=size(M);
for i=1:1:k 
	[h,ii]=ismember(i,I1);
	aa=M(i,:);
	alpha=aa*dd;
	if h==0 && alpha>0
		beta=(b(i)-aa*v')/(aa*dd);
		if beta<betamin
			betamin=beta;
			iimin=i;
		end;
	end;
end;
v = v+betamin*dd';
if betamin == 10^10
    flag = 1;
end
return;

% computeI takes a vertex v and the tight inequality ineq at v, and outputs
% the set I of tight inequalities at v and the position of ineq in I.
function [I,s]=computeI(v,ineq);
global M b epsilon;
[k l] = size(M);
I=[];
s=0;
for i=1:1:k
    %M(i,:)*v'
	if abs(M(i,:)*v'- b(i)) < epsilon
		I = [I i];
	end;
	if i==ineq
		[dummy,s]=size(I);
	end;
end;
return;
