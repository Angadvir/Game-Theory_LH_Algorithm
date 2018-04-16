function [y,z]=LH2_part3__2() % returns Nash equilibrium found by the Lemke-Howson algorithm starting with (0,0) and label k, 1 <= k <= m+n
global epsilon M b m n en count flg;  % global variables
  count = 0; 
  % en = 6;%int16(6*rand(1)); %YOU CAN CHANGE 'EN' = 2 OR 3 OR 4, FOR TESTING THE CODE!
    U = 10.*rand(en,en);
    Mon = 10.*rand(1,en);
    sum_mon = sum(Mon);
    [m n] = size(U);
    %Let's first write all the linear constraints {Ay <=1; x^TB <= 1; x>=0; y>=0 as Mz <= b  where z=(x y)
    Dum1 = zeros(m,(m*m)+1);
    for i=1:m
        W(i,:) = Mon(1,i)/sum_mon;
    end
    
    W=repmat(W,1,n);  %MIGHT NEED TO MODIFY THIS
    
    % 1st equation:
    iye = -1.*eye(m);
    Con = [iye Dum1];
    Con1 = repmat(Con,n,1);
    for in=1:n*n
        x=in/n;
        x=ceil(x);
        if in> n
            H=in-(x-1)*n;
            Up(in,x)=U(x,H);
        else
            Up(in,x)=U(x,in);
        end
    end
    Up;
    Upp1 = [Up Con1];
    
    % 2nd Equation:
    Be = zeros(n,n);
    iye2 = repmat(eye(n,n),1,n);
    neg_iye2= -1.*eye(n,n);
    SecZ = zeros(n,1);

    Upp2 = [Be neg_iye2 iye2 SecZ];
    
    % 3rd Equation:
    Part1=zeros(n,n);
    iden3=-1*eye(n,n);
  %  iden4=repmat(iden3,1,n);
    A=zeros(n,n*n);
    temp=n; tog =0;
    for i =1:n
        if tog==1
            A(:,(temp+1):i*n) = repmat(iye2(:,i),1,n); 
        else
            A(:,i:i*n) = repmat(iye2(:,i),1,n);
            tog = 1;
        end
        temp=i*n;
    end
    neg_ones = -1.*ones(n,1);
    Upp3=[Part1 W -1.*A neg_ones];

    % 4th Equation:
    zer1 = zeros(n*n,2*n);
    iye1 = -1.*eye(n*n,n*n);
    lastzer = zeros(n*n,1);
   
    Down1=[zer1 iye1 lastzer];   

    % 5th Equation:
    zer2 = zeros(n,n);
    iye11 = -1.*(eye(n,n));
    zer11 = repmat(zeros(n,n),1,n);
    lastzer1 = zeros(n,1);
    Down2=[zer2 iye11 zer11 lastzer1];    

    % 6th Equation:
    Down3=[iye11 zer11 zer2 lastzer1]; 
    
    % 7th Equation:
    finzer = zeros(1,(n*n)+(2*n));
    Down4=[finzer -1];
   
    M = [Upp1;Upp2;Upp3;Down1;Down2;Down3;Down4]; %BIG M MATRIX 
    
    for i =1:n
     sumrow_W(i)=sum(W(i,:));
    end
    
    % b vector: 
    b = [ones(n*n,1); ones(n,1); -1.*sumrow_W'; zeros((n*n)+2*n,1); 0]; %PLEASE CHECK THIS AGAIN
    epsilon = 10^(-6); % For floating pt. computation
    
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
    v = [zeros(1,(2*n)+(n*n)) z_ini]; % starting vertex
  
    %iimin%it gives the next vertex and the new label acquired when moving from v and leaving hyperplane corresponding to k
    count=1;  % it will store the number of pivotings
    
    while abs(v((2*n+n*n)+1))> epsilon && flg ~=1
        eta = (2*n+n*n);
        count=count + 1
        if iimin<=eta
           [v,iimin]=computenextvertex1(v,iimin+eta);
        else
            [v,iimin]=computenextvertex1(v,iimin-eta);
        end
    end
    y=v(1:(2*n)+(n*n));
    z=v(1,(2*n)+(n*n)+1);
    return;	

% computenextvertex takes a vertex v and a tight inequality ineq at v,
% and outputs another vertex v after leaving ineq at v and iimin - which
% captures the newly tight inequality at the new v.
function [v,iimin]=computenextvertex1(v,ineq);
global M  b m n; 
[I,s1]=computeI1(v,ineq);
[k,l]=size(I);
iimin=0;
cond = (2*n+n*n+1);
if l > cond
	fprintf(1, 'Degeneracy\n');	
	return; 
elseif l < cond
    fprintf(1, 'Not a vertex\n');
	return;
end;
dirs = -inv(M(I,:));
dd=dirs(:,s1);
I1=I;
betamin=10^10;
[k,l]=size(M);
k=k-1;
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
return;

% computeI takes a vertex v and the tight inequality ineq at v, and outputs
% the set I of tight inequalities at v and the position of ineq in I.
function [I,s]=computeI1(v,ineq);
global M b epsilon;
[k l] = size(M);
I=[];
s=0;
for i=1:1:k
	if abs(M(i,:)*v'- b(i)) < epsilon
		I = [I i];
	end;
	if i==ineq
		[dummy,s]=size(I);
	end;
end;
return;