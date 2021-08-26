function [U,P]=WCLeastSquaresCurve(Q,r,Wq,D,s,I,Wd,n,p)
%r starts from 1
%s starts from 1
%n starts from 0
ru=0;
rc=0;
for i=1:r
    if Wq(i)>0
        ru=ru+1;
    else
        rc=rc+1;
    end
end
su=0;
sc=0;

for j=1:s
    if Wd(j)>0
        su=su+1;
    else 
        sc=sc+1;
    end
end

mu=ru+su+1;
mc=rc+sc+1;

if mc>=n+2 || mc+n>=mu+1
    disp('error')
end

d=0;
for i=2:size(Q,1)
    d=d+norm(Q(i,:)-Q(i-1,:));
end

ub(1)=0;
ub(size(Q,1))=1;
for i=2:size(Q,1)-1
    ub(i)=ub(i-1)+(norm(Q(i,:)-Q(i-1,:))/d);
end

dd=r/(n-p+1);
U(1,1:p+1)=0;
for j=1:n-p
    i=floor(j*dd);
    alpha=j*dd-i;
    U(p+j+1)=(1-alpha)*ub(i)+alpha*ub(i+1);
end
U(1,p+j+2:p+j+2+p)=1;

%% set up arrays N,W,S,T,M

j=1;
mu2=1;
mc2=1;

for i=1:r
    span=FindSpan(n+1,p,ub(i),U);
    dflag=0;
    if j<=s
        if i==I(j)
            dflag=1;
        end
    end
    
    if dflag==00
        funs=BasisFuns(span,ub(i),p,U);
    else
        funs=DersBasisFuns(span,ub(i),1,U);
    end
    
    if Wq(i)>0
        W(mu2,mu2)=Wq(i);
        N(mu2,span-p:span)=funs(1,:);
        S(mu2,:)=W(mu2,mu2)*Q(i,:);
        mu2=mu2+1;
    else
        M(mc2,span-p:span)=funs(1,:);
        T(mc2)=Q(i,:);
        mc2=mc2+1;
    end
    
    if dflag==1
        if Wd(j)>0
            W(mu2,mu2)=Wd(j);
            N(mu2,span-p:span)=funs(2,:);
            S(mu2,:)=W(mu2)*D(j,:);
            mu2=mu2+1;
        else
            M(mc2,span-p:span)=funs(2,:);
            T(mc2,:)=D(j,:);
            mc2=mc2+1;
        end
        j=j+1;
    end
end

NTWN=transpose(N)*W*N;
NTWS=transpose(N)*W*S;
if mc<2
    P=inv(NTWN)*NTWS;
else
A=inv(M*(inv(NTWN)*transpose(M))) * (M*(inv(NTWN))*NTWS-T);
P=inv(NTWN)*NTWS - inv(NTWN) * transpose(M) *A;
end

  
        
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
