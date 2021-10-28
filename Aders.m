function [aders]=Aders(n,p,U,Pw2,u,d)

du = min(d,p);
for k=p+1+1:d+1
    aders(k,1:3)=0;
end
span = FindSpan(n,p,u,U);
nders = DersBasisFuns(span,u,p,du,U);
for k=1:du+1
    aders(k,1:3) = 0;
    for j=1:p+1
        aders(k,1:3) = aders(k,1:3) + nders(k,j)*Pw2(span-p+j-1,:);
    end
end
end
    