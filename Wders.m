function [wders]=Wders(n,p,U,u,d,weight)

du = min(d,p);
for k=p+1:d+1
    wders(k,1)=0;
end
span = FindSpan(n,p,u,U);
nders = DersBasisFuns(span,u,p,du,U);
for k=1:du+1
    wders(k,1) = 0;
    for j=1:p+1
        wders(k,1) = wders(k,1) + nders(k,j)*weight(span-p+j-1);
    end
end
end
    