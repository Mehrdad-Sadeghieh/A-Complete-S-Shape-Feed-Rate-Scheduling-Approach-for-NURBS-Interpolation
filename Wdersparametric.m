function [wdersparametric]=Wdersparametric(p,U,u,d,weight,spannumber)

syms u
wdersparametric(1,1)=u;
wdersparametric(1,1)=wdersparametric(1,1)-u;

du = min(d,p);
for k=p+1:d+1
    wdersparametric(k,1)=0;
end
span = spannumber;
nders = DersBasisFunsparametric(span,u,p,du,U);
for k=1:du+1
    wdersparametric(k,1) = 0;
    for j=1:p+1
        wdersparametric(k,1) = wdersparametric(k,1) + nders(k,j)*weight(span-p+j-1);
    end
end

end