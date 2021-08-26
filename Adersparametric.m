function [adersparametric]=Adersparametric(p,U,Pw2,u,d,spannumber)

syms u
adersparametric(1,1:3)=u;
adersparametric(1,1:3)=adersparametric(1,1:3)-u;
du = min(d,p);
for k=p+1+1:d+1
    adersparametric(k,1:3)=0;
end
span = spannumber;
nders = DersBasisFunsparametric(span,u,p,du,U);
for k=1:du+1
    adersparametric(k,1:3) = 0;
    for j=1:p+1
        adersparametric(k,1:3) = adersparametric(k,1:3) + nders(k,j)*Pw2(span-p+j-1,:);
    end
end
end