function [C]= CurvePoint(n,p,U,Pw1,u)

span=FindSpan(n,p,u,U);
N = BasisFuns(span,u,p,U);
Cw=0;
for j=1:p+1
    Cw=Cw+N(j)*Pw1(span-p+j-1,:);
end
C=Cw/Cw(1,4);
end
    
