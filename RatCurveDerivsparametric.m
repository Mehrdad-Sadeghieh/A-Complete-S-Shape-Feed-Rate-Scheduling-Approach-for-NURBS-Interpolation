function CKparametric = RatCurveDerivsparametric(d,p,U,Pw2,u,weight,spannumber)

% CK=zeros(d+1,3);
for k=1:d+1
aders = Adersparametric(p,U,Pw2,u,(k-1),spannumber);
v=aders(k,:);
for i=2:k
wders=Wdersparametric(p,U,u,(i-1),weight,spannumber);
v = v - nchoosek(k-1,i-1) *wders(i)*CKparametric(k-i+1,1:3);
end
CKparametric( k,1:3)= v/Wdersparametric(p,U,u,0,weight,spannumber);
end
end
