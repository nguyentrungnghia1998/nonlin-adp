[~,~,yy]=sim('DATA.slx');
x=yy(:,1:2)';
u=yy(:,3)';
W0=[0;3;3/2];
O=zeros(1000,3);
V=zeros(1000,1);
um=zeros(size(u));
tg=zeros(10001,3);
tr=zeros(10001,1);
intg=zeros(1000,3);
for k=1:16
for i=1:10001
    um(i)=-1/2*[0 sin(x(1,i))]*dsig(x(:,i))'*W0;
    tr(i)=-x(:,i)'*x(:,i)+um(i)^2-2*(u(i)-um(i))'*um(i);
end
for i=1:1000
    V(i,:)=0.0005*(0.5*tr(10*i-9,:)+0.5*tr(10*i+1,:)+sum(tr(10*i-8:10*i,:)));
end
for i=1:1000
    O(i,:)=sig(x(:,10*i+1))'-sig(x(:,10*i-9))';
end
Wnew=O\V;
if norm(Wnew-W0)<1e-5
    break;
end
W0=Wnew;
end