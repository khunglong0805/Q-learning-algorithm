N=1000;
u=zeros(1,N);
y= zeros(1,N);
C=1;
sigma=zeros(2,N);
omega=zeros(2,N);
z=zeros(5,N);
zbar=zeros(15,N);


%%dùng dynamic của hệ để lấy ra input
x=randn;
for k= 1:N
    y(k)= C*x;
    u(k)=randn;
    x= 0.8*x+u(k);
end    
 
for k=3:N-1
    sigma(:,k)=[u(k-1);u(k-2)];
    omega(:,k)=[y(k-1);y(k-2)];
    z(:,k)=[sigma(:,k); omega(:,k); u(k)];
end   

%%chuyển từ z sang zbar
for k=3:N
zk = z(:,k);

M = zk*zk.';        
zbar(:,k) = M(triu(true(length(zk))));
end

H=zeros(5);
H_old=zeros(5);
u_new=zeros(1,N);

for iter=1:100
    Phi=[];
    c=[];
    row=1;
    if iter>1
        if( norm(H-H_old)<1e-4)
            break;
        end
    end
    H_old=H;


%%giải pt bellman
row=1;
for k=3:N-1
    phi= zbar(:,k)-zbar(:,k+1);
    Phi(row,:)= phi';
    c(row) = y(k)^2 + u(k)^2;
    row=row+1;
 end    
Hbar= Phi \ c';


%chuyển từ Hbar về H
l = 5;
H = zeros(l,l);

idx = 1;
for i = 1:l
    for j = i:l
        H(i,j) = Hbar(idx);
        H(j,i) = Hbar(idx);   % vì H đối xứng
        idx = idx + 1;
    end
end

%%update policy mới
H_usigma= H(5,1:2);
H_uomega= H(5,3:4);
H_uu=H(5,5);
    for k=3:N-1
        sigma_k = z(1:2,k);
        omega_k= z(3:4,k);
        u_new(k)=-pinv(H_uu)*(H_usigma*sigma_k+ H_uomega*omega_k);
    end
    K = pinv(H_uu)*[H_usigma H_uomega];

    u=u_new;
    
H_error(iter) = norm(H - H_old);
end


s=zeros(4,N);
for k= 3:N
    s(:,k)=[u(k-1);u(k-2);y(k-1);y(k-2)];
end    


figure(1)
plot(H_error,'LineWidth',2)
xlabel('Iteration')
ylabel('||H_k - H_{k-1}||')
title('Convergence of H')
grid on

figure(2)
scatter(K*s(:,3:N),u(3:N),'.')
xlabel('K*s(k)')
ylabel('u(k)')
title('Check linear control law')
grid on


    
