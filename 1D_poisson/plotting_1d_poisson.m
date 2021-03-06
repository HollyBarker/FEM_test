load('glob_Jacobian_matrix.txt');
load('glob_residual_vec.txt');
Y=-glob_Jacobian_matrix\glob_residual_vec;
U_init=zeros(99,1);

U_new=U_init+Y;
U_new=[0;U_new;-1];
X_guess=0:0.01:1;
X_exact=0:0.01:1;
U_exact=(sin(sqrt(30))-1).*X_exact-sin(sqrt(30).*X_exact);

plot(X_guess,U_new,'*');
hold on ;
plot(X_exact,U_exact);
legend('FEM solution','Exact solution');

[val,index]=max(abs(transpose(U_exact)-U_new));
max_pc_difference=abs((val*100)/U_exact(index));

title(['maximum percentage difference is ',num2str(max_pc_difference),'%'])
