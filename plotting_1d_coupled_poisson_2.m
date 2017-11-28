%load('glob_Jacobian_matrix.txt')
load('glob_residual_vec.txt')

U_init=[zeros(1,198)]

Y=-transpose(glob_Jacobian_matrix\glob_residual_vec);


lilJac=glob_Jacobian_matrix(100:198,100:198)
lilres=glob_residual_vec(100:198)
lilY=-transpose(lilJac\lilres)

v_init=zeros(1,99)
v_new=v_init+lilY

U_new=U_init+Y;

u_new=[0,U_new(1:99),1]
v_new=[0,v_new,-1]


X_guess=0:0.01:1;

plot(X_guess,u_new,'k')
hold on 
plot(X_guess,v_new,'r')

