load('glob_Jacobian_matrix.txt')
load('glob_residual_vec.txt')

U_init=[zeros(1,198)];

Y=-transpose(glob_Jacobian_matrix\glob_residual_vec);


lilJac=glob_Jacobian_matrix(100:198,100:198);
lilres=glob_residual_vec(100:198);
lilY=-transpose(lilJac\lilres);

v_init=zeros(1,99);
v_new=v_init+lilY;

U_new=U_init+Y;

u_new=[0,U_new(1:99),1];
v_new=[0,v_new,-1];

X_guess=0:0.01:1;

u_exact=(sin(sqrt(30))-1).*X_guess.^3./6+(1/30)*sin(sqrt(30).*X_guess)+(7/6-1/5*sin(sqrt(30)))*X_guess;
v_exact=(sin(sqrt(30))-1).*X_guess-sin(sqrt(30).*X_guess);

[valu,indexu]=max(abs(u_exact-u_new));
max_pc_differenceu=abs((valu*100)/u_exact(indexu));

[valv,indexv]=max(abs(v_exact-v_new));
max_pc_differencev=abs((valv*100)/v_exact(indexv));



plot(X_guess,u_new,'m*')
hold on 
plot(X_guess,v_new,'r*')
plot(X_guess,u_exact,'g')
plot(X_guess,v_exact,'k')

title(['maximum percentage difference is ',num2str(max_pc_differenceu),'% [u] and ',num2str(max_pc_differencev),'% [v]'])