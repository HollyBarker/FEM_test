load('glob_Jacobian_matrix_U.txt'); %These need to be loaded after one iteration has already been performed
load('glob_Jacobian_matrix_V.txt'); %already done with U_init is zeros and v_init zeros
load('glob_residual_vec_U.txt');    %the solution for which is read into cpp file with readFile function
load('glob_residual_vec_V.txt'); %after pasting the new u and v into a txt file

U_init=load('U_vector3.txt'); %Should be the resultant after 1 iteration already completed
V_init=load('V_vector3.txt');
X_guess=0:0.01:1;

YU=-transpose(glob_Jacobian_matrix_U\glob_residual_vec_U);
YV=-transpose(glob_Jacobian_matrix_V\glob_residual_vec_V);

YU=[0,YU,0];
YV=[0,YV,0];

U_new=U_init+YU;
V_new=V_init+YV;



plot(X_guess,U_new,'g');
hold on
plot(X_guess,V_new,'y');