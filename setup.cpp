#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include<iomanip>
#include "MVector.h"
#include "MMatrix.h"




int equation_no_func(int glob_node_no, int no_nodes) //returns the equation number from the global node number and the number of nodes.
{
	if (glob_node_no==1 || glob_node_no==no_nodes) return -1;
	else if (glob_node_no<1 || glob_node_no>no_nodes) 
	{
		std::cout<<"ERROR IN equation_no_func :: Global node number out of range."<<std::endl; 
		exit(1);
	}
	else return glob_node_no-1;
}

int glob_node_no_func(int loc_node_no, int element_no) //returns the global node number from the local node number and element nnumber
{	
	return element_no+loc_node_no-1;
}

int glob_equation_no_func(int local_equation_no, int element_no, int no_elements, int no_nodes)
{
	if (element_no==1) return 1;
	else if (element_no==no_elements) return no_nodes-2;
	else return element_no+local_equation_no-2;
}
double phi1(double x) {return 0.5*(1-x);}
double phi2(double x) {return 0.5*(1+x);}

double dphi1() {return -0.5;}
double dphi2() {return 0.5;}

double source_func(double x) {return 30*std::sin(5.5*x);}

int main()
{
	//Phase 1
	//Phase 1a
	int no_elements=100, no_nodes_per_element=2, no_nodes=101;	
	std::vector<double> gauss_position(3), gauss_weight(3), shape_function1(3), shape_function2(3);
	double deriv_shape_function1=-0.5, deriv_shape_function2=0.5;
	gauss_position[0]=-sqrt(3.0/5.0); gauss_position[1]=0.0; gauss_position[2]=sqrt(3.0/5.0);
	gauss_weight[0]=5.0/9.0; gauss_weight[1]=8.0/9.0; gauss_weight[2]=5.0/9.0;
	for (int i=0;i<3;i++)
	{
		shape_function1[i]=0.5*(1-gauss_position[i]);
		shape_function2[i]=0.5*(1+gauss_position[i]);
	}
	
	//Phase 1b
	double X_0 = 0, X_N=1, domain_length=X_N-X_0, element_length=domain_length/no_elements;
	std::vector<double> X(no_nodes); X[0]=X_0; X[no_nodes-1]=X_N;
	//std::cout<<X[0]<<"  ";
	for (int j=1; j<no_nodes-1; j++)
	{
		X[j]=X[j-1]+element_length;
		//std::cout<<X[j]<<"  ";
	}
	//std::cout<<X[no_nodes-1]<<std::endl;
	
	//Phase 1c (inside equation_number)
	//Phase 1d
	std::vector<double> U_initial_guess (no_nodes,0);
	U_initial_guess[0]=0; U_initial_guess[no_nodes-1]=-1;
	
	//Phase 1e
	std::vector<double> equation_no_vec(no_nodes);
	for (int j=0;j<no_nodes; j++)
	{
		equation_no_vec[j]=equation_no_func(j+1,no_nodes);
		//std::cout<<equation_no_vec[j];
	}

	//Phase 1f
	int j_dof,glob_node_no,equation_no;
	std::vector<double> glob_equation_no_vec(no_elements*no_nodes_per_element,0),loc_equation_no_vec(no_elements*no_nodes_per_element,0), N_dof_vec(no_elements);
	for (int e=1;e<=no_elements;e++)
	{
		j_dof=0;
		for(int j=1; j<=no_nodes_per_element; j++)
		{
			glob_node_no=glob_node_no_func(j,e);
			equation_no=equation_no_func(glob_node_no,no_nodes);
			if (equation_no!=-1)
			{
				j_dof++;
				glob_equation_no_vec[(e-1)*no_nodes_per_element+j-1]=equation_no;
				loc_equation_no_vec[(e-1)*no_nodes_per_element+j-1]=j_dof;
			}
			else loc_equation_no_vec[(e-1)*no_nodes_per_element+j-1]=-1;
		}
		N_dof_vec[e-1]=j_dof;
	}
	/*for (int i=0;i<no_elements*no_nodes_per_element;i++)
	std::cout<<glob_equation_no_vec[i]<<" , "<<loc_equation_no_vec[i]<<std::endl;
	for (int i=0; i<no_elements;i++) std::cout<<N_dof_vec[i]<<std::endl;
	std::cout<<glob_equation_no_vec[199]<<" , "<<loc_equation_no_vec[199]<<std::endl;
	std::cout<<glob_equation_no_vec[300]<<" , "<<loc_equation_no_vec.size()<<" , "<<N_dof_vec[105]<<std::endl;*/

	//Phase2
	//Phase 2a
	
	int n_dof,count=0;
	double X_point=0.0,f_X_point=0.0, dX_ds_point=0.0, dU_ds_point=0.0, dU_dX_point=0.0;
	std::vector<double> global_residual_vec(glob_equation_no_vec.size(),0),element_residual_vec(1,0), U=U_initial_guess;
	MMatrix global_Jacobian_matrix(glob_equation_no_vec.size(),glob_equation_no_vec.size(),0), element_Jacobian_matrix;
	
	for (int e=1;e<=no_elements;e++) //loops over the elements
	{
		n_dof=N_dof_vec[e-1]; //gets n_dof for the current element
 		element_residual_vec.resize(n_dof, 0); //resets the element residual vector length to the length number of degrees of freedom for the element
		element_Jacobian_matrix.setRows(n_dof); element_Jacobian_matrix.setCols(n_dof); //resets the size of the Jacobian matrix to n_dof by n_dof
		for (int i=0; i<gauss_position.size(); i++)
		{
			//gauss_weight_point=gauss_weight[i];
			X_point=0;
			dX_ds_point=0;
			dU_ds_point=0;

			count++;
			glob_node_no=glob_node_no_func(1,e);
			//std::cout<<glob_node_no<<std::endl;
			X_point+=X[glob_node_no-1]*phi1(gauss_position[i]);
			dX_ds_point+=X[glob_node_no-1]*dphi1();
			dU_ds_point+=U[glob_node_no-1]*dphi1();

			glob_node_no=glob_node_no_func(2,e);
			X_point+=X[glob_node_no-1]*phi2(gauss_position[i]);
			dX_ds_point+=X[glob_node_no-1]*dphi2();
			dU_ds_point+=U[glob_node_no-1]*dphi2();

 
			dU_dX_point=(dU_ds_point)/(dX_ds_point);
			//std::cout<<count<<" , "<<X_point<<" , "<<dU_dX_point<<std::endl;

			f_X_point=source_func(X_point);
			
		}
	}

	
	return 0;
}
