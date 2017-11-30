#include<iostream>
#include<vector>
#include<cmath>
#include<fstream>
#include<algorithm>
#include<iomanip>
#include "MVector.h"
#include "MMatrix.h"




int eq_no_func(int glob_node_no, int no_nodes) //returns the eq number from the global node number and the number of nodes.
{
	if (glob_node_no==1 || glob_node_no==no_nodes) return -1;
	else if (glob_node_no<1 || glob_node_no>no_nodes) 
	{
		std::cout<<"ERROR IN eq_no_func :: Global node number out of range."<<std::endl; 
		exit(1);
	}
	else return glob_node_no-1;
}

int glob_node_no_func(int loc_node_no, int element_no) //returns the global node number from the local node number and element nnumber
{	
	return element_no+loc_node_no-1;
}

int glob_eq_no_func(int local_eq_no, int element_no, int no_elements, int no_nodes) //is only used in phase 2a, when (1,1,100,101) is put in, 1 is returned. but in phase 2a, this is what is wanted since we are looping over n_dof, instead of local_eq_no.
{
	if (element_no==1) return 1;
	else if (element_no==no_elements) return no_nodes-2;
	else return element_no+local_eq_no-2;
}

int loc_eq_no_func(int loc_node_no, int element_no, int no_nodes)
{
	int j_dof=0, eq_no=0, glob_node_no=0, loc_eq_no=0;
	for (int j=1; j<=loc_node_no; j++)
	{
		glob_node_no=glob_node_no_func(j,element_no);
		eq_no=eq_no_func(glob_node_no,no_nodes);
		if (eq_no!=-1)
		{
			j_dof++;
			loc_eq_no=j_dof;
		}
		else loc_eq_no=-1;
		
	}
	return loc_eq_no;
}


double phi1(double x) {return 0.5*(1-x);}
double phi2(double x) {return 0.5*(1+x);}

double dphi1_ds() {return -0.5;}
double dphi2_ds() {return 0.5;}

double source_func(double x) {return 30*std::sin(sqrt(30)*x);}

int main()
{
	//Phase 1
	//Phase 1a
	int no_elements=100, no_nodes_per_element=2, no_nodes=101; // MAYBE THESE SHOULD BE GLOBAL VARIABLES: WOULDNT NEED TO PASS TO FUCNTIONS		
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
	
	//Phase 1c (inside eq_number)
	//Phase 1d
	std::vector<double> U_initial_guess (no_nodes,0);
	U_initial_guess[0]=0; U_initial_guess[no_nodes-1]=-1;
	
	
	//Phase 1e
	std::vector<double> eq_no_vec(no_nodes);
	for (int j=0;j<no_nodes; j++)
	{
		eq_no_vec[j]=eq_no_func(j+1,no_nodes);
		//std::cout<<eq_no_vec[j];
	}




	//Phase 1f Loc eq no is a fucntion also. Might not need this whole thing
	int j_dof,glob_node_no,eq_no;
	std::vector<double> glob_eq_no_vec(no_elements*no_nodes_per_element,0),loc_eq_no_vec(no_elements*no_nodes_per_element,0), N_dof_vec(no_elements);
	for (int e=1;e<=no_elements;e++)
	{
		j_dof=0;
		for(int j=1; j<=no_nodes_per_element; j++)
		{
			glob_node_no=glob_node_no_func(j,e);
			eq_no=eq_no_func(glob_node_no,no_nodes);
			if (eq_no!=-1)
			{
				j_dof++;
				glob_eq_no_vec[(e-1)*no_nodes_per_element+j-1]=eq_no;
				loc_eq_no_vec[(e-1)*no_nodes_per_element+j-1]=j_dof;
			}
			else loc_eq_no_vec[(e-1)*no_nodes_per_element+j-1]=-1;
		}
		N_dof_vec[e-1]=j_dof;
	}
	/*for (int i=0;i<no_elements*no_nodes_per_element;i++)
	std::cout<<glob_eq_no_vec[i]<<" , "<<loc_eq_no_vec[i]<<std::endl;
	for (int i=0; i<no_elements;i++) std::cout<<N_dof_vec[i]<<std::endl;
	std::cout<<glob_eq_no_vec[199]<<" , "<<loc_eq_no_vec[199]<<std::endl;
	std::cout<<glob_eq_no_vec[300]<<" , "<<loc_eq_no_vec.size()<<" , "<<N_dof_vec[105]<<std::endl;*/

	//Phase2
	//Phase 2a
	
	//Initialising 
	int n_dof,count=0;
	double X_point=0.0,f_X_point=0.0, dX_ds_point=0.0, dU_ds_point=0.0, dU_dX_point=0.0, Jacobianxs=0.0, i_loc_eq_no, j_loc_eq_no, i_glob_eq_no, j_glob_eq_no;
	double phi1_point=0.0, phi2_point=0.0, dphi1_dx_point=0.0, dphi2_dx_point=0.0, MULTIPLIER;
	std::vector<double> U=U_initial_guess;
	MVector element_residual_vec, glob_residual_vec(no_nodes-2,0);
	MMatrix glob_Jacobian_matrix(no_nodes-2,no_nodes-2,0);
	std::cout<<glob_residual_vec.size()<<std::endl;
	
	
	for (int e=1;e<=no_elements;e++) //loop over the elements
	{
		n_dof=N_dof_vec[e-1]; //gets n_dof for the current element
 		element_residual_vec.resize(n_dof, 0); //resets the element residual vector length to the length number of degrees of freedom for the element
		MMatrix element_Jacobian_matrix(n_dof,n_dof,0);

		//element_Jacobian_matrix.setRows(n_dof+1); element_Jacobian_matrix.setCols(n_dof+1); //resets the size of the Jacobian matrix to n_dof by n_dof
		//element_Jacobian_matrix=0;
		for (int i=0; i<gauss_position.size(); i++) //Loop over the gauss points within the element
		{
			X_point=0; //Reset these to zero so they can be incremented.
			dX_ds_point=0;
			dU_ds_point=0;

			//count++;
			glob_node_no=glob_node_no_func(1,e); //This is the function J(1,e), the global node no for the first node of this element

			X_point+=X[glob_node_no-1]*phi1(gauss_position[i]); // This is finding X_{J(1,e)} * phi 1
			dX_ds_point+=X[glob_node_no-1]*dphi1_ds(); // This is X_{(1,e)}* d phi_1/ds
			dU_ds_point+=U[glob_node_no-1]*dphi1_ds(); // This is U_{(1,e)}* d phi_1/ds

			glob_node_no=glob_node_no_func(2,e); //This is the function J(2,e), the global node no for the second node of this element
			X_point+=X[glob_node_no-1]*phi2(gauss_position[i]);
			dX_ds_point+=X[glob_node_no-1]*dphi2_ds(); // This is X_{(2,e)}* d phi_2/ds
			dU_ds_point+=U[glob_node_no-1]*dphi2_ds(); // This is U_{(2,e)}* d phi_2/ds

			dU_dX_point=(dU_ds_point)/(dX_ds_point); // This is dU/dx 
			Jacobianxs=dX_ds_point; //This is the Jacobian J^=dx/ds

			
			//The shape functions and their derivatives
			phi1_point=phi1(gauss_position[i]);
			phi2_point=phi2(gauss_position[i]);
			dphi1_dx_point=dphi1_ds()/Jacobianxs;
			dphi2_dx_point=dphi2_ds()/Jacobianxs;

			//std::cout<<count<<" , "<<X_point<<" , "<<dU_dX_point" , "<<dU_dX_point<<std::endl;
			//when we look at this output, we have 300 points with x from 0 to 1, dU/dx is zero for all the points, except the last 3, which have 					-100. I think this is right because U first guess is 0 across the whole way and -1 at X=1 for the BC. With these shape 					functions, the gradient is 0 everywhere and -100 between 0.99 and 1.

			f_X_point=source_func(X_point);
			
			for (int j=1; j<=no_nodes_per_element;j++) //Loop over local nodes
			{
				glob_node_no=glob_node_no_func(j,e); //Get global node number J(j,e)
				eq_no=eq_no_func(glob_node_no,no_nodes); //Get the equation number for this node
				if (eq_no!=-1) //If equaton no not equal to -1
				{
					i_loc_eq_no=loc_eq_no_func(j,e, no_nodes); //get local equation number of node L(j,e)
					//Increment the residual vector:
					if (j==1) 
					{
						element_residual_vec[i_loc_eq_no-1]+=(dU_dX_point*dphi1_dx_point+f_X_point*phi1_point)*Jacobianxs*gauss_weight[i];
					}
					if (j==2)
					{
						element_residual_vec[i_loc_eq_no-1]+=(dU_dX_point*dphi2_dx_point+f_X_point*phi2_point)*Jacobianxs*gauss_weight[i];
					}
					for (int k=1;k<=no_nodes_per_element;k++) //Loop over local nodes
					{
						glob_node_no=glob_node_no_func(k,e); //Get global node number J(k,e)
						eq_no=eq_no_func(glob_node_no,no_nodes); //Get the eq number for this node
						if (eq_no!=-1) // If eq no not -1
						{
							j_loc_eq_no=loc_eq_no_func(k,e,no_nodes); //Get local eq no L(k,e)
							//Get the correct dphi_k/dx*dphi_j/dx combo
							if (j==1 && k==1){MULTIPLIER=dphi1_dx_point*dphi1_dx_point;}
							else if (j==2 && k==2){MULTIPLIER=dphi2_dx_point*dphi2_dx_point;}
							else {MULTIPLIER=dphi1_dx_point*dphi2_dx_point;}
							//Increment the element Jacobian:
							element_Jacobian_matrix(i_loc_eq_no-1,j_loc_eq_no-1)+=MULTIPLIER*Jacobianxs*gauss_weight[i];
						}
					}
					
				}
			}
			//std::cout<<n_dof<<std::endl;
			//std::cout<<element_Jacobian_matrix<<std::endl;
			for(int i_dof=1; i_dof<=n_dof;i_dof++) //Loop over number of degrees of freedom for this element
			{
				i_glob_eq_no=glob_eq_no_func(i_dof,e,no_elements,no_nodes);
				std::cout<<i_glob_eq_no<<std::endl;
				glob_residual_vec[i_glob_eq_no-1]+=element_residual_vec[i_dof-1];
				//std::cout<<glob_residual_vec.size()<<std::endl;
				for (int j_dof=1; j_dof<=n_dof;j_dof++)
				{
					j_glob_eq_no=glob_eq_no_func(j_dof,e,no_elements,no_nodes);
					
					glob_Jacobian_matrix(i_glob_eq_no-1, j_glob_eq_no-1)+=element_Jacobian_matrix(i_dof-1,j_dof-1);
					//std::cout<<i_glob_eq_no<<" , "<<j_glob_eq_no<<" , "<<glob_Jacobian_matrix.Rows()<<" , "<<glob_Jacobian_matrix.Cols()<<std::endl;
				}
			}
				
		}
	//std::cout<<"element= "<<e<<"res="<<element_residual_vec<<"Jac="<<std::endl<<element_Jacobian_matrix<<std::endl;
	
	
	




	}

//std::cout<<glob_residual_vec<<std::endl;
//std::cout<<glob_Jacobian_matrix<<std::endl;
	std::ofstream filenameres("glob_residual_vec.txt");
	if(!filenameres){return 1;}
	filenameres<<glob_residual_vec;
	std::ofstream filenameJac("glob_Jacobian_matrix.txt");
	if(!filenameJac){return 1;}
	filenameJac<<glob_Jacobian_matrix;
	return 0;
}
