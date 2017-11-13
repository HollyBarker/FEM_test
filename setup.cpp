#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>

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

int glob_equation_no_func(int local_equation_no, int element_no, int no_elements)
{
	if (element_no==1) return 1;
	else if (element_no==no_elements) return no_nodes-2;
	else return element_no+local_equation_no-2;
}

int main()
{
	//Phase 1
	//Phase 1a
	int no_elements=100, no_nodes_per_element=2, no_nodes=101;	
	std::vector<double> gauss_position(3), gauss_weight(3), shape_function1(3), shape_function2(3);
	double deriv_shape_function1=-1/2, deriv_shape_function2=1/2;
	gauss_position[0]=-sqrt(3/5); gauss_position[1]=0; gauss_position[2]=sqrt(3/5);
	gauss_weight[0]=5/9; gauss_weight[1]=8/9; gauss_weight[2]=5/9;
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
	U_initial_guess[0]=0; U_initial_guess[no_nodes-1]=1;
	
	//Phase 1e
	std::vector<double> equation_no_vec(no_nodes);
	for (int j=0;j<no_nodes; j++)
	{
		equation_no_vec[j]=equation_no_func(j+1,no_nodes);
		//std::cout<<equation_no_vec[j];
	}

	//Phase 1f LOOK AT THIS	- DO WE NEED A MAP? FUNCTION? MATRIX? OTHER LOOK-UP SCHEMES ALSO NEED CHECKING
	int j_dof, glob_node_no;
	std::vector<double> glob_equation_no_vec(no_nodes+no_nodes-no_elements_per_node), loc_equation_no_vec(no_nodes+no_nodes-no_elements_per_node),N_dof_vec(no_elements);
	for(int e=1;e<=no_elements;e++)
	{
		j_dof=0;
		for(int j=1;j<=no_nodes_per_element;j++)
		{
			glob_node_no=glob_node_no_func(j,e);
			if (equation_no_vec[e+j-2]!=-1) 
			{
				j_dof++;
				glob_equation_no_vec[e*no_nodes_per_element+j-4]=equation_no_vec[e+j-2];
				loc_equation_no_vec[e*no_nodes_per_element-1]=j; //when i put j_dof, it would set the local equation number to 1, then the next one to 2, but when this node is reviewed in the next element, its local equation number is set back to 1
				std::cout<<loc_equation_no_vec[e+j-2]<<std::endl;
			}
			else loc_equation_no_vec[e+j-2]=-1;
		}
		N_dof_vec[e-1]=j_dof;
		//std::cout<<N_dof_vec[e-1]<<std::endl;
	}
	for (int i=0;i<no_nodes;i++)
	std::cout<<glob_equation_no_vec[i]<<" , "<<loc_equation_no_vec[i]<<std::endl;
	//for (int i=0;i<no_elements;i++)
	//std::cout<<N_dof_vec[i]<<std::endl;
	//std::cout<<N_dof_vec.size()<<std::endl;
	return 0;
}
