#include<iostream>
#include<fstream>
#include<string>
#include<vector>


void readFile(std::vector<double> &vec, std::string strFile) //Read the file into the vector function definition
{
	double number; //First and last name
	std::ifstream vecfile;
	vecfile.open(strFile.c_str()); //Opens file
	while (vecfile >> number) //While the file is copying into the first and last names 
	{
		vec.push_back(number); //Push the names onto the back of the of the vector
	}

	vecfile.close(); //Close the input file
}

int main()
{
	std::vector<double> numbers;
	std::ifstream vecfile;
	readFile(numbers,"readfile.txt");
	for (int i=0;i<numbers.size();i++) {std::cout<<numbers[i]<<" , "<<std::endl;}
return 0;
}
