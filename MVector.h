#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2

#include <vector>
#include<iostream>
#include<cmath>
#include<algorithm>
#include<fstream>
#include<string>

// Class that represents a mathematical vector
class MVector
{

//Operator overload for << like a  coordinate
/*friend std::ostream& operator<< (std::ostream& os, const MVector& output)
{
	os<<"(";
	for (int i=0;i<output.size()-1;i++)
	{
		os<<output[i]<<",";
	}
	os<<output[output.size()-1]<<")";
	return os;
}*/

//Operator overload for << like a col of values
friend std::ostream& operator<< (std::ostream& os, const MVector& output)
{
	for (int i=0;i<output.size();i++)
	{
		os.precision(8);
		os.width(12);
		os<<output[i]<<std::endl;
	}
	return os;
}

public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	//the explicit keyword prevents implicit conversions from type integer to MVector. without the explicit keyword (and using the commented out section at the bottom) then a is a length 2 vector when added
	//to y because of the implicit conversion. Adding the two does not change the first 2 cpts, so z[0]=y[0], and z[1]=y[1] i.e. where there is an a[i] value, but  z[2] is very strange.
	MVector(int n, double x) : v(n, x) {}

	// access element (lvalue) (see example sheet 5, q5.6)
	double &operator[](int index) { return v[index]; }

	// access element (rvalue) (see example sheet 5, q5.7)
	double operator[](int index) const { return v[index]; }

	int size() const { return v.size(); } // number of elements
	void resize(int newsize) // change number of elements
	{ 
		v.resize(newsize); 
	} 
	void resize(int newsize, double value) // change number of elements and set entries to a value
	{ 
		v.resize(newsize);
		for(int i=0;i<newsize;i++){v[i]=value;} 
	} 

	double LInfNorm() const
	{
		double maxAbs = 0;
		std::size_t s = size();
		for (int i=0; i<s; i++)
		{
			maxAbs = std::max(std::abs(v[i]), maxAbs);
		}
		return maxAbs;
	}
	
	double L2Norm() const
	{
		double rootsumsquares=0;
		for (int i=0;i<size();i++)
		{
			rootsumsquares+= v[i]*v[i];
		}
		rootsumsquares=pow(rootsumsquares,0.5);
		return rootsumsquares;
	}

private:
	std::vector<double> v;
};

// Operator overload for "scalar * vector"
inline MVector operator*(const double& lhs,const MVector& rhs)
{
	MVector temp(rhs);
	for (int i=0;i<temp.size();i++) temp[i]*=lhs;
	return temp;
}
// Operator overload for "vector * scalar"
inline MVector operator*(const MVector& lhs,const double& rhs)
{
	MVector temp(lhs);
	for (int i=0;i<temp.size();i++) temp[i]*=rhs;
	return temp;
}
//Operator overload for "vector + vector"
inline MVector operator+ (const MVector& lhs, const MVector& rhs)
{
	if (lhs.size()==rhs.size())
	{
		MVector temp(lhs);
		for (int i=0; i<temp.size();i++) temp[i]+=rhs[i];
		return temp;
	}
	else 
	{
		std::cout<<"ERROR: Attempted addition of two vectors of different length."<<std::endl;
		exit(1);
	}
} 
//Operator overload for "vector - vector"
inline MVector operator- (const MVector& lhs, const MVector& rhs)
{
	if (lhs.size()==rhs.size())
	{
		MVector temp(lhs);
		for (int i=0; i<temp.size();i++) temp[i]-=rhs[i];
		return temp;
	}
	else
	{
		std::cout<<"ERROR: Attempted subtraction of two vectors of different length."<<std::endl;
		exit(1);
	}
}
//Operator overload for "vector / scalar"
inline MVector operator/ (const MVector& lhs, const double& rhs)
{
	MVector temp(lhs);
	for (int i=0;i<temp.size();i++) temp[i]/=rhs;
	return temp;
}

double dot(const MVector& lhs, const MVector& rhs)
{
	double dotproduct=0;
	if (lhs.size()==rhs.size())
	{
		for (int i=0; i<lhs.size(); i++)
		{
			dotproduct+=lhs[i]*rhs[i];
		}
		return dotproduct;
	}
	else 
	{
		std::cout<<"ERROR: Attempted dot product of two vectors of different length."<<std::endl;
		exit(1);
	}
}
#endif
