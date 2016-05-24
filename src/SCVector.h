#include"SCMatrix.h"

class SCVector{
private:
	int dimension;
	double* data;

public:
	SCVector(int dim);
	SCVector(const SCVector& v);
	SCVector(int col, const SCMatrix& A);
	~SCVector();

	int Dimension() const;
	double Length();  /* Euclidean Norm of the Vector */
	void Normalize();

	double Norm_l1();
	double Norm_l2();
	double Norm_linf();
	//double MaxMod();
	//int MaxModindex();
	int operator==(const SCVector& v) const;
	int operator!=(const SCVector& v) const;
	SCVector& operator=(const SCVector& v);
	
	double operator()(const int i) const;
	double& operator()(const int i);

	void Print() const;
	void Initialize(double a);
	void Initialize(double *v);
};
