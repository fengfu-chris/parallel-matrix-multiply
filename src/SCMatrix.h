class SCMatrix{
private:
	int rows;
	int cols;
	double* data;

public:
	SCMatrix();
	SCMatrix(int m);
	SCMatrix(int m, int n);
	SCMatrix(const SCMatrix& A);
	~SCMatrix();

	int getRows() const;
	int getCols() const;

	int operator==(const SCMatrix& A) const;
	int operator!=(const SCMatrix& A) const;
	SCMatrix& operator=(const SCMatrix& A);

	double* operator&() const;

	SCMatrix operator*(const SCMatrix& A);
	SCMatrix transpose() const;	

	double operator()(const int i, const int j) const;
	double& operator()(const int i, const int j);	

	void print() const;
	
	void init_with_RV();
	void init_with_Zero(); // default
};
