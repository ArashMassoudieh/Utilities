// Matrix.h: interface for the CMatrix class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include <iostream>
#include <math.h>
#define ARMA_DONT_PRINT_ERRORS
#ifdef _arma
#include "armadillo"
#include "Matrix_arma.h"
#include "Vector_arma.h"
#include "Matrix_arma_sp.h"
using namespace arma;
#endif
class QVariant;
//class QString;
//class QList;
#ifdef QT_version
#include <QMap>
#endif // QT_version

using namespace std;

class CVector;
class CMatrix
{
friend class D5Matrix;
private:
	int numrows;
	int numcols;
	int range(int);
public:
	vector<CVector> matr;
	CMatrix(int, int);
	CMatrix(int);
    static CMatrix Diag(int n);
	CMatrix();
	CMatrix(string filename);
	CMatrix(const CMatrix&);
#ifdef _arma
    CMatrix(CMatrix_arma_sp&);
    CMatrix(const CMatrix_arma &M);
    CMatrix& operator=(mat&);
#endif
	CMatrix(const CVector&);
	CVector& operator[](int);
    double & operator()(int i, int j);
    double & operator()(unsigned int i, unsigned int j);
	int getnumrows() const;
	int getnumcols() const;
	virtual ~CMatrix();
	CMatrix& operator=(const CMatrix&);
	CMatrix & operator=(const double & m);
	CMatrix& operator+=(const CMatrix&);
	CMatrix& operator-=(const CMatrix &);

	friend CMatrix mult(CMatrix&, CMatrix&);
	friend CVector mult(CMatrix&, CVector&);
	friend CVector mult(CVector&, CMatrix&);
	friend void triangulate(CMatrix&, CVector&);
	friend void backsubst(CMatrix&, CVector&, CVector&);
	friend CVector gauss0(CMatrix, CVector);
    friend CVector diag(const CMatrix&);
    CVector maxelements() const;
	friend CMatrix Cholesky_factor(CMatrix &M);
	friend CMatrix LU_decomposition(CMatrix &M);
    CMatrix LU_decomposition();
    CMatrix Cholesky_factor();
    double det();
    void Print(FILE *FIL);
    void print(string s);
    void setval(double a);
    void setvaldiag(double a);
    void writetofile(FILE *f);
    void writetofile(string filename);
    void writetofile_app(string filename);
	friend void write_to_file(vector<CMatrix> M, string filename);
	friend CMatrix Average(vector<CMatrix> M);
    CVector diag_ratio();
    vector<vector<bool> > non_posdef_elems(double tol = 1);
    CMatrix non_posdef_elems_m(double tol = 1);
    CMatrix Preconditioner(double tol = 1);
    vector<string> toString(string format = "", vector<string> columnHeaders = vector<string>(), vector<string> rowHeaders = vector<string>()) const;
	vector<string> toHtml(string format = "", vector<string> columnHeaders = vector<string>(), vector<string> rowHeaders = vector<string>());
    void setnumcolrows();
	void ScaleDiagonal(double x);
    void setcol(int i,  const CVector &V);
    void setrow(int j,  const CVector &V);
#ifdef QT_version
    QMap<QString, QVariant> compact() const;
    static CMatrix unCompact(QMap<QString, QVariant>);
#endif // QT_version

};

double det(const CMatrix &);
CMatrix Log(const CMatrix &M1);
CMatrix Exp(const CMatrix &M1);
CMatrix Sqrt(const CMatrix &M1);
CMatrix operator+(const CMatrix&, const CMatrix&);
CMatrix operator+(double, CMatrix);
CMatrix operator+(CMatrix, double);
CMatrix operator-(double d, CMatrix m1);
CMatrix operator+(CMatrix m1, double d);
CMatrix operator-(CMatrix m1,double d);
CMatrix operator/(CMatrix m1,double d);
CMatrix operator/(double d, CMatrix m1);
CMatrix operator-(const CMatrix&, const CMatrix&);
CMatrix operator*(CMatrix, CMatrix);
CVector operator*(CMatrix, CVector);
CMatrix operator*(CVector, CMatrix);
CMatrix operator*(double, CMatrix);
CVector operator/(CVector&, CMatrix&);
CVector operator/(const CVector &V, const CMatrix &M);
CMatrix Transpose(CMatrix &M1);
#ifdef  arma
CMatrix Invert(CMatrix M1);
CVector SpareSolve(CMatrix, CVector);
#endif //  arma


CMatrix oneoneprod(CMatrix &m1, CMatrix &m2);
#ifdef _arma
CVector solve_ar(const CMatrix&, const CVector&);
#endif
CMatrix inv(CMatrix&);

CMatrix normalize_diag(CMatrix&, CMatrix&);
CVector normalize_diag(CVector&, CMatrix&);
CVector normalize_diag(const CVector &V, const CVector &D);
CVector normalize_diag(const CVector &V, const CMatrix &M2);
CMatrix normalize_max( const CMatrix &M1, const CMatrix &M2);
CVector normalize_max( const CVector &V, const CMatrix &M2);
CVector normalize_max( const CVector &V, const CVector &D);

CMatrix Identity(int rows);


