#ifndef CPOINTSET_H
#define CPOINTSET_H

#include <vector>
#include <string>
#include "cpoint.h"
#include "cpoint3d.h"

using namespace std;

enum class dims {d2, d3};

template<class T>
class CPointSet:public vector<T>
{
public:
    CPointSet();
    CPointSet(const CPointSet &RHS);
    CPointSet(const string &fileName, int val_column=4);
    CPointSet<T>&  operator = (const CPointSet &RHS);
    void SetDimentions(dims _dim) {dimentions=_dim;}
    void WriteToVtp2D(const string &filename);
    void WriteToVtp3D(const string &filename);
    void WriteToPointsVtp(const string &filename, vector<double> limits = vector<double>());
    CPointSet<CPoint> MapToCylindrical(const double &x, const double &y);
    CPointSet<CPoint> MapToGrid(const double &_dx, const double &_dy, vector<double> span);
    CPointSet<T> Range(); //provide the range of pointset
    double KernelSmoothValue(const T &point, vector<double> span);
    double x(int i);
    double y(int i);
    double Value(int i);
private:
    dims dimentions = dims::d3;
};



#include "cpointset.hpp"

#endif // CPOINTSET_H
