#include "cpoint3d.h"

CPoint3d::CPoint3d()
{
    coordinate.resize(3);
}

CPoint3d::CPoint3d(const CPoint3d& P)
{
    coordinate = P.coordinate;
    values = P.values;
}
CPoint3d& CPoint3d::operator = (const CPoint3d &P)
{
    coordinate = P.coordinate;
    values = P.values;
    return *this;
}

CPoint3d::CPoint3d(const double &x,const double &y, const double &z)
{
    coordinate.resize(3);
    coordinate[0] = x;
    coordinate[1] = y;
    coordinate[2] = z;
}

double CPoint3d::z() const
{
    return coordinate[2];
}

void CPoint3d::setz(const double &val)
{
    coordinate[2] = val;
}
