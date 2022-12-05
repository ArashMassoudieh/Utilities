#include "cpoint.h"

CPoint::CPoint()
{
    coordinate.resize(2);
}

CPoint::CPoint(const CPoint& P)
{
    coordinate = P.coordinate;
    values = P.values;
}
CPoint& CPoint::operator = (const CPoint &P)
{
    coordinate = P.coordinate;
    values = P.values;
    return *this;
}

CPoint::CPoint(const double &x,const double &y)
{
    coordinate.resize(2);
    coordinate[0] = x;
    coordinate[1] = y;
}

double CPoint::x() const
{
    return coordinate[0];
}
double CPoint::y() const
{
    return coordinate[1];
}
void CPoint::setx(const double &val)
{
    coordinate[0] = val;
}
void CPoint::sety(const double &val)
{
    coordinate[1] = val;
}

double CPoint::Value(int i)
{
    if (i<values.size())
    {
        return values[i];
    }
    else
        return 0;
}
