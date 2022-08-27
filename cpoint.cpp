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
