#ifndef CPOINT3D_H
#define CPOINT3D_H

#include <cpoint.h>

class CPoint3d : public CPoint
{
public:
    CPoint3d();
    CPoint3d(const CPoint3d& P);
    CPoint3d& operator = (const CPoint3d &P);
    CPoint3d(const double &x,const double &y, const double &z);
    double z() const;
    void setz(const double &val);

};

#endif // CPOINT3D_H
