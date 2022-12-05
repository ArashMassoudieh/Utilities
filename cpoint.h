#ifndef CPOINT_H
#define CPOINT_H

#include <vector>

enum class range {max,min};
enum class dir {x,y};

using namespace std;
class CPoint
{
public:
    CPoint();
    CPoint(const CPoint&);
    CPoint(const double &x,const double &y);
    CPoint& operator = (const CPoint &C);
    double x() const;
    double y() const;
    void setx(const double &val);
    void sety(const double &val);
    void AppendValue(const double &val) {values.push_back(val);}
    double Value(int i);

protected:
    vector<double> coordinate;
    vector<double> values;
};

#endif // CPOINT_H
