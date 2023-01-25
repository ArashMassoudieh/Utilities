#include "cpointset.h"
#include "fstream"
#include "Utilities.h"
#include "cpoint3d.h"
#include "cpoint.h"
#include "VTK.h"
#include "vtkUnstructuredGrid.h"

template<class T>
CPointSet<T>::CPointSet():vector<T> ()
{

}

template<class T>
CPointSet<T>::CPointSet(const CPointSet &RHS):vector<T> (RHS)
{

}

template<class T>
CPointSet<T>& CPointSet<T>::operator = (const CPointSet &RHS)
{
    vector<T>::operator = (RHS);
    return *this;
}

template<class T>
CPointSet<T>::CPointSet(const string &fileName)
{
    ifstream file(fileName);
    if (!file.good()) return;
    int counter=0;
    if (dimentions == dims::d3)
    {
        while (!file.eof())
        {   vector<string> vals = aquiutils::getline(file);
            if (counter>0 && vals.size()>=5)
            {
                vector<double> vals_d = aquiutils::ATOF(vals);
                CPoint3d P(vals_d[1],vals_d[2],vals_d[3]);
                P.AppendValue(vals_d[4]);
                vector<T>::push_back(P);
            }
            counter++;
        }
    }

}

template<class T>
void CPointSet<T>::WriteToVtp3D(const string &filename)
{
    if (dimentions==dims::d2) return;
    vtkSmartPointer<vtkPoints> points_3 =
        vtkSmartPointer<vtkPoints>::New();

    double xx, yy, zz;
    vtkSmartPointer<vtkFloatArray> values =
        vtkSmartPointer<vtkFloatArray>::New();

    values->SetNumberOfComponents(1);

    values->SetName("Moisture Content");

    for (unsigned int i = 0; i < vector<CPoint3d>::size(); i++)
    {
        xx = vector<CPoint3d>::at(i).x();
        yy = vector<CPoint3d>::at(i).y();
        zz = vector<CPoint3d>::at(i).z();
        float t[1] = { float(vector<CPoint3d>::at(i).Value(0)) };
        points_3->InsertNextPoint(xx, yy, zz);
        values->InsertNextTupleValue(t);

    }

    // Add the grid points to a polydata object
    vtkSmartPointer<vtkPolyData> inputPolyData =
        vtkSmartPointer<vtkPolyData>::New();
    inputPolyData->SetPoints(points_3);
    vtkUnstructuredGrid* outputPolyData;
    // Triangulate the grid points

    vtkSmartPointer<vtkDelaunay3D> delaunay =
        vtkSmartPointer<vtkDelaunay3D>::New();
#if VTK_MAJOR_VERSION <= 5
    delaunay->SetInput(inputPolyData);
#else
    delaunay->SetInputData(inputPolyData);
#endif
    delaunay->Update();
    outputPolyData = delaunay->GetOutput();

    outputPolyData->GetPointData()->SetScalars(values);

    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(filename.c_str());
    writer->SetInputData(outputPolyData);
    writer->Write();

}

template<class T>
void CPointSet<T>::WriteToVtp2D(const string &filename)
{
    if (dimentions==dims::d3) return;
    vtkSmartPointer<vtkPoints> points_3 =
        vtkSmartPointer<vtkPoints>::New();

    double xx, yy, zz;
    vtkSmartPointer<vtkFloatArray> values =
        vtkSmartPointer<vtkFloatArray>::New();

    values->SetNumberOfComponents(1);

    values->SetName("Moisture Content");

    for (unsigned int i = 0; i < vector<CPoint>::size(); i++)
    {
        xx = vector<CPoint>::at(i).x();
        yy = vector<CPoint>::at(i).y();
        zz = 0;
        float t[1] = { float(vector<CPoint>::at(i).Value(0)) };
        points_3->InsertNextPoint(xx, yy, zz);
        values->InsertNextTupleValue(t);

    }

    // Add the grid points to a polydata object
    vtkSmartPointer<vtkPolyData> inputPolyData =
        vtkSmartPointer<vtkPolyData>::New();
    inputPolyData->SetPoints(points_3);
    vtkPolyData* outputPolyData;
    // Triangulate the grid points

    vtkSmartPointer<vtkDelaunay2D> delaunay =
        vtkSmartPointer<vtkDelaunay2D>::New();
#if VTK_MAJOR_VERSION <= 5
    delaunay->SetInput(inputPolyData);
#else
    delaunay->SetInputData(inputPolyData);
#endif
    delaunay->Update();
    outputPolyData = delaunay->GetOutput();

    outputPolyData->GetPointData()->SetScalars(values);


    //Append the two meshes
    vtkSmartPointer<vtkAppendPolyData> appendFilter =
        vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
    appendFilter->AddInputConnection(input1->GetProducerPort());
    appendFilter->AddInputConnection(input2->GetProducerPort());
#else
    //appendFilter->AddInputData(polydata);
    //appendFilter->AddInputData(polydata_1);
    appendFilter->AddInputData(outputPolyData);
#endif
    appendFilter->Update();


    // Visualization
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
    mapper->SetInputConnection(polydata->GetProducerPort());
#else
    mapper->SetInputConnection(appendFilter->GetOutputPort());
    //mapper->SetInputData(polydata_1);
#endif


vtkSmartPointer<vtkActor> actor =
        vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(5);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(mapper->GetInput());
    // This is set so we can see the data in a text editor.
    writer->SetDataModeToAscii();
    writer->Write();


}


template<class T>
void CPointSet<T>::WriteToPointsVtp(const string &filename, vector<double> limits)
{
    if (limits.size()==0)
    {
        limits.resize(6);
        limits[0]=-10000;
        limits[1]=10000;
        limits[2]=-10000;
        limits[3]=10000;
        limits[4]=-10000;
        limits[5]=-10000;
    }
    double xx, yy, zz;
    vtkSmartPointer<vtkFloatArray> values =
        vtkSmartPointer<vtkFloatArray>::New();

    values->SetNumberOfComponents(1);

    values->SetName("Moisture Content");

    vtkNew<vtkPoints> points;


    // Create the topology of the point (a vertex)
    vtkNew<vtkCellArray> vertices;
    // We need an an array of point id's for InsertNextCell.
    int counter=0;
    for (unsigned int i = 0; i < vector<T>::size(); i++)
    {
        if (vector<T>::at(i).x()>limits[0] && vector<T>::at(i).x()<limits[1] && vector<T>::at(i).y()>limits[2] && vector<T>::at(i).y()<limits[3] && vector<T>::at(i).z()>limits[4] && vector<T>::at(i).z()<limits[5])
        {
            counter++;
        }
    }

    vtkIdType* pid = new vtkIdType[counter];
    counter=0;
    for (unsigned int i = 0; i < vector<T>::size(); i++)
    {
        if (vector<T>::at(i).x()>limits[0] && vector<T>::at(i).x()<limits[1] && vector<T>::at(i).y()>limits[2] && vector<T>::at(i).y()<limits[3] && vector<T>::at(i).z()>limits[4] && vector<T>::at(i).z()<limits[5])
        {   const float p[3] = {vector<T>::at(i).x(), vector<T>::at(i).y(), vector<T>::at(i).z()};
            pid[counter] = points->InsertNextPoint(p);
            float t[1] = { float(vector<T>::at(i).Value(0)) };
            values->InsertNextTupleValue(t);
            counter++;
        }
    }

    vertices->InsertNextCell(counter, pid);
    // Create a polydata object
      vtkNew<vtkPolyData> point;

    // Set the points and vertices we created as the geometry and topology of the
    // polydata
    point->SetPoints(points);
    point->SetVerts(vertices);

    // Add the grid points to a polydata object
    vtkSmartPointer<vtkPolyData> inputPolyData =
        vtkSmartPointer<vtkPolyData>::New();
    inputPolyData->SetPoints(points);

    point->GetPointData()->SetScalars(values);


    //Append the two meshes
    vtkSmartPointer<vtkAppendPolyData> appendFilter =
        vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
    appendFilter->AddInputConnection(input1->GetProducerPort());
    appendFilter->AddInputConnection(input2->GetProducerPort());
#else
    //appendFilter->AddInputData(polydata);
    //appendFilter->AddInputData(polydata_1);
    appendFilter->AddInputData(point);
#endif
    appendFilter->Update();


    // Visualization
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
    mapper->SetInputConnection(polydata->GetProducerPort());
#else
    mapper->SetInputConnection(appendFilter->GetOutputPort());
    //mapper->SetInputData(polydata_1);
#endif


vtkSmartPointer<vtkActor> actor =
        vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(5);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(mapper->GetInput());
    // This is set so we can see the data in a text editor.
    writer->SetDataModeToAscii();
    writer->Write();


}

template<class T>
CPointSet<CPoint> CPointSet<T>::MapToCylindrical(const double &_x, const double &_y)
{
    CPointSet<CPoint> out;
    for (int i=0; i<vector<T>::size(); i++)
    {
        CPoint P(sqrt(pow(vector<T>::at(i).x()-_x,2)+pow(vector<T>::at(i).y()-_y,2)),pow(vector<T>::at(i).z()));
        P.AppendValue(vector<T>::at(i).Value(0));
        out.push_back(P);

    }
    return out;
}

template<class T>
double CPointSet<T>::KernelSmoothValue(const T &point, vector<double> span)
{
    double sumnumerator=0;
    double sumdenominator=0;
    if (dimentions==dims::d2)
    {
        for (int i=0; i<vector<T>::size(); i++)
        {
            sumnumerator += exp(-pow(point.x()-x(i),2)/(2*span[0]*span[0])-pow(point.y()-y(i),2)/((2*span[1]*span[1])))*Value(i);
            sumdenominator += exp(-pow(point.x()-x(i),2)/(2*span[0]*span[0])-pow(point.y()-y(i),2)/((2*span[1]*span[1])));
        }
        return  sumnumerator/sumdenominator;
    }
}


template<class T>
double CPointSet<T>::x(int i)
{
    return vector<T>::at(i).x();
}


template<class T>
double CPointSet<T>::y(int i)
{
    vector<T>::at(i).y();
}

template<class T>
double CPointSet<T>::Value(int i)
{
    return vector<T>::at(i).Value()[0];
}


