#pragma once

#include <Eigen/Dense>
#include <QString>

using namespace std;
using namespace Eigen;

class Mesh
{
public:
    Mesh(const MatrixXf &vertices, const MatrixXi &faces);
    Mesh(const QString &path);
    ~Mesh();

    void saveToFile(const QString& path);

    MatrixXf &getVertices()
    {
        return this->m_vertices;
    }

    void setVertices(const Eigen::MatrixXf &vertices);

    MatrixXi &getFaces()
    {
        return this->m_faces;
    }

private:
    void init(const MatrixXf &vertices, const MatrixXi &faces);

    MatrixXf m_vertices;
    MatrixXi m_faces;
};
