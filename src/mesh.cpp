#include "mesh.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "util/tiny_obj_loader.h"
#include <QString>
#include <QFileInfo>
#include <Eigen/Dense>
#include <iostream>
#include <QTextStream>

using namespace std;
using namespace Eigen;

void Mesh::init(const MatrixXf &vertices, const MatrixXi &faces)
{
    this->m_vertices = vertices;
    this->m_faces = faces;
}

Mesh::Mesh(const MatrixXf &vertices, const MatrixXi &faces)
{
    this->init(vertices, faces);
}

Mesh::Mesh(const QString &path)
{
    tinyobj::attrib_t attrib;
    vector<tinyobj::shape_t> shapes;
    vector<tinyobj::material_t> materials;

    QFileInfo info(path);
    string err;
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err,
                                info.absoluteFilePath().toStdString().c_str(), (info.absolutePath().toStdString() + "/").c_str(), true);
    if (!err.empty())
    {
        cerr << err << endl;
    }

    if (!ret)
    {
        cerr << "Failed to load/parse .obj file" << endl;
        return;
    }

    MatrixXf vertices(3, attrib.vertices.size() / 3);
    for (size_t i = 0; i < attrib.vertices.size() / 3; i++)
    {
        vertices.col(i) = Vector3f(attrib.vertices[3*i], attrib.vertices[3*i + 1], attrib.vertices[3*i + 2]);
    }

    size_t nFaces = 0;
    for (size_t s = 0; s < shapes.size(); s++)
    {
        nFaces += shapes[s].mesh.num_face_vertices.size();
    }
    MatrixXi faces(3, nFaces);

    size_t faceI = 0;
    for (size_t s = 0; s < shapes.size(); s++)
    {
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++)
        {
            unsigned int fv = shapes[s].mesh.num_face_vertices[f];

            Vector3i face;
            for (size_t v = 0; v < fv; v++)
            {
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];

                face[v] = idx.vertex_index;
            }
            faces.col(faceI) = face;
            faceI += 1;
            index_offset += fv;
        }
    }

    this->init(vertices, faces);
    cout << "Loaded " << faces.cols() << " faces and " << vertices.cols() << " vertices" << endl;
}

Mesh::~Mesh()
{
}

void Mesh::saveToFile(const QString &path)
{
    QFile outfile(path);
    if (!outfile.open(QIODevice::ReadWrite)) {
        cerr << "Failed to open file " << path.toStdString() << endl;
    }
    QTextStream stream(&outfile);

    for (size_t i = 0; i < this->m_vertices.cols(); i++)
    {
        stream << "v " << this->m_vertices(0, i) << " " << this->m_vertices(1, i) << " " << this->m_vertices(2, i) << "\n";
    }

    for (size_t i = 0; i < this->m_faces.cols(); i++)
    {
        stream << "f " << this->m_faces(0, i) + 1 << " " << this->m_faces(1, i) + 1 << " " << this->m_faces(2, i) + 1 << "\n";
    }

    outfile.close();
}
