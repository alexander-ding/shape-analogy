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

void Mesh::center() {
    Eigen::VectorXf centroid = this->m_vertices.rowwise().mean().transpose();

    // Step 2: Translate vertices to center them around the origin
    Eigen::MatrixXf centeredVertices = this->m_vertices.colwise() - centroid;

    // Step 3: Scale vertices to fit within [-1, 1]
    float maxAbs = centeredVertices.cwiseAbs().maxCoeff(); // Get the maximum absolute value
    centeredVertices /= maxAbs; // Scale to fit within the range

    this->m_vertices = centeredVertices;
}

Mesh::Mesh(const MatrixXf &vertices, const MatrixXi &faces, bool center)
{
    this->init(vertices, faces);
    if (center) this->center();
}

Mesh::Mesh(const QString &path, bool center)
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
    if (center) this->center();
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

MatrixXf Mesh::computeFaceNormals()
{
    MatrixXf faceNormals(3, this->m_faces.cols());
    for (size_t i = 0; i < this->m_faces.cols(); ++i)
    {
        // Assuming each face has three indices: idx1, idx2, idx
        Vector3i indices = this->m_faces.col(i);

        // Get the vertices from the indices
        Vector3f v1 = this->m_vertices.col(indices[0]);
        Vector3f v2 = this->m_vertices.col(indices[1]);
        Vector3f v3 = this->m_vertices.col(indices[2]);

        // Compute the two edges from the vertices
        Vector3f edge1 = v2 - v1;
        Vector3f edge2 = v3 - v1;

        // Cross product of the two edges gives the normal
        Vector3f normal = edge1.cross(edge2);
        normal.normalize();

        // Store the computed normal
        faceNormals.col(i) = normal;
    }
    return faceNormals;
}


MatrixXf Mesh::computeVertexNormals()
{
    // Create a MatrixXf to store vertex normals. Initially, fill it with zeros.
    MatrixXf vertexNormals(3, this->m_vertices.cols());
    vertexNormals.setZero();

    // Compute or retrieve face normals (assuming computeFaceNormals is available and correct)
    MatrixXf faceNormals = this->computeFaceNormals();

    // Accumulate normals for each vertex from each face
    for (size_t i = 0; i < this->m_faces.cols(); ++i)
    {
        Vector3i indices = this->m_faces.col(i);
        Vector3f normal = faceNormals.col(i);

        // Add this normal to each of the vertices that make up the face
        for (int j = 0; j < 3; ++j)  // Assuming triangular faces
        {
            vertexNormals.col(indices[j]) += normal;
        }
    }

    // Normalize all vertex normals to convert them into unit vectors
    for (size_t i = 0; i < vertexNormals.cols(); ++i)
    {
        vertexNormals.col(i).normalize();
    }

    return vertexNormals;
}
