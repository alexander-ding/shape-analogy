#include "analogy.h"
#include <iostream>
#include <limits>

Analogy::Analogy(Mesh aPrime, Mesh b)
    : m_aPrime(aPrime), m_b(b), m_bPrime(b)
{

}

void Analogy::computeBPrime(float lambda)
{
    Eigen::MatrixXf& vertices = this->m_b.getVertices();
    Eigen::MatrixXf& newVertices = this->m_bPrime.getVertices();
    Eigen::MatrixXf prevVertices = newVertices;
    this->m_bTargetNormals = this->computeTargetNormals();
    std::vector<std::unordered_map<size_t, float> > edgeWeights = this->computeEdgeWeights();
    Eigen::VectorXf voronoiAreas = this->computeVoronoiAreas();
    std::unique_ptr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>>> solver = this->buildL(edgeWeights);
    Eigen::MatrixXf normals = this->m_b.computeVertexNormals();

    size_t nVertices = edgeWeights.size();
    MatrixXf b(newVertices.cols(), 3);
    std::vector<Eigen::Matrix3f> Rs(edgeWeights.size(), Eigen::Matrix3f::Identity());

    float prevError = std::numeric_limits<float>::max();
    for (size_t iter = 0; iter < 50; iter++) {
        for (size_t i = 0; i < nVertices; i++) {
            Matrix3f S = Matrix3f::Zero();
            for (auto [j, weight] : edgeWeights[i]) {
                S += weight * (vertices.col(i) - vertices.col(j)) * (newVertices.col(i) - newVertices.col(j)).transpose();
            }
            S += lambda * voronoiAreas[i] * normals.col(i) * this->m_bTargetNormals.col(i).transpose();
            auto SVD = S.jacobiSvd(ComputeFullU | ComputeFullV);
            Matrix3f R = SVD.matrixV() * SVD.matrixU().transpose();

            if (R.determinant() < 0) {
                // std::cout << "Negative:" << R.determinant() << std::endl;
                size_t smallestCol = 0;
                for (size_t j = 1; j < 3; j++) {
                    if (SVD.singularValues()[j] < SVD.singularValues()[smallestCol]) {
                        smallestCol = j;
                    }
                }
                Matrix3f U(SVD.matrixU());
                U.col(smallestCol) = -U.col(smallestCol);
                R = SVD.matrixV() * U.transpose();
                assert(R.determinant() >= 0);
            }
            Rs[i] = R;
        }

        b.setZero();
        for (size_t i = 0; i < nVertices; i++) {
            for (auto [j, weight] : edgeWeights[i]) {
                b.row(i) += weight / 2.f * ((Rs[i] + Rs[j]) * (vertices.col(i) - vertices.col(j)));
            }
        }

        prevVertices = newVertices;
        newVertices = solver->solve(b).transpose();

        float error = (newVertices - prevVertices).norm();
        // std::cout << error << std::endl;
        // threshold
        if (error < 1e-2) {
            break;
        }
        prevError = error;
    }
}

MatrixXf Analogy::computeTargetNormals()
{
    // n x 3
    Eigen::MatrixXf aPrimeNormals = this->m_aPrime.computeFaceNormals().transpose();
    // 3 x m
    Eigen::MatrixXf bNormals = this->m_b.computeVertexNormals();
    // n x m
    Eigen::MatrixXf similarities = aPrimeNormals * bNormals;
    // m x n
    for (size_t i = 0; i < bNormals.cols(); i++) {
        size_t maxIndex;
        similarities.col(i).maxCoeff(&maxIndex);
        bNormals.col(i) = aPrimeNormals.row(maxIndex);
    }
    return bNormals;
}

inline float triangleArea(const Eigen::Vector3f& v1, const Eigen::Vector3f& v2, const Eigen::Vector3f& v3) {
    return 0.5 * (v2 - v1).cross(v3 - v1).norm();
}

VectorXf Analogy::computeVoronoiAreas()
{
    MatrixXf& vertices = this->m_b.getVertices();
    MatrixXi& faces = this->m_b.getFaces();


    VectorXf voronoiAreas = VectorXf::Zero(vertices.cols());

    for (int i = 0; i < faces.cols(); i++) {
        // Get the indices of the vertices of the triangle
        int idx1 = faces(0, i);
        int idx2 = faces(1, i);
        int idx3 = faces(2, i);

        // Retrieve the vertices of the triangle
        Eigen::Vector3f v1 = vertices.col(idx1);
        Eigen::Vector3f v2 = vertices.col(idx2);
        Eigen::Vector3f v3 = vertices.col(idx3);

        // Compute the area of the triangle
        float area = triangleArea(v1, v2, v3);

        // Compute the angles of the triangle
        float angle1 = acos((v2 - v1).dot(v3 - v1) / ((v2 - v1).norm() * (v3 - v1).norm()));
        float angle2 = acos((v1 - v2).dot(v3 - v2) / ((v1 - v2).norm() * (v3 - v2).norm()));
        float angle3 = acos((v1 - v3).dot(v2 - v3) / ((v1 - v3).norm() * (v2 - v3).norm()));

        // Distribute the triangle's area to the vertices based on the angles
        voronoiAreas[idx1] += (area / 8.0) * (tan(angle2/2.0) + tan(angle3/2.0));
        voronoiAreas[idx2] += (area / 8.0) * (tan(angle1/2.0) + tan(angle3/2.0));
        voronoiAreas[idx3] += (area / 8.0) * (tan(angle1/2.0) + tan(angle2/2.0));
    }

    return voronoiAreas;
}

std::vector<std::unordered_map<size_t, float> > Analogy::computeEdgeWeights()
{
    MatrixXf& vertices = this->m_b.getVertices();
    MatrixXi& faces = this->m_b.getFaces();
    std::vector<std::unordered_map<size_t, float> > edgeWeights(vertices.cols());
    for (int j = 0; j < faces.cols(); j++) {
        VectorXi face = faces.col(j);
        for (int i = 0; i < 3; i++) {
            int v = face[i];
            int edgeA = face[(i + 1) % 3];
            int edgeB = face[(i + 2) % 3];
            Vector3f a = vertices.col(edgeA) - vertices.col(v);
            Vector3f b = vertices.col(edgeB) - vertices.col(v);

            float edgeWeight = std::abs(a.dot(b) / a.cross(b).norm() / 2.f);
            assert(edgeWeight >= 0);
            edgeWeights[edgeA].emplace(edgeB, 0.f);
            edgeWeights[edgeB].emplace(edgeA, 0.f);
            edgeWeights[edgeA][edgeB] += edgeWeight;
            edgeWeights[edgeB][edgeA] += edgeWeight;
        }
    }
    return edgeWeights;
}

std::unique_ptr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>>> Analogy::buildL(const std::vector<std::unordered_map<size_t, float> >& edgeWeights)
{
    MatrixXf& vertices = this->m_b.getVertices();
    SparseMatrix<float> L(vertices.cols(), vertices.cols());
    for (size_t i = 0; i < vertices.cols(); i++) {
        for (auto [j, weight] : edgeWeights[i]) {
            L.coeffRef(i, i) += weight;
            L.coeffRef(i, j) -= weight;
        }
    }

    MatrixXf dense = L.toDense();
    // precompute the decomposition
    std::unique_ptr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>>> solver = std::make_unique<Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>>>();
    solver->compute(L);
    // assert(solver->info()==Success);
    return solver;
}
