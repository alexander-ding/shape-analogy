#include "analogy.h"
#include <iostream>
#include <limits>

Analogy::Analogy(Mesh aPrime, Mesh b, bool sphereDeform)
    : m_aPrime(aPrime), m_b(b), m_bPrime(b), m_aPrimeCache(aPrime), m_bCache(b), m_tessellatedSphereNormals(m_aPrimeCache.computeFaceNormals().transpose()),
    m_computeTargetNormalsFunc(sphereDeform ? &Analogy::computeTargetNormalsMapping : &Analogy::computeTargetNormalsSnapping),
    m_prevFrameAPrimeVerts(aPrime.getVertices()), m_prevFrameBPrimeVerts(b.getVertices())
    {}

void Analogy::computeBPrime(float lambda)
{
    Eigen::MatrixXf& vertices = this->m_b.getVertices();
    Eigen::MatrixXf& newVertices = this->m_bPrime.getVertices();
    Eigen::MatrixXf prevVertices = newVertices;
    this->m_bTargetNormals = (this->*m_computeTargetNormalsFunc)();
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

/**
 * @brief barycentricInterpolation get normal at projected point p using barycentric interpolation
 * @return interpolated normal at point p if it lies inside the triangle, else the triangle face's normal
 */
Eigen::VectorXf barycentricInterpolation(Vector3f& p, Vector3f& vert1, Vector3f& vert2, Vector3f& vert3,
                                        Vector3f& n1, Vector3f& n2, Vector3f& n3, Eigen::VectorXf& faceNormal) {
    Vector3f v0 = vert2 - vert1;
    Vector3f v1 = vert3 - vert1;
    Vector3f v2 = p - vert1;
    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);
    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.f - v - w;
    // if projected point lies inside the triangle, use barycentric interpolation
    if ((v >= 0) && (w >= 0) && (v + w <= 1)) {
        Vector3f interpolated_normal = (u * n1 + v * n2 + w * n3).normalized();
        return interpolated_normal;
    } else { // else just approximate with the face normal
        return faceNormal;
    }
}

/**
 * @brief projectAndInterpolate orthogonal projection onto plane defined by triangle face, then barycentric interpolation
 * @return interpolated normal at point projected onto triangle face plane
 */
Eigen::VectorXf projectAndInterpolate(Eigen::VectorXf& p, Eigen::Vector3i& indices, Eigen::MatrixXf& vertices,
                                Eigen::MatrixXf& vertNormals, Eigen::VectorXf& faceNormal) {
    Eigen::Vector3f v0 = vertices.col(indices[0]);
    Eigen::Vector3f v1 = vertices.col(indices[1]);
    Eigen::Vector3f v2 = vertices.col(indices[2]);

    // Calculate the normal vector of the plane
    Vector3f n = -((v1 - v0).cross(v2 - v0).normalized());
    // vector from point on plane to point on our unit sphere
    Vector3f P = p - v0;
    // project the vector onto plane normal
    Vector3f proj = (P.dot(n) / n.squaredNorm()) * n;
    // project the point onto the plane
    Vector3f projectedPoint = p - proj;
    // vert normals
    Vector3f n0 = vertNormals.col(indices[0]);
    Vector3f n1 = vertNormals.col(indices[1]);
    Vector3f n2 = vertNormals.col(indices[2]);

    return barycentricInterpolation(projectedPoint, v0, v1, v2, n0, n1, n2, faceNormal);
}

MatrixXf Analogy::computeTargetNormalsMapping() // between B and initial sphere
{
    // n x 3
    Eigen::MatrixXf aPrimeNormals = m_tessellatedSphereNormals; // cached
    // 3 x m
    Eigen::MatrixXf bNormals = this->m_b.computeVertexNormals();
    // n x m
    Eigen::MatrixXf similarities = aPrimeNormals * bNormals;
    // 1 x m
    Eigen::MatrixXf aPrimeDeformedNormals = this->m_aPrime.computeFaceNormals().transpose();

    Eigen::MatrixXi faces = this->m_aPrimeCache.getFaces();
    Eigen::MatrixXf vertNormals = this->m_aPrime.computeVertexNormals();
    Eigen::MatrixXf vertices = this->m_aPrimeCache.getVertices();

    // m x n
    for (size_t i = 0; i < bNormals.cols(); i++) {
        size_t maxIndex;
        similarities.col(i).maxCoeff(&maxIndex);
        Vector3i indices = faces.col(maxIndex); // closest face vert
        VectorXf p = bNormals.col(i); // point on unit sphere
        VectorXf faceNormal = aPrimeDeformedNormals.row(maxIndex); // use if projection is not inside triangle
        bNormals.col(i) = projectAndInterpolate(p, indices, vertices, vertNormals, faceNormal);
    }
    return bNormals;
}

MatrixXf Analogy::computeTargetNormalsSnapping()
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

void Analogy::undo() {
    m_aPrime.setVertices(m_prevFrameAPrimeVerts);
    m_bPrime.setVertices(m_prevFrameBPrimeVerts);
}
