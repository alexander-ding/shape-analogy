#ifndef ANALOGY_H
#define ANALOGY_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <vector>
#include <unordered_map>

#include "mesh.h"

class Analogy
{
public:
    Analogy(Mesh aPrime, Mesh b,  bool sphereDeform);

    Mesh& getAPrime() {return this->m_aPrime;}
    void setAPrime(const Mesh& aPrime) {this->m_aPrime = aPrime;}

    Mesh& getB() {return this->m_b;}
    Mesh& getBPrime() {return this->m_bPrime;}
    Eigen::MatrixXf& getBTargetNormals() { return this->m_bTargetNormals; }
    void computeBPrime(float lambda = 100);

    void setBPrime(Mesh mesh) { this->m_bPrime = mesh; }
    Mesh& getAPrimeCache() { return this->m_aPrimeCache; }

private:
    Eigen::MatrixXf computeTargetNormalsMapping();
    Eigen::MatrixXf computeTargetNormalsSnapping();
    Eigen::VectorXf computeVoronoiAreas();
    std::vector<std::unordered_map<size_t, float>> computeEdgeWeights();
    std::unique_ptr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>>> buildL(const std::vector<std::unordered_map<size_t, float> >& edgeWeights);

    Eigen::MatrixXf m_bTargetNormals;
    Mesh m_bPrime;
    Mesh m_aPrime;
    Mesh m_b;
    Mesh m_aPrimeCache;
    Mesh m_bCache;
    Eigen::MatrixXf m_tessellatedSphereNormals;
    Eigen::MatrixXf (Analogy::*m_computeTargetNormalsFunc)();
};

#endif // ANALOGY_H
