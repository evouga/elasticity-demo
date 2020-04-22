#include "GooHook.h"
#include <igl/triangle/triangulate.h>

using namespace Eigen;

void GooHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputDouble("Timestep",  &params_.timeStep);        
        ImGui::InputDouble("Lame Alpha", &params_.lameAlpha);
        ImGui::InputDouble("Lame Beta", &params_.lameBeta);
        ImGui::InputDouble("Density", &params_.density);
    }
    if (ImGui::CollapsingHeader("Forces", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Checkbox("Gravity Enabled", &params_.gravityEnabled);
        ImGui::InputDouble("  Gravity g", &params_.gravityG);        
    }    
}

void GooHook::updateRenderGeometry()
{
    renderF = restF;
    int nverts = restV.rows();
    Eigen::Vector3d color(0.5, 0.5, 0.0);
    renderC.resize(nverts, 3);
    renderQ.resize(nverts, 3);
    for (int i = 0; i < nverts; i++)
    {
        for (int j = 0; j < 2; j++)
            renderQ(i, j) = Q[2 * i + j];
        renderQ(i, 2) = 0;
        renderC.row(i) = color;
    }    
}


void GooHook::initSimulation()
{
    time_ = 0;    

    Eigen::MatrixXd planarV(10, 2);
    planarV << -.5, .5,
        .5, .5,
        .5, -.5,
        -.5, -.5,
        -.5, .3,
        -.4, .3,
        -.4, -.4,
        0.4, -.4,
        .4,.4,
        -.5,.4;
    Eigen::MatrixXi planarE(10, 2);
    planarE << 0, 1,
        1, 2,
        2, 3,
        3, 4,
        4, 5,
        5, 6,
        6, 7,
        7, 8,
        8, 9,
        9, 0;


    std::string flags("za0.002q");

    Eigen::MatrixXi H(0,2);
    
    igl::triangle::triangulate(planarV, planarE, H, flags, restV, restF);

    int nverts = restV.rows();
    Q.resize(2 * nverts);
    for (int i = 0; i < nverts; i++)
        Q.segment<2>(2 * i) = restV.row(i);
    Qdot.resize(2 * nverts);
    Qdot.setZero();
}

void GooHook::tick()
{

}

bool GooHook::simulateOneStep()
{
    numericalIntegration(Q, Qdot);
    
    time_ += params_.timeStep;
    return false;
}

void GooHook::numericalIntegration(VectorXd &q, VectorXd &v)
{
    VectorXd F;
    SparseMatrix<double> Minv;

    computeMassInverse(Minv);
        
    q += params_.timeStep*v;
    computeForce(q, F);
    v += params_.timeStep*Minv*F;    
}

void GooHook::computeMassInverse(Eigen::SparseMatrix<double> &Minv)
{
    int ndofs = 2*restV.rows();

    Minv.resize(ndofs, ndofs);
    Minv.setZero();
    std::vector<Eigen::Triplet<double> > Minvcoeffs;

    Eigen::SparseMatrix<double> M;
    computeMass(M);
    for (int i = 0; i < ndofs; i++)
    {
        Minvcoeffs.push_back({ i, i, 1.0/M.coeff(i,i) });
    }
    
    Minv.setFromTriplets(Minvcoeffs.begin(), Minvcoeffs.end());
}

void GooHook::computeMass(Eigen::SparseMatrix<double> &M)
{
    int ndofs = 2*restV.rows();

    M.resize(ndofs, ndofs);
    M.setZero();
    std::vector<Eigen::Triplet<double> > Mcoeffs;

    int nfaces = restF.rows();
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector2d v0 = restV.row(restF(i, 0)).transpose();
        Eigen::Vector2d v1 = restV.row(restF(i, 1)).transpose();
        Eigen::Vector2d v2 = restV.row(restF(i, 2)).transpose();
        Eigen::Vector2d e1 = v1 - v0;
        Eigen::Vector2d e2 = v2 - v0;
        double area = 0.5 * std::fabs(e1[0] * e2[1] + e1[1] * e2[0]);
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                Mcoeffs.push_back({ 2 * restF(i, j) + k, 2 * restF(i, j) + k, params_.density * area / 3.0 });
            }
        }
    }

    M.setFromTriplets(Mcoeffs.begin(), Mcoeffs.end());
}

void GooHook::computeForce(const VectorXd &q, Eigen::VectorXd &F)
{
    F.resize(q.size());
    F.setZero();
    
    std::vector<Eigen::Triplet<double> > Hcoeffs;
    if(params_.gravityEnabled)
        processGravityForce(F);
    processElasticForce(q, F);

    F[0] = 0;
    F[1] = 0;
}

void GooHook::processGravityForce(VectorXd &F)
{
    int nverts = restV.rows();
    Eigen::SparseMatrix<double> M;
    computeMass(M);
    Eigen::VectorXd gF(2*nverts);
    gF.setZero();
    for(int i=0; i<nverts; i++)
    {
        gF[2*i+1] += params_.gravityG;        
    }
    F += M * gF;
}

Vector4d GooHook::dPhibar(const Eigen::VectorXd& q, int face)
{
    Eigen::Matrix2d Bbar;
    Bbar.col(0) = restV.row(restF(face, 1)).transpose() - restV.row(restF(face, 0)).transpose();
    Bbar.col(1) = restV.row(restF(face, 2)).transpose() - restV.row(restF(face, 0)).transpose();
    Eigen::Matrix2d Bbarinv = Bbar.inverse();
    Eigen::Vector2d du = Bbarinv.row(0);
    Eigen::Vector2d dv = Bbarinv.row(1);
    Eigen::Vector2d e1(1, 0);
    Eigen::Vector2d e2(0, 1);
    Eigen::Vector2d q0 = q.segment<2>(2 * restF(face, 0));
    Eigen::Vector2d q1 = q.segment<2>(2 * restF(face, 1));
    Eigen::Vector2d q2 = q.segment<2>(2 * restF(face, 2));
    Vector4d result;
    result[0] = e1.dot(q1 - q0) * du.dot(e1) + e1.dot(q2 - q0) * dv.dot(e1);
    result[1] = e1.dot(q1 - q0) * du.dot(e2) + e1.dot(q2 - q0) * dv.dot(e2);
    result[2] = e2.dot(q1 - q0) * du.dot(e1) + e2.dot(q2 - q0) * dv.dot(e1);
    result[3] = e2.dot(q1 - q0) * du.dot(e2) + e2.dot(q2 - q0) * dv.dot(e2);

    return result;
}

Eigen::Matrix<double, 4, 6> GooHook::ddPhibar(const Eigen::VectorXd& q, int face)
{
    Eigen::Matrix2d Bbar;
    Bbar.col(0) = restV.row(restF(face, 1)).transpose() - restV.row(restF(face, 0)).transpose();
    Bbar.col(1) = restV.row(restF(face, 2)).transpose() - restV.row(restF(face, 0)).transpose();
    Eigen::Matrix2d Bbarinv = Bbar.inverse();
    Eigen::Vector2d du = Bbarinv.row(0);
    Eigen::Vector2d dv = Bbarinv.row(1);
    Eigen::Vector2d e1(1, 0);
    Eigen::Vector2d e2(0, 1);
    Eigen::Vector2d q0 = q.segment<2>(2 * restF(face, 0));
    Eigen::Vector2d q1 = q.segment<2>(2 * restF(face, 1));
    Eigen::Vector2d q2 = q.segment<2>(2 * restF(face, 2));

    Matrix<double, 4, 6> result;
    result.setZero();
    result.block(0, 0, 1, 2) = -du.transpose() * e1 * e1.transpose() - dv.transpose() * e1 * e1.transpose();
    result.block(0, 2, 1, 2) = du.transpose() * e1 * e1.transpose();
    result.block(0, 4, 1, 2) = dv.transpose() * e1 * e1.transpose();
    result.block(1, 0, 1, 2) = -du.transpose() * e2 * e1.transpose() - dv.transpose() * e2 * e1.transpose();
    result.block(1, 2, 1, 2) = du.transpose() * e2 * e1.transpose();
    result.block(1, 4, 1, 2) = dv.transpose() * e2 * e1.transpose();
    result.block(2, 0, 1, 2) = -du.transpose() * e1 * e2.transpose() - dv.transpose() * e1 * e2.transpose();
    result.block(2, 2, 1, 2) = du.transpose() * e1 * e2.transpose();
    result.block(2, 4, 1, 2) = dv.transpose() * e1 * e2.transpose();
    result.block(3, 0, 1, 2) = -du.transpose() * e2 * e2.transpose() - dv.transpose() * e2 * e2.transpose();
    result.block(3, 2, 1, 2) = du.transpose() * e2 * e2.transpose();
    result.block(3, 4, 1, 2) = dv.transpose() * e2 * e2.transpose();
    return result;
}

Eigen::Vector4d GooHook::Sbar(const Eigen::VectorXd& q, int face)
{
    Vector4d dphibar = dPhibar(q, face);
    Vector4d result;
    result[0] = (dphibar[0] * dphibar[0] + dphibar[2] * dphibar[2] - 1) / 2;
    result[1] = (dphibar[0] * dphibar[1] + dphibar[2] * dphibar[3]) / 2;
    result[2] = (dphibar[0] * dphibar[1] + dphibar[2] * dphibar[3]) / 2;
    result[3] = (dphibar[1] * dphibar[1] + dphibar[3] * dphibar[3] - 1) / 2;
    return result;
}

Eigen::Matrix<double, 4, 6> GooHook::dSbar(const Eigen::VectorXd& q, int face)
{
    Vector4d dphibar = dPhibar(q, face);
    Matrix<double, 4, 6> ddphibar = ddPhibar(q, face);
    Eigen::Matrix4d frontmat;
    frontmat.setZero();
    frontmat(0, 0) = dphibar[0];
    frontmat(0, 2) = dphibar[2];
    frontmat(1, 0) = 0.5*dphibar[1];
    frontmat(1, 1) = 0.5*dphibar[0];
    frontmat(1, 2) = 0.5*dphibar[3];
    frontmat(1, 3) = 0.5*dphibar[2];
    frontmat(2, 0) = 0.5*dphibar[1];
    frontmat(2, 1) = 0.5*dphibar[0];
    frontmat(2, 2) = 0.5*dphibar[3];
    frontmat(2, 3) = 0.5*dphibar[2];
    frontmat(3, 1) = dphibar[1];
    frontmat(3, 3) = dphibar[3];
    return frontmat * ddphibar;
}

Eigen::Matrix<double, 1, 6> GooHook::dW(const Eigen::VectorXd& q, int face)
{
    Vector4d sbar = Sbar(q, face);
    Matrix<double, 4, 6> dsbar = dSbar(q, face);
    Eigen::Vector4d tracevec(1, 0, 0, 1);

    Matrix<double, 1, 6> result = params_.lameAlpha * tracevec.dot(sbar) * sbar.transpose() * dsbar;
    Vector4d transsbar(sbar[0], sbar[2], sbar[1], sbar[3]);
    result += 2.0 * params_.lameBeta * transsbar.transpose() * dsbar;
    return result;
}

void GooHook::processElasticForce(const Eigen::VectorXd &q, VectorXd& F)
{
    int nfaces = restF.rows();
    for (int i = 0; i < nfaces; i++)
    {
        Matrix<double, 1, 6> dw = dW(q, i);
        Eigen::Vector2d v0 = restV.row(restF(i, 0)).transpose();
        Eigen::Vector2d v1 = restV.row(restF(i, 1)).transpose();
        Eigen::Vector2d v2 = restV.row(restF(i, 2)).transpose();

        Eigen::Vector2d e0 = v1 - v0;
        Eigen::Vector2d e1 = v2 - v0;
        double area = 0.5 * std::fabs(e0[1] * e1[0] - e0[0] * e1[1]);
        Matrix<double, 1, 6> localF = -area * dw;

        F.segment<2>(2 * restF(i, 0)) += localF.segment<2>(0);
        F.segment<2>(2 * restF(i, 1)) += localF.segment<2>(2);
        F.segment<2>(2 * restF(i, 2)) += localF.segment<2>(4);
    }
}