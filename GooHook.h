#include "PhysicsHook.h"
#include "SceneObjects.h"
#include <deque>
#include "SimParameters.h"
#include <Eigen/Sparse>
#include <Eigen/StdVector>

class GooHook : public PhysicsHook
{
public:
    GooHook() : PhysicsHook() {}

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu);

    virtual void initSimulation();

    virtual void mouseClicked(double x, double y, int button)
    {

    }

    virtual void updateRenderGeometry();

    virtual void tick();

    virtual bool simulateOneStep();

    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
    {
        viewer.data().clear();
        viewer.data().set_mesh(renderQ, renderF);
        viewer.data().set_colors(renderC);
    }

private:
    SimParameters params_;
    double time_;
    

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd renderC;

    Eigen::MatrixXd restV;
    Eigen::MatrixXi restF;

    Eigen::VectorXd Q;
    Eigen::VectorXd Qdot;

    void computeMass(Eigen::SparseMatrix<double> &M);
    void computeMassInverse(Eigen::SparseMatrix<double> &Minv);
    void numericalIntegration(Eigen::VectorXd &q, Eigen::VectorXd &v);

    void computeForce(const Eigen::VectorXd &q, Eigen::VectorXd &F);
    void processGravityForce(Eigen::VectorXd &F);

    Eigen::Vector4d dPhibar(const Eigen::VectorXd& q, int face);
    Eigen::Matrix<double, 4, 6> ddPhibar(const Eigen::VectorXd& q, int face);

    Eigen::Vector4d Sbar(const Eigen::VectorXd& q, int face);
    Eigen::Matrix<double, 4, 6> dSbar(const Eigen::VectorXd& q, int face);

    Eigen::Matrix<double, 1, 6> dW(const Eigen::VectorXd& q, int face);

    void processElasticForce(const Eigen::VectorXd &q, Eigen::VectorXd &F);    
};
