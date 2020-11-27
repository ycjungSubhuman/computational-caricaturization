#include "compcari.hpp"

#include <cassert>

#include <igl/doublearea.h>
#include <igl/grad.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/principal_curvature.h>
#include <igl/average_onto_faces.h>

namespace compcari
{
    namespace
    {
        void preprocess_problem_common(
            const Eigen::MatrixXd &V,
            const MeshPrep &prep_mesh,
            ProblemPrep &prep_problem)
        {
            prep_problem.b.resize(3);
            const auto &L = prep_mesh.L;
            prep_problem.A = L.block(1, 1, L.rows() - 1, L.cols() - 1);
            prep_problem.left = L.block(1, 0, L.rows() - 1, 1);
            prep_problem.dec.compute(prep_problem.A);
            if (prep_problem.dec.info() != Eigen::Success)
            {
                throw std::runtime_error("Solve failed");
            }
            prep_problem.x0 = V.row(0).transpose();
        }

        void divergence(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                        Eigen::SparseMatrix<double> &Div)
        {
            // Div implementation taken from libigl's heat_geodesics.cpp implementation
            Eigen::VectorXd dblA;
            igl::doublearea(V, F, dblA);
            Eigen::SparseMatrix<double> G;
            igl::grad(V, F, G);
            Div = -0.5 * G.transpose() * dblA.colwise().replicate(3).asDiagonal();
        }

        void preprocess_problem_edit_grad(
            const Eigen::MatrixXd &V,
            const Eigen::VectorXd &scale,
            const MeshPrep &prep_mesh,
            ProblemPrep &prep_problem)
        {
            Eigen::MatrixXd lap[3];
            for (int k = 0; k < 3; k++)
            {
                lap[k] = prep_mesh.L * V.col(k);
            }

            for (int k = 0; k < 3; k++)
            {
                Eigen::VectorXd GV = prep_mesh.G * V.col(k);
                Eigen::VectorXd sGV(GV.rows());
                assert(GV.rows() % 3 == 0);
                assert(GV.rows() / 3 == scale.rows());
                for (Eigen::Index i = 0; i < GV.rows() / 3; i++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        sGV(GV.rows() / 3 * l + i) = scale(i) * GV(GV.rows() / 3 * l + i);
                    }
                }
                prep_problem.b[k] = prep_mesh.Div * sGV;
                prep_problem.b[k] *= lap[k].norm() / prep_problem.b[k].norm();
            }
        }
    } // namespace

    void preprocess_mesh(
        const Eigen::MatrixXd &V, const Eigen::MatrixXi F,
        MeshPrep &prep_mesh)
    {
        prep_mesh.F = F;
        igl::cotmatrix(V, F, prep_mesh.L);
        divergence(V, F, prep_mesh.Div);
        igl::grad(V, F, prep_mesh.G);
        igl::principal_curvature(V, F, prep_mesh.dir_max_principal,
                                 prep_mesh.dir_min_principal, prep_mesh.C_max,
                                 prep_mesh.C_min);
        prep_mesh.C_mean = (prep_mesh.C_max + prep_mesh.C_min) / 2.0;
        prep_mesh.C_gaussian =
            (prep_mesh.C_max.array() * prep_mesh.C_min.array()).matrix();
    }

    void preprocess_problem_caricature_noref(
        const Eigen::MatrixXd &V,
        const MeshPrep &prep_mesh,
        const float gamma,
        ProblemPrep &prep_problem)
    {
        preprocess_problem_common(V, prep_mesh, prep_problem);
        Eigen::VectorXd K;
        igl::average_onto_faces(prep_mesh.F, prep_mesh.C_gaussian, K);
        Eigen::VectorXd scale(K.rows());
        for (Eigen::Index i = 0; i < scale.rows(); i++)
        {
            scale(i) = std::pow(std::abs(K(i)), gamma);
        }
        preprocess_problem_edit_grad(V, scale, prep_mesh, prep_problem);
    }

    void preprocess_problem_caricature_ref(
        const Eigen::MatrixXd &V_ref,
        const Eigen::MatrixXd &V,
        const MeshPrep &prep_mesh,
        const float beta,
        ProblemPrep &prep_problem)
    {
        preprocess_problem_common(V, prep_mesh, prep_problem);
        Eigen::VectorXd K;
        igl::average_onto_faces(prep_mesh.F, prep_mesh.C_gaussian, K);

        Eigen::VectorXd dblA, dblA_ref;
        igl::doublearea(V, prep_mesh.F, dblA);
        igl::doublearea(V_ref, prep_mesh.F, dblA_ref);
        assert(dblA.rows() == K.rows());

        Eigen::VectorXd scale(K.rows());
        for (Eigen::Index i = 0; i < scale.rows(); i++)
        {
            scale(i) =
                std::pow(std::abs(K(i)), beta * -std::log(dblA(i) / dblA_ref(i)));
        }
        preprocess_problem_edit_grad(V, scale, prep_mesh, prep_problem);
    }

    void solve_caricature_problem(
        const ProblemPrep &prep_problem, Eigen::MatrixXd &V_result)
    {
        const int num_vertex = prep_problem.b.at(0).rows();
        V_result.resize(num_vertex, 3);

        const auto &b = prep_problem.b;
        const auto &left = prep_problem.left;
        const auto &x0 = prep_problem.x0;
        const auto &dec = prep_problem.dec;

        for (int k = 0; k < 3; k++)
        {
            Eigen::VectorXd bsub = b.at(k).segment(1, num_vertex - 1);
            V_result(0, k) = x0(k);
            V_result.block(1, k, num_vertex - 1, 1) = dec.solve(bsub - x0(k) * left);
            if (dec.info() != Eigen::Success)
            {
                throw std::runtime_error("Solve failed in dimension " +
                                         std::to_string(k));
            }
        }
    }
} // namespace compcari