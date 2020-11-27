#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace compcari
{
    /**
     * Mesh pre-processing results
     */
    typedef struct MeshPrep
    {
        Eigen::MatrixXi F;                 // Face vertex indices
        Eigen::SparseMatrix<double> L;     // Laplacian operator
        Eigen::SparseMatrix<double> Div;   // Divergence operator
        Eigen::SparseMatrix<double> G;     // Gradient operator
        Eigen::MatrixXd dir_min_principal; // Principal minimum curvature directions
                                           // for each vertex
        Eigen::MatrixXd dir_max_principal; // Principal maximum curvature directions
                                           // for each vertex
        Eigen::VectorXd C_min;             // Minimum curvature for each vertex
        Eigen::VectorXd C_max;             // Maximum curcature for each vertex
        Eigen::VectorXd C_mean;            // Mean curvature for each vertex
        Eigen::VectorXd C_gaussian;        // Gaussian curvature for each vertex
    } MeshPrep;

    using Decomposition = Eigen::SparseLU<Eigen::SparseMatrix<double>>;

    // Problem pre-precessing results.
    typedef struct ProblemPrep
    {
        Eigen::SparseMatrix<double> A;
        Eigen::VectorXd left;           // Leftmost column of Laplacian without its first element
        Decomposition dec;              // Matrix decomposition of the problem
        std::vector<Eigen::VectorXd> b; // Differential coordinate

        Eigen::Vector3d x0; // fixed position of 0-th vertex
    } ProblemPrep;

    /**
     * Pre-process data associated to the mesh for other operations
     *
     * Input
     * @param V                 vertex positions. (#Vx3)
     * @param F                 face vertex indices. (#Fx3)
     *
     * Output
     * @param prep_mesh         preprocessing result.
     */
    void preprocess_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi F,
                         MeshPrep &prep_mesh);
    /**
     * Pre-process the caricature problem without reference shape.
     *
     * Call this function along with other problem-specific preprocessing. i.e.
     * preprocess_problem_laplacian or preprocess_problem_caricature
     * 
     * Note)
     *  Because this is a gradient-domain scaling, when given extreme parameters, 
     *  the scale of the mesh is usually not preserved.
     * 
     * Input
     * @param V                 vertex positions. (#Vx3)
     * @param prep_mesh         the result of 'preprocess_mesh'
     *
     * Output
     * @param prep_problem      preprocessing result.
     */
    void preprocess_problem_caricature_noref(const Eigen::MatrixXd &V,
                                             const MeshPrep &prep_mesh,
                                             const float gamma,
                                             ProblemPrep &prep_problem);

    /**
     * Pre-process the caricature problem without reference shape.
     *
     * Call this function along with other problem-specific preprocessing. i.e.
     * preprocess_problem_laplacian or preprocess_problem_caricature
     * 
     * Note)
     *  When given the reference mesh, the scale of the mesh is better preserved.
     *  I recommend using this function.
     * 
     * Input
     * @param V_ref             vertex positions of a reference shape (#Vx3)
     * @param V                 vertex positions. (#Vx3)
     * @param prep_mesh         the result of 'preprocess_mesh'
     *
     * Output
     * @param prep_problem      preprocessing result.
     */
    void preprocess_problem_caricature_ref(const Eigen::MatrixXd &V_ref,
                                           const Eigen::MatrixXd &V,
                                           const MeshPrep &prep_mesh,
                                           const float beta,
                                           ProblemPrep &prep_problem);

    /**
     * Solve the problem
     *
     * Input
     * @param prep_problem      preprocessing result after
     *                          'preprecess_problem_common' and a problem-specific
     *                          preprocessing.
     *
     * Output
     * @para V_result           output vertex position
     */
    void solve_caricature_problem(const ProblemPrep &prep_problem, Eigen::MatrixXd &V_result);
} // namespace compcari