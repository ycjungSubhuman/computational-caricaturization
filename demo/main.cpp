#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <Eigen/Dense>
#include <string>
#include <stdexcept>
#include <iostream>
#include <cstdlib>

#include "compcari.hpp"

using namespace compcari;

void print_help_and_die()
{
    std::cout << "Usage 1: <executable> ref <reference .obj file> <target .obj file> <beta (default 0.2)> <output .obj path (default: output.obj)>" << std::endl;
    std::cout << "Usage 2: <executable> noref <target .obj file> <gamma (default 0.2)> <output .obj path (default: output.obj)>" << std::endl;
    std::exit(1);
}

enum class Mode
{
    NOREF,
    REF,
};

void take_args(int argc, char **argv, Mode &mode, std::string &path_ref, std::string &path_target, std::string &path_output, double &factor)
{
    path_output = "output.obj";
    factor = 0.2;
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h")
        {
            print_help_and_die();
        }
    }

    if (argc <= 1)
    {
        print_help_and_die();
    }

    if (std::string(argv[1]) == "ref")
    {
        mode = Mode::REF;
    }
    else if (std::string(argv[1]) == "noref")
    {
        mode = Mode::NOREF;
    }
    else
    {
        print_help_and_die();
    }

    if (mode == Mode::NOREF)
    {
        if (argc < 3)
        {
            print_help_and_die();
        }
        path_target = argv[2];

        if (argc >= 4)
        {
            factor = std::atof(argv[3]);
        }

        if (argc >= 5)
        {
            path_output = argv[4];
        }
    }
    else if (mode == Mode::REF)
    {
        if (argc < 4)
        {
            print_help_and_die();
        }
        path_ref = argv[2];
        path_target = argv[3];

        if (argc >= 5)
        {
            factor = std::atof(argv[4]);
        }

        if (argc >= 6)
        {
            path_output = argv[5];
        }
    }
    else
    {
        throw std::runtime_error("Should not reach here");
    }
}

int main(int argc, char **argv)
{
    std::string path_ref, path_target, path_output;
    double factor;
    Mode mode;

    take_args(argc, argv, mode, path_ref, path_target, path_output, factor);

    Eigen::MatrixXd V_result;
    Eigen::MatrixXi F;

    if (mode == Mode::NOREF)
    {
        Eigen::MatrixXd V;
        igl::readOBJ(path_target, V, F);

        MeshPrep prep_mesh;
        preprocess_mesh(V, F, prep_mesh);

        ProblemPrep prep_problem;
        preprocess_problem_caricature_noref(V, prep_mesh, factor, prep_problem);

        solve_caricature_problem(prep_problem, V_result);
    }
    else if (mode == Mode::REF)
    {
        Eigen::MatrixXd V, V0;
        igl::readOBJ(path_ref, V0, F);
        igl::readOBJ(path_target, V, F);

        MeshPrep prep_mesh;
        preprocess_mesh(V, F, prep_mesh);

        ProblemPrep prep_problem;
        preprocess_problem_caricature_ref(V0, V, prep_mesh, factor, prep_problem);

        solve_caricature_problem(prep_problem, V_result);
    }

    igl::writeOBJ(path_output, V_result, F);

    return 0;
}