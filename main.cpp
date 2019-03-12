#define CGAL_EIGEN3_ENABLED

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/IO/read_xyz_points.h>

#include <vector>
#include <fstream>
#include <string>

#include "STL.h"

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Pwn;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_Mesh;

int main(void)
{
    STL_manager<Kernel> stl;

    std::vector<Pwn> points;
    std::ifstream stream("data/kitten.xyz");
    if (!stream ||
        !CGAL::read_xyz_points_and_normals(stream,
                                           std::back_inserter(points),
                                           CGAL::First_of_pair_property_map<Pwn>(),
                                           CGAL::Second_of_pair_property_map<Pwn>()))
    {
        std::cerr << "Error: cannot read file data/kitten.xyz" << std::endl;
        return EXIT_FAILURE;
    }
    Polyhedron output_mesh;

    double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(points.begin(), points.end(), CGAL::First_of_pair_property_map<Pwn>(), 6);
    if (CGAL::poisson_surface_reconstruction_delaunay(points.begin(),
                                                      points.end(),
                                                      CGAL::First_of_pair_property_map<Pwn>(),
                                                      CGAL::Second_of_pair_property_map<Pwn>(),
                                                      output_mesh, average_spacing))
    {
        Surface_Mesh smesh;
        CGAL::copy_face_graph(output_mesh, smesh);
        std::string fo = "out.stl";
        std::ofstream outfile(fo);
        CGAL::write_STL(smesh, outfile);
    }
    else
        return EXIT_FAILURE;
    return EXIT_SUCCESS;
}