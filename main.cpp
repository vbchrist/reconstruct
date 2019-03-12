#define CGAL_EIGEN3_ENABLED

//http://cgal-discuss.949826.n4.nabble.com/Getting-error-quot-BOOST-PARAMETER-MAX-ARITY-must-be-at-least-12-for-CGAL-Mesh-3-quot-in-CGAL-4-9-td4662304.html
#define BOOST_PARAMETER_MAX_ARITY 12

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/IO/read_xyz_points.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
//#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

// IO
//#include <CGAL/IO/Polyhedron_iostream.h>

#include <vector>
#include <fstream>
#include <string>

#include "STL.h"
#include "output_to_vtu.h"



// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Pwn;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_Mesh;

// Domain
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Kernel> Mesh_domain;
//#ifdef CGAL_CONCURRENT_MESH_3
//typedef CGAL::Parallel_tag Concurrency_tag;
//#else
typedef CGAL::Sequential_tag Concurrency_tag;
//#endif
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

int main(void)
{
    // To avoid verbose function and named parameters call
    using namespace CGAL::parameters;

    STL_manager<Kernel> stl;

    std::vector<Pwn> points;
    std::ifstream stream("data/pointData.xyz");
    if (!stream ||
        !CGAL::read_xyz_points_and_normals(stream,
                                           std::back_inserter(points),
                                           CGAL::First_of_pair_property_map<Pwn>(),
                                           CGAL::Second_of_pair_property_map<Pwn>()))
    {
        std::cerr << "Error: cannot read file." << std::endl;
        return EXIT_FAILURE;
    }
    Polyhedron output_mesh;

    double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(points.begin(), points.end(), CGAL::First_of_pair_property_map<Pwn>(), 5);
    if (CGAL::poisson_surface_reconstruction_delaunay(points.begin(),
                                                      points.end(),
                                                      CGAL::First_of_pair_property_map<Pwn>(),
                                                      CGAL::Second_of_pair_property_map<Pwn>(),
                                                      output_mesh, 
                                                      average_spacing))
    {
        /*
        Surface_Mesh smesh;
        CGAL::copy_face_graph(output_mesh, smesh);
        std::string fo = "out.stl";
        std::ofstream outfile(fo);
        CGAL::write_STL(smesh, outfile);
*/
        //http://cgal-discuss.949826.n4.nabble.com/CGAL-ERROR-assertion-violation-td4661475.html
        //CGAL::Polygon_mesh_processing::triangulate_faces(output_mesh);

        // Create domain
        Mesh_domain domain(output_mesh);
        
        // Mesh criteria (no cell_size set)
        Mesh_criteria criteria( facet_angle=25, 
                                facet_size=0.005, 
                                facet_distance=0.008,
                                cell_radius_edge_ratio=3
                                );
                                
        // Mesh generation
        C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

        // Output
        //std::ofstream medit_file("out_1.mesh");
        std::ofstream vtu_file("out_1.vtu");
        CGAL::output_to_vtu(vtu_file, c3t3);
        vtu_file.close();
/*
        // Set tetrahedron size (keep cell_radius_edge_ratio), ignore facets
        Mesh_criteria new_criteria(cell_radius_edge_ratio=3, cell_size=0.003);

        // Mesh refinement
        CGAL::refine_mesh_3(c3t3, domain, new_criteria);

        // Output
        medit_file.open("out_2.mesh");
        c3t3.output_to_medit(medit_file);
        */
    }
    else
        return EXIT_FAILURE;
    return EXIT_SUCCESS;
}