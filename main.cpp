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

//NEW
#include <CGAL/trace.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/compute_average_spacing.h>

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

//NEW
typedef Kernel::FT FT;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef std::vector<Point_with_normal> PointList;
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;

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
/* OLD
    std::vector<Pwn> points;

    PointList points; //NEW
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
*/
// --------------------------------------------------------------------


    // Poisson options
    FT sm_angle = 10.0; // Min triangle angle in degrees.
    FT sm_radius = 0.5; // Max triangle size w.r.t. point set average spacing.
    FT sm_distance = 0.2; // Surface Approximation error w.r.t. point set average spacing.

    // Reads the point set file in points[].
    // Note: read_xyz_points_and_normals() requires an iterator over points
    // + property maps to access each point's position and normal.
    // The position property map can be omitted here as we use iterators over Point_3 elements.
    PointList points;
    std::ifstream stream("data/pointData.xyz");
    if (!stream ||
        !CGAL::read_xyz_points_and_normals(
                              stream,
                              std::back_inserter(points),
                              CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type())))
    {
      std::cerr << "Error: cannot read file data/kitten.xyz" << std::endl;
      return EXIT_FAILURE;
    }

    // Creates implicit function from the read points using the default solver.
    // Note: this method requires an iterator over points
    // + property maps to access each point's position and normal.
    // The position property map can be omitted here as we use iterators over Point_3 elements.
    Poisson_reconstruction_function function(points.begin(), points.end(),
                                             CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type()) );
    // Computes the Poisson indicator function f()
    // at each vertex of the triangulation.
    if ( ! function.compute_implicit_function() ) 
      return EXIT_FAILURE;

      // Computes average spacing
    FT average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(points.begin(), points.end(), 6 /* knn = 1 ring */);
    // Gets one point inside the implicit surface
    // and computes implicit function bounding sphere radius.
    Point inner_point = function.get_inner_point();
    Sphere bsphere = function.bounding_sphere();
    FT radius = std::sqrt(bsphere.squared_radius());
    // Defines the implicit surface: requires defining a
    // conservative bounding sphere centered at inner point.
    FT sm_sphere_radius = 5.0 * radius;
    FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance
    Surface_3 surface(function,
                      Sphere(inner_point,sm_sphere_radius*sm_sphere_radius),
                      sm_dichotomy_error/sm_sphere_radius);
    // Defines surface mesh generation criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
                                                        sm_radius*average_spacing,  // Max triangle size
                                                        sm_distance*average_spacing); // Approximation error
    // Generates surface mesh with manifold option
    STr tr; // 3D Delaunay triangulation for surface mesh generation
    C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
    CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
                            surface,                              // implicit surface
                            criteria,                             // meshing criteria
                            CGAL::Manifold_with_boundary_tag());  // require manifold mesh
    if(tr.number_of_vertices() == 0)
      return EXIT_FAILURE;

    Polyhedron output_mesh;
    CGAL::output_surface_facets_to_polyhedron(c2t3, output_mesh);
    Surface_Mesh smesh;
    CGAL::copy_face_graph(output_mesh, smesh);
    std::string fo = "out.stl";
    std::ofstream outfile(fo);
    CGAL::write_STL(smesh, outfile);
// --------------------------------------------------------------------





       /* 
    double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(points.begin(), points.end(), CGAL::First_of_pair_property_map<Pwn>(), 5);
    if (CGAL::poisson_surface_reconstruction_delaunay(points.begin(),
                                                      points.end(),
                                                      CGAL::First_of_pair_property_map<Pwn>(),
                                                      CGAL::Second_of_pair_property_map<Pwn>(),
                                                      output_mesh, 
                                                      average_spacing))
    {
   
        Surface_Mesh smesh;
        CGAL::copy_face_graph(output_mesh, smesh);
        std::string fo = "out.stl";
        std::ofstream outfile(fo);
        CGAL::write_STL(smesh, outfile);

    
        // Create domain
        Mesh_domain domain(output_mesh);
        
        // Mesh criteria (no cell_size set)
        Mesh_criteria criteria( facet_angle=25, 
                                facet_size=0.005, 
                                facet_distance=0.008,
                                cell_radius_edge_ratio=3,
                                cell_size=average_spacing
                                );
                                
        // Mesh generation
        //C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());
        C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
        // Output
        //std::ofstream medit_file("out_1.mesh");
        std::ofstream vtu_file("out_1.vtu");
        CGAL::output_to_vtu(vtu_file, c3t3);
        vtu_file.close();
    

    }
    else
        return EXIT_FAILURE;
    return EXIT_SUCCESS;
        */
}