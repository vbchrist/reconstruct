// std
#include <iostream>
#include <fstream>
#include <string>

// CGAL
#include <CGAL/IO/STL_reader.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/license/Polyhedron.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

// local
#include "STL_writer.h"

class STL_manager
{
    typedef CGAL::Simple_cartesian<double> Kernel;
    typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_Mesh;

  public:
    bool import_STL(const std::string &file, Surface_Mesh &mesh)
    {
        typedef CGAL::Simple_cartesian<double> Kernel;
        typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_Mesh;

        std::vector<CGAL::cpp11::array<double, 3>> points;
        std::vector<CGAL::cpp11::array<int, 3>> triangles;

        std::ifstream file_stream(file, std::ios::binary);

        std::cout << "Reading STL from " << file << " ...\n";
        CGAL::read_STL(file_stream, points, triangles);
        if (CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles))
        {
            CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, triangles, mesh);
            if (!mesh.is_valid() || mesh.is_empty())
            {
                std::cout << "Invalid or empty mesh.\n";
                std::cout << "Cannot proceed.\n";
                return false;
            }
            else
            {
                std::cout << "Surface mesh sucessfully read.\n";
                std::cout << "  Points: " << points.size() << "\n";
                std::cout << "  Triangles: " << triangles.size() << "\n";
                std::cout << "  Area: " << CGAL::Polygon_mesh_processing::area(mesh)<< "\n";
                std::cout << "  Volume: " << CGAL::Polygon_mesh_processing::volume(mesh) << "\n";                
                std::cout << "Mesh read successfully.\n\n";
                return true;
            }
        }
        return false;
    }
    bool export_STL(const std::string &file, const Surface_Mesh &mesh)
    {
        std::cout << "Writing new mesh to " << file << " ...\n";
        std::ofstream outfile(file);
        CGAL::write_STL(mesh, outfile);
        outfile.close();
        std::cout << "Mesh write complete.\n";
    }
};