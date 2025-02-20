#include <iostream>
#include <filesystem>
#include <sys/stat.h>
#include <Types.h>
#include <CGAL/draw_polygon_2.h>
#include <Random.h>
#include <stdlib.h>
#include <set>
#include <vector>
#include <cmath>
#include <CGAL/draw_surface_mesh.h>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>  // For color conversion


#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random.h>

#include <Partition.h>
#include <IO.h>

#define PI 3.14159265358979323846

typedef CGAL::Creator_uniform_2<int, CGL::Point2>        Creator;
typedef CGAL::Random_points_in_square_2<CGL::Point2, Creator> Point_generator;


typedef CGL::Mesh Mesh;
typedef CGL::Point3 Point3;
typedef CGL::Point2 Point2;

typedef CGL::Mesh::Vertex_index Vertex_index;
typedef CGL::Mesh::Halfedge_index Halfedge_index;
typedef CGL::Mesh::Face_index Face_index;
typedef CGL::Mesh::Edge_index Edge_index;



float getDeterminant(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3);
float getDeterminant(std::vector<CGL::Point2> pts);
double getDeterminant4x4(double matrix[4][4], int n);

bool getOrientationTriangle(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3);
//Returns CW if collinear:
bool getOrientationTriangle(CGL::Point3 p1, CGL::Point3 p2,CGL::Point3 p3);
float getTriangleArea(std::vector<CGL::Point2> pts);


bool isInsideTriangle(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3, CGL::Point2 point);
bool isInsideTriangle(Point3 v1, Point3 v2, Point3 v3, Point3 pt) ;
bool isOnTriangle(Point3 v1, Point3 v2, Point3 v3, Point3 pt);

std::vector<CGL::Mesh::Vertex_index> get_face_vertices(CGL::Mesh::Face_index face, CGL::Mesh polygon);
std::vector<CGL::Mesh::Halfedge_index> get_face_halfedges(CGL::Mesh::Face_index face, CGL::Mesh polygon);

//Delaunay
bool inCircle(CGL::Point3 a, CGL::Point3 b,CGL::Point3 c,CGL::Point3 d);
std::vector<CGL::Point2>  createSuperTriangle(const std::vector<CGL::Point2>& points);
CGL::Mesh delaunay_triangulation(std::vector<CGL::Point2> mesh_vertices);
void insert_point(CGL::Mesh& mesh, CGL::Point2 vi);
std::vector<Halfedge_index> insert_triangles(Halfedge_index halfedge, Vertex_index new_vertex, CGL::Mesh& mesh);
std::vector<Halfedge_index> insert_point_on_edge(Halfedge_index halfedge, Vertex_index new_vertex, CGL::Mesh& mesh) ;
void legalize_edge(Halfedge_index hedge, CGL::Mesh& mesh);
void remove_super_triangle(CGL::Mesh& mesh, const std::vector<Vertex_index>& super_indexes);
void add_border(CGL::Mesh& mesh);

std::vector<Point2> select_mesh_points(const cv::Mat& img_pgm, unsigned long long amount = 1500, bool add_random_points = false);
void addBorderPoints(const cv::Mat& img_pgm, std::vector<Point2>& point_vec, int border_spacing = 16);

void computeBarycentrics(const Point3& p, const Point3& v0, const Point3& v1, const Point3& v2,
    double &lambda0, double &lambda1, double &lambda2);
uchar bilinearInterpolate(const cv::Mat& img, float x, float y);
cv::Mat rasterizeMesh(const Mesh& mesh, const cv::Mat& inputImg);

cv::Mat convert_to_pgm(const std::string& image_path, const std::string& output_path);

using namespace cv;


int main(int argc, char* argv[]) {
    try {
        // Check that an image path and pixel count were provided.
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " <image_path> <pixel_count>" << std::endl;
            return EXIT_FAILURE;
        }

        // Retrieve command-line arguments.
        std::string image_path = argv[1];
        unsigned long long pixel_amount = 0;
        try {
            pixel_amount = std::stoi(argv[2]);
        } catch (const std::exception& e) {
            std::cerr << "Invalid pixel count provided: " << argv[2] << std::endl;
            return EXIT_FAILURE;
        }

        // Define the output path for the PGM file.
        std::string output_path = "./original_converted.pgm";

        // Convert the input image to PGM.
        cv::Mat img_pgm = convert_to_pgm(image_path, output_path);
        if (img_pgm.empty()) {
            std::cout << "Image conversion failed." << std::endl;
            return EXIT_FAILURE;
        }
        std::cout << "Image conversion successful." << std::endl;

        std::vector<Point2> point_vec = select_mesh_points(img_pgm,pixel_amount,true);

        // Compute the Delaunay triangulation to obtain a mesh.
        Mesh mesh = delaunay_triangulation(point_vec);

        // Rasterize the mesh: for each triangle, interpolate colors from the input image.
        cv::Mat rasterized = rasterizeMesh(mesh, img_pgm);

        // Display the original and rasterized images.
        cv::imshow("Original Image", img_pgm);
        cv::imshow("Rasterized Mesh", rasterized);
        cv::waitKey(0);

        if (!cv::imwrite("rasterized.jpg", rasterized)) {
            std::cerr << "Failed to save the rasterized image." << std::endl;
            return EXIT_FAILURE;
        } else {
            std::cout << "Rasterized image saved successfully as 'rasterized.jpg'." << std::endl;
        }

    } catch (const std::exception& e) {
        std::cout << "Exception in main: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

std::vector<Point2> createSuperTriangle(const std::vector<Point2>& points) {
    // Find bounding box
    double minX = INFINITY;
    double maxX = -INFINITY;
    double minY = INFINITY;
    double maxY = -INFINITY;
    
    for (auto p = points.begin(); p != points.end(); ++p) {
        minX = std::min(minX, p->x());
        maxX = std::max(maxX, p->x());
        minY = std::min(minY, p->y());
        maxY = std::max(maxY, p->y());
    }
    
    // Calculate bounding box dimensions
    double dx = maxX - minX;
    double dy = maxY - minY;
    double dmax = std::max(dx, dy);
    double midx = (minX + maxX) / 2.0;
    double midy = (minY + maxY) / 2.0;
    
    // Create super triangle vertices
    std::vector<CGL::Point2>  superTriangle;

    Point2 p1(midx - 2 * dmax, midy - dmax);
    Point2 p2(midx + 2 * dmax, midy - dmax);
    Point2 p3(midx, midy + 2 * dmax);

    //CCW order
    superTriangle.push_back(p1);
    superTriangle.push_back(p2);
    superTriangle.push_back(p3);

    return superTriangle;
}

CGL::Mesh delaunay_triangulation(std::vector<CGL::Point2> vertices){

    CGL::Mesh delaunay_mesh;

    std::vector<CGL::Mesh::Vertex_index> super_indexes;
    std::vector<CGL::Point2> super_triangle = createSuperTriangle(vertices);

    // Add vertices to the mesh
    for (auto p = super_triangle.begin(); p != super_triangle.end(); ++p) {
        CGL::Point3 p3(p->x(), p->y(), 0.0);  // Convert 2D point to 3D 
        CGL::Mesh::Vertex_index vi = delaunay_mesh.add_vertex(p3);
        super_indexes.push_back(vi);
    }

    CGAL::Euler::add_face(super_indexes,delaunay_mesh);
    if(!delaunay_mesh.is_valid(false))
        printf("deu errado");

    // Gets the halfedges of the super triangle
    std::vector<Halfedge_index> super_hedges;
    for(Vertex_index vertex : super_indexes){
        // For some reason, delaunay_mesh.halfedge(vertex) returns the halfedge in the CW order
        super_hedges.push_back(delaunay_mesh.opposite(delaunay_mesh.halfedge(vertex))); 
    }


    unsigned i = 0;
    for (auto p = vertices.begin(); p != vertices.end(); ++p, i++) {
        insert_point(delaunay_mesh,*p);
        if(i%500 == 0)
            std::cout << "i = " << i << std::endl;
    }

    std::cout << "Ended inserting points" << std::endl;

    // CGL::print_mesh(delaunay_mesh);
    for(i = 0;i < 5;i++){
        for(Halfedge_index h : delaunay_mesh.halfedges())
            legalize_edge(h,delaunay_mesh);
    }

    remove_super_triangle(delaunay_mesh,super_indexes);

    add_border(delaunay_mesh);
    add_border(delaunay_mesh);
    add_border(delaunay_mesh);

    delaunay_mesh.collect_garbage();

    for(i = 0;i < 5;i++){
        for(Halfedge_index h : delaunay_mesh.halfedges())
            legalize_edge(h,delaunay_mesh);
    }

    if (!delaunay_mesh.is_valid(false))
        printf("Deu ruim");

    CGAL::draw(delaunay_mesh);

    return delaunay_mesh;
}

void insert_point(CGL::Mesh& mesh, CGL::Point2 vi){

    Point3 p3(vi.x(), vi.y(), 0.0);  // Convert 2D point to 3D 
    Vertex_index new_vertice = mesh.add_vertex(p3);

    Halfedge_index point_halfedge;
    //Get face the point is inside
    Face_index boundary_triangle;
    bool isInside = false;
    for(Face_index f : mesh.faces()){
        Halfedge_index hedge = mesh.halfedge(f);
        Halfedge_index next_he = mesh.next(hedge);
        Halfedge_index next_next_he = mesh.next(next_he);

        // Get vertices
        Point3 t1 = mesh.point(mesh.target(hedge));
        Point3 t2 = mesh.point(mesh.target(next_he));
        Point3 t3 = mesh.point(mesh.target(next_next_he));


        if(isInsideTriangle(t1, t2, t3, p3)){
            boundary_triangle = f;
            isInside = true;
            break;
        }
        else if(isOnTriangle(t1, t2, t3, p3)){

            if(CGL::Collinear(t1,t2,p3))
                point_halfedge = next_he;
            else if(CGL::Collinear(t2,t3,p3))
                point_halfedge = next_next_he;
            else if(CGL::Collinear(t3,t1,p3))
                point_halfedge = hedge;

            break;
        }
    }

    if(!isInside){
        std::vector<Halfedge_index> new_hedges = insert_point_on_edge(point_halfedge,new_vertice,mesh);

        legalize_edge(mesh.prev(new_hedges[0]), mesh);
        legalize_edge(mesh.prev(new_hedges[1]), mesh);
        legalize_edge(mesh.prev(new_hedges[2]), mesh);
        legalize_edge(mesh.prev(new_hedges[3]), mesh);
    }
    else{
        Halfedge_index halfedge = mesh.halfedge(boundary_triangle);

        // std::vector<Vertex_index> face_vertexes = get_face_vertices(boundary_triangle,mesh);
        //add_face do cgal n funciona pq ele recebe vertice como parametro, 
        //e os vertices aponta para as twins

        // Legalize the three new edges
        std::vector<Halfedge_index> new_hedges = insert_triangles(halfedge,new_vertice,mesh);

        legalize_edge(mesh.prev(new_hedges[0]), mesh);
        legalize_edge(mesh.prev(new_hedges[1]), mesh);
        legalize_edge(mesh.prev(new_hedges[2]), mesh);
    }

    // CGL::print_mesh(mesh);
    // CGAL::draw(mesh);

}
std::vector<Halfedge_index> insert_triangles(Halfedge_index halfedge, Vertex_index new_vertex, CGL::Mesh& mesh) {
    // // Debug output
    // std::cout << "Starting insert_triangles operation" << std::endl;
    // std::cout << "Initial halfedge: " << halfedge << ", New vertex: " << new_vertex << std::endl;

    try {
        // Store original face and its halfedges
        Face_index old_face = mesh.face(halfedge);
        Halfedge_index he1 = halfedge;
        Halfedge_index he2 = mesh.next(he1);
        Halfedge_index he3 = mesh.next(he2);

        // Store original vertices
        Vertex_index v1 = mesh.source(he1);
        Vertex_index v2 = mesh.target(he1);
        Vertex_index v3 = mesh.target(he2);

        // Create new halfedges connecting to new_vertex
        Halfedge_index new_he1 = mesh.add_edge(v2, new_vertex);
        Halfedge_index new_he2 = mesh.add_edge(v3, new_vertex);
        Halfedge_index new_he3 = mesh.add_edge(v1, new_vertex);

        // Get opposite halfedges
        Halfedge_index new_he1_opp = mesh.opposite(new_he1);
        Halfedge_index new_he2_opp = mesh.opposite(new_he2);
        Halfedge_index new_he3_opp = mesh.opposite(new_he3);

        // Create new faces
        Face_index new_triangle1 = mesh.add_face();
        Face_index new_triangle2 = mesh.add_face();
        // old_face will be reused for the third triangle

        // Set faces for halfedges - First triangle
        mesh.set_face(he1, old_face);
        mesh.set_face(new_he1, old_face);
        mesh.set_face(new_he3_opp, old_face);

        // Second triangle
        mesh.set_face(he2, new_triangle1);
        mesh.set_face(new_he2, new_triangle1);
        mesh.set_face(new_he1_opp, new_triangle1);

        // Third triangle
        mesh.set_face(he3, new_triangle2);
        mesh.set_face(new_he3, new_triangle2);
        mesh.set_face(new_he2_opp, new_triangle2);

        // Set next relationships - First triangle
        mesh.set_next(he1, new_he1);
        mesh.set_next(new_he1, new_he3_opp);
        mesh.set_next(new_he3_opp, he1);

        // Second triangle
        mesh.set_next(he2, new_he2);
        mesh.set_next(new_he2, new_he1_opp);
        mesh.set_next(new_he1_opp, he2);

        // Third triangle
        mesh.set_next(he3, new_he3);
        mesh.set_next(new_he3, new_he2_opp);
        mesh.set_next(new_he2_opp, he3);

        // Set halfedge for faces
        mesh.set_halfedge(old_face, he1);
        mesh.set_halfedge(new_triangle1, he2);
        mesh.set_halfedge(new_triangle2, he3);

        // Set halfedge for the new vertex
        mesh.set_halfedge(new_vertex, new_he1);

        // // Debug output
        // std::cout << "Created new halfedges: " << new_he1 << ", " << new_he2 << ", " << new_he3 << std::endl;
        // std::cout << "Created new faces: " << new_triangle1 << ", " << new_triangle2 << std::endl;

        // Debug verification
        if (mesh.is_valid(false)) {
            // std::cout << "Triangulation completed successfully" << std::endl;
        } else {
            std::cout << "Mesh validation failed after triangulation" << std::endl;
        }
        std::vector<Halfedge_index> new_halfedges = {new_he1,new_he2,new_he3};

        return  new_halfedges;


    } catch (const std::exception& e) {
        std::cerr << "Error during triangulation: " << e.what() << std::endl;
        throw;
    }
    
}
std::vector<Halfedge_index> insert_point_on_edge(Halfedge_index halfedge, Vertex_index new_vertex, CGL::Mesh& mesh) {
    try {
        Halfedge_index he1 = CGAL::Euler::join_face(halfedge,mesh);
        Face_index old_face = mesh.face(he1);
        Halfedge_index he2 = mesh.next(he1);
        Halfedge_index he3 = mesh.next(he2);
        Halfedge_index he4 = mesh.next(he3);

        // Store original vertices starting from he1
        Vertex_index v1 = mesh.source(he1);
        Vertex_index v2 = mesh.target(he1);
        Vertex_index v3 = mesh.target(he2);
        Vertex_index v4 = mesh.target(he3);

        // Create new halfedges in CCW order around new_vertex
        Halfedge_index new_he1 = mesh.add_edge(v2, new_vertex);
        Halfedge_index new_he2 = mesh.add_edge(v3, new_vertex);
        Halfedge_index new_he3 = mesh.add_edge(v4, new_vertex);
        Halfedge_index new_he4 = mesh.add_edge(v1, new_vertex);

        // Get opposite halfedges
        Halfedge_index new_he1_opp = mesh.opposite(new_he1);
        Halfedge_index new_he2_opp = mesh.opposite(new_he2);
        Halfedge_index new_he3_opp = mesh.opposite(new_he3);
        Halfedge_index new_he4_opp = mesh.opposite(new_he4);

        // Create new faces
        Face_index new_triangle1 = mesh.add_face();
        Face_index new_triangle2 = mesh.add_face();
        Face_index new_triangle3 = mesh.add_face();

        // Set faces for halfedges in CCW order
        // First triangle (reusing old_face)
        mesh.set_face(he1, old_face);
        mesh.set_face(new_he1, old_face);
        mesh.set_face(new_he4_opp, old_face);

        // Second triangle
        mesh.set_face(he2, new_triangle1);
        mesh.set_face(new_he2, new_triangle1);
        mesh.set_face(new_he1_opp, new_triangle1);

        // Third triangle
        mesh.set_face(he3, new_triangle2);
        mesh.set_face(new_he3, new_triangle2);
        mesh.set_face(new_he2_opp, new_triangle2);

        // Fourth triangle
        mesh.set_face(he4, new_triangle3);
        mesh.set_face(new_he4, new_triangle3);
        mesh.set_face(new_he3_opp, new_triangle3);

        // Set next relationships maintaining CCW order
        // First triangle
        mesh.set_next(he1, new_he1);
        mesh.set_next(new_he1, new_he4_opp);
        mesh.set_next(new_he4_opp, he1);

        // Second triangle
        mesh.set_next(he2, new_he2);
        mesh.set_next(new_he2, new_he1_opp);
        mesh.set_next(new_he1_opp, he2);

        // Third triangle
        mesh.set_next(he3, new_he3);
        mesh.set_next(new_he3, new_he2_opp);
        mesh.set_next(new_he2_opp, he3);

        // Fourth triangle
        mesh.set_next(he4, new_he4);
        mesh.set_next(new_he4, new_he3_opp);
        mesh.set_next(new_he3_opp, he4);

        // Set halfedge for faces
        mesh.set_halfedge(old_face, he1);
        mesh.set_halfedge(new_triangle1, he2);
        mesh.set_halfedge(new_triangle2, he3);
        mesh.set_halfedge(new_triangle3, he4);

        // Set halfedge for the new vertex
        mesh.set_halfedge(new_vertex, new_he1);

        // Debug verification
        if (mesh.is_valid(false)) {
            // std::cout << "Triangulation completed successfully on edge" << std::endl;
        } else {
            std::cout << "Mesh validation failed after triangulation ON EDGE" << std::endl;
        }

        return {new_he1, new_he2, new_he3, new_he4};
    } catch (const std::exception& e) {
        std::cerr << "Error during point insertion on edge: " << e.what() << std::endl;
        throw;
    }
}

void legalize_edge(Halfedge_index hedge, CGL::Mesh& mesh) {
    try {

        // Check if the edge is on the boundary
        if (mesh.is_border(hedge) || mesh.is_border(mesh.opposite(hedge))) {
            return;
        }

        // Get the two triangles sharing the edge
        Face_index face1 = mesh.face(hedge);
        Face_index face2 = mesh.face(mesh.opposite(hedge));

        // Get vertices of both triangles
        Vertex_index a = mesh.source(hedge);
        Vertex_index b = mesh.target(hedge);
        Vertex_index c = mesh.target(mesh.next(hedge));  // vertex in triangle 1
        Vertex_index d = mesh.target(mesh.next(mesh.opposite(hedge)));  // vertex in triangle 2


        // Check if edge needs to be flipped (Delaunay criterion)
        if (inCircle(mesh.point(a), mesh.point(b), mesh.point(c), mesh.point(d)) ||
            inCircle(mesh.point(a), mesh.point(d), mesh.point(b), mesh.point(c)) ) {
            // std::cout << "Flipping edge between vertices " << a << " and " << b << std::endl;

            // CGAL::draw(mesh);
            // Flip the edge
            CGAL::Euler::flip_edge(hedge, mesh);
            // CGAL::draw(mesh);

            Halfedge_index h1,h2,h3,h4;
            h1 = mesh.next(hedge);
            h2 = mesh.prev(hedge);
            h3 = mesh.next(mesh.opposite(hedge));
            h4 = mesh.prev(mesh.opposite(hedge));

            // Recursively legalize the affected edges
            legalize_edge(h1, mesh);
            legalize_edge(h2, mesh);
            legalize_edge(h3, mesh);
            legalize_edge(h4, mesh);

            // Debug verification after flip
            if (mesh.is_valid(false)) {
                // std::cout << "Edge flip completed successfully" << std::endl;
            } else {
                std::cout << "Mesh validation failed after flip" << std::endl;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error during edge legalization: " << e.what() << std::endl;
        throw;
    }
}


//This function receives as parameter three points
float getDeterminant(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3){
    // Formula to get determinant with the points
    return ((p2.x() - p1.x())*(p3.y() - p1.y()) - (p2.y() - p1.y())*(p3.x() - p1.x()));
}

//This function receives as parameter a triangle in a Point2 Vector
float getDeterminant(std::vector<CGL::Point2> pts){
    int x[3] = {0}, y[3] = {0};

    for(int i = 0;i < 3;i++){
        x[i] = pts[i].x();
        y[i] = pts[i].y();
    }
    // Formula to get determinant with the points
    return ((x[1] - x[0])*(y[2] - y[0]) - (y[1] - y[0])*(x[2] - x[0]));
}

bool inCircle(CGL::Point3 a, CGL::Point3 b,CGL::Point3 c,CGL::Point3 d){
    double matrix[4][4] = {{a.x(),a.y(), (a.x()*a.x()) + (a.y()*a.y()),1},
                           {b.x(),b.y(), (b.x()*b.x()) + (b.y()*b.y()),1},
                           {c.x(),c.y(), (c.x()*c.x()) + (c.y()*c.y()),1},
                           {d.x(),d.y(), (d.x()*d.x()) + (d.y()*d.y()),1}};

    //Not includind if the point is in the circumference
    if(getDeterminant4x4(matrix, 4) > 0)
        return true;
    return false;
}

double getDeterminant4x4(double matrix[4][4], int n = 4){
    double det = 0;
    double submatrix[4][4];
     
    if(n == 1){
        return matrix[0][0];
    }
    else if(n == 2){
        return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
    }
    else{
        for(int x = 0; x < n; x++){
            int subi = 0;
            for(int i = 1; i < n; i++){
                int subj = 0;
                for(int j = 0; j < n; j++){
                    if(j == x){
                        continue;
                    }
                    submatrix[subi][subj] = matrix[i][j];
                    subj++;
                }
                subi++;
            }
            det = det + (pow(-1, x) * matrix[0][x] * getDeterminant4x4(submatrix, n-1));
        }
    }
    return det;
}

//This function receives as parameter three points
bool getOrientationTriangle(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3){
    return (getDeterminant(p1,p2,p3) >= 0);
}

bool getOrientationTriangle(CGL::Point3 p1, CGL::Point3 p2,CGL::Point3 p3){
    bool left1 = CGL::Left(p1, p2, p3);
    bool left2 = CGL::Left(p2, p3, p1);
    bool left3 = CGL::Left(p3, p1, p2);
    return (left1 && left2 && left3);
}



// Returns true if point is inside triangle
bool isInsideTriangle(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3, CGL::Point2 point){
    bool orientation1 = getOrientationTriangle(p1,p2,point);
    bool orientation2 = getOrientationTriangle(p2,p3,point);
    bool orientation3 = getOrientationTriangle(p3,p1,point);

    return (orientation1 == orientation2 && orientation2 == orientation3);
}


std::vector<CGL::Mesh::Vertex_index> get_face_vertices(CGL::Mesh::Face_index face, CGL::Mesh polygon) {
    CGL::Mesh::Halfedge_index first_hedge = polygon.halfedge(face);
    CGL::Mesh::Halfedge_index actual_hedge = first_hedge;

    std::vector<CGL::Mesh::Vertex_index> triangle_vertices; //Vector that stores all vertices of given triangle

    if(face == polygon.null_face())
        printf("\n\nAQUI\n\n");

        do {
            //saves on triangle_vertices the vertices
            triangle_vertices.push_back(polygon.source(actual_hedge));

            CGL::Mesh::Halfedge_index next_hedge = polygon.next(actual_hedge);

            //Gets the next edge of the triangle
            while(polygon.face(actual_hedge) != polygon.face(next_hedge)) {
                next_hedge = polygon.next(polygon.opposite(actual_hedge));
            }
            actual_hedge = next_hedge; // goes to next hedge

        }while(actual_hedge != first_hedge);

        for (size_t i = 0; i < triangle_vertices.size(); ++i) {
            CGL::Point3 pt = polygon.point(triangle_vertices[i]); 
        }
    return triangle_vertices;
}


std::vector<CGL::Mesh::Halfedge_index> get_face_halfedges(CGL::Mesh::Face_index face, CGL::Mesh polygon) {
    CGL::Mesh::Halfedge_index first_hedge = polygon.halfedge(face);
    CGL::Mesh::Halfedge_index actual_hedge = first_hedge;

    std::vector<CGL::Mesh::Halfedge_index> face_halfedges; // Vector that stores all halfedges of the given face

    do {
        // Save the current halfedge
        face_halfedges.push_back(actual_hedge);
    
        CGL::Mesh::Halfedge_index next_hedge = polygon.next(actual_hedge);

        // Get the next halfedge in the face
        while (polygon.face(actual_hedge) != polygon.face(next_hedge)) {
            next_hedge = polygon.next(polygon.opposite(actual_hedge));
        }
        actual_hedge = next_hedge; // Move to the next halfedge

    } while (actual_hedge != first_hedge);

    return face_halfedges;
}

bool isInsideTriangle(Point3 v1, Point3 v2, Point3 v3, Point3 pt) {
    // Check if pt is to the left of all three edges
    // For CCW triangle, point must be to the left of all edges
    bool left1 = CGL::Left(v1, v2, pt);
    bool left2 = CGL::Left(v2, v3, pt);
    bool left3 = CGL::Left(v3, v1, pt);

    // Point is inside if it's to the left of all edges
    return (left1 && left2 && left3);
}

bool isOnTriangle(Point3 v1, Point3 v2, Point3 v3, Point3 pt) {
    // Check if pt is to the left of all three edges
    // For CCW triangle, point must be to the left of all edges
    bool left1 = CGL::LeftOn(v1, v2, pt);
    bool left2 = CGL::LeftOn(v2, v3, pt);
    bool left3 = CGL::LeftOn(v3, v1, pt);

    // Point is inside if it's to the left of all edges
    return (left1 && left2 && left3);
}

void remove_super_triangle(CGL::Mesh& mesh, const std::vector<Vertex_index>& super_indexes) {
    // First, collect all faces and edges to remove
    std::vector<Face_index> faces_to_remove;

    // Identify faces connected to super triangle vertices
    for (Face_index f : mesh.faces()) {
        bool connected_to_super = false;
        for (Vertex_index v : mesh.vertices_around_face(mesh.halfedge(f))) {
            if (std::find(super_indexes.begin(), super_indexes.end(), v) != super_indexes.end()) {
                connected_to_super = true;
                break;
            }
        }
        if (connected_to_super) {
            faces_to_remove.push_back(f);
            Halfedge_index aux_hedge = mesh.halfedge(f);
            mesh.set_face(aux_hedge, mesh.null_face());
            mesh.set_face(mesh.next(aux_hedge), mesh.null_face());
            mesh.set_face(mesh.prev(aux_hedge), mesh.null_face());
        }
    }

    // Remove faces
    for (Face_index f : faces_to_remove) {
        mesh.remove_face(f);
    }

    // Join faces where both sides of an edge are null
    for (Halfedge_index h : mesh.halfedges()) {
        if (mesh.face(h) == mesh.null_face() && 
            mesh.face(mesh.opposite(h)) == mesh.null_face()) {
            CGAL::Euler::join_face(h, mesh);
        }
    }

    // Remove super triangle vertices
    for (Vertex_index v : super_indexes) {
        mesh.remove_vertex(v);
    }

    // Clean up the mesh
    mesh.collect_garbage();
}

cv::Mat convert_to_pgm(const std::string& image_path, const std::string& output_path) {
    // Check if file exists
    struct stat buffer;
    if (stat(image_path.c_str(), &buffer) != 0) {
        std::cout << "File does not exist: " << image_path << std::endl;
        return cv::Mat();
    }

    // Read image in color
    cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
    if (img.empty()) {
        std::cout << "Could not read the image: " << image_path << std::endl;
        return cv::Mat();
    }

    std::cout << "Original loaded image size: " << img.size() << std::endl;

     // Convert to grayscale
     cv::Mat gray_img;
     cv::cvtColor(img, gray_img, cv::COLOR_BGR2GRAY);
 
    // Resize to 512x512 if it's not already that size
    if (gray_img.size() != cv::Size(512, 512)) {
        cv::Mat resized;
        cv::resize(gray_img, resized, cv::Size(512, 512));
        gray_img = resized;
        std::cout << "Resized image to 512x512" << std::endl;
    }

    // Save the grayscale image in PGM format
    if (cv::imwrite(output_path, gray_img)) {
        std::cout << "PGM file saved: " << output_path << std::endl;
    } else {
        std::cout << "Failed to save PGM file: " << output_path << std::endl;
    }

    return gray_img;
}

void add_border(CGL::Mesh& mesh) {
    // Find a border halfedge (one with a null face)
    Halfedge_index first_hedge;
    bool found_border = false;
    for (Halfedge_index h : mesh.halfedges()) {
        if (mesh.face(h) == mesh.null_face()) {
            first_hedge = h;
            found_border = true;
            break;
        }
    }
    if (!found_border) {
        std::cout << "No border halfedge found" << std::endl;
        return;
    }

    // Use a set to record visited border vertices (using target vertices)
    std::set<Vertex_index> visited;
    Vertex_index start_vertex = mesh.target(first_hedge);
    visited.insert(start_vertex);

    Halfedge_index actual_hedge = first_hedge;

    try {
        while (true) {
            // Get the next halfedges in the current border cycle
            Halfedge_index next_hedge = mesh.next(actual_hedge);
            Halfedge_index next_next_hedge = mesh.next(next_hedge);
            Halfedge_index next3_hedge = mesh.next(next_next_hedge);

            // Get the three consecutive border points
            Point3 p1 = mesh.point(mesh.target(actual_hedge));
            Point3 p2 = mesh.point(mesh.target(next_hedge));
            Point3 p3 = mesh.point(mesh.target(next_next_hedge));

            // Check if these points form a valid triangle
            if (getOrientationTriangle(p1, p2, p3)) {
                try {
                    // Create new face and edge
                    Halfedge_index new_he = CGAL::Euler::add_face_to_border(actual_hedge,next_next_hedge,mesh);
                    Halfedge_index new_he_opp = mesh.opposite(new_he);

                    // Continue with the new border edge
                    actual_hedge = new_he_opp;
                } catch (const std::exception& e) {
                    std::cerr << "Error creating face: " << e.what() << std::endl;
                    break;
                }
            } else {
                // Move along the border
                actual_hedge = next_hedge;
            }

            // Check for cycle completion
            Vertex_index current_vertex = mesh.target(actual_hedge);
            if (visited.find(current_vertex) != visited.end()) {
                break;
            }
            visited.insert(current_vertex);

            // Validate mesh after each iteration
            if (!mesh.is_valid(false)) {
                std::cerr << "Invalid mesh state" << std::endl;
                break;
            }
        }


    } catch (const std::exception& e) {
        std::cerr << "Error in add_border: " << e.what() << std::endl;
    }

    // Final mesh cleanup
    mesh.collect_garbage();
}


// Compute barycentric coordinates for point p in triangle (v0,v1,v2)
void computeBarycentrics(const Point3& p, const Point3& v0, const Point3& v1, const Point3& v2,
                         double &lambda0, double &lambda1, double &lambda2)
{
    double denom = (v1.x() - v0.x()) * (v2.y() - v0.y()) - (v2.x() - v0.x()) * (v1.y() - v0.y());

    if (std::abs(denom) < 1e-10) { // Prevent division by zero
        lambda0 = lambda1 = lambda2 = -1.0;
        return;
    }
    lambda0 = ((v1.x() - p.x()) * (v2.y() - p.y()) - (v2.x() - p.x()) * (v1.y() - p.y())) / denom;
    lambda1 = ((v2.x() - p.x()) * (v0.y() - p.y()) - (v0.x() - p.x()) * (v2.y() - p.y())) / denom;
    lambda2 = 1.0 - lambda0 - lambda1;
}


uchar bilinearInterpolate(const cv::Mat& img, float x, float y) {
    // Get the four surrounding pixel coordinates
    int x0 = std::floor(x);
    int x1 = x0 + 1;
    int y0 = std::floor(y);
    int y1 = y0 + 1;

    // Calculate interpolation weights
    float wx = x - x0;
    float wy = y - y0;

    // image boundaries
    x0 = std::min(std::max(x0, 0), img.cols - 1);
    x1 = std::min(std::max(x1, 0), img.cols - 1);
    y0 = std::min(std::max(y0, 0), img.rows - 1);
    y1 = std::min(std::max(y1, 0), img.rows - 1);

    // Get pixel values
    uchar p0 = img.at<uchar>(y0, x0);
    uchar p1 = img.at<uchar>(y0, x1);
    uchar p2 = img.at<uchar>(y1, x0);
    uchar p3 = img.at<uchar>(y1, x1);

    // Perform bilinear interpolation
    float value = (1 - wx) * (1 - wy) * p0 +
                 wx * (1 - wy) * p1 +
                 (1 - wx) * wy * p2 +
                 wx * wy * p3;

    return static_cast<uchar>(std::round(value));
}
cv::Mat rasterizeMesh(const Mesh& mesh, const cv::Mat& inputImg)
{
    // Create output image
    cv::Mat output = cv::Mat::zeros(inputImg.size(), inputImg.type());
    
    // Create a debug image for visualization
    cv::Mat debug_image = inputImg.clone();
    cv::cvtColor(debug_image, debug_image, cv::COLOR_GRAY2BGR);

    for (Face_index f : mesh.faces())
    {
        std::vector<Vertex_index> verts = get_face_vertices(f, mesh);
        if (verts.size() != 3) continue;

        Point3 v0 = mesh.point(verts[0]);
        Point3 v1 = mesh.point(verts[1]);
        Point3 v2 = mesh.point(verts[2]);

        double min_x = std::min({ v0.x(), v1.x(), v2.x() });
        double max_x = std::max({ v0.x(), v1.x(), v2.x() });
        double min_y = std::min({ v0.y(), v1.y(), v2.y() });
        double max_y = std::max({ v0.y(), v1.y(), v2.y() });

        int box_min_x = std::max(0, (int)std::floor(min_x));
        int box_max_x = std::min(inputImg.cols - 1, (int)std::ceil(max_x));
        int box_min_y = std::max(0, (int)std::floor(min_y));
        int box_max_y = std::min(inputImg.rows - 1, (int)std::ceil(max_y));

        uchar c0 = bilinearInterpolate(inputImg, v0.x(), v0.y());
        uchar c1 = bilinearInterpolate(inputImg, v1.x(), v1.y());
        uchar c2 = bilinearInterpolate(inputImg, v2.x(), v2.y());

        std::vector<cv::Point> triangle_pts = {
            cv::Point((int)v0.x(), (int)v0.y()),
            cv::Point((int)v1.x(), (int)v1.y()),
            cv::Point((int)v2.x(), (int)v2.y())
        };

        cv::Scalar color = cv::Scalar(255, 0, 0);
        for (Halfedge_index h : mesh.halfedges_around_face(mesh.halfedge(f))) {
            if (mesh.is_border(mesh.opposite(h))) {
                color = cv::Scalar(0, 0, 255);
                break;
            }
        }

        cv::line(debug_image, triangle_pts[0], triangle_pts[1], color, 1);
        cv::line(debug_image, triangle_pts[1], triangle_pts[2], color, 1);
        cv::line(debug_image, triangle_pts[2], triangle_pts[0], color, 1);

        for (int y = box_min_y; y <= box_max_y; ++y) {
            for (int x = box_min_x; x <= box_max_x; ++x) {
                Point3 p((double)x + 0.5, (double)y + 0.5, 0.0);
                double lambda0, lambda1, lambda2;
                computeBarycentrics(p, v0, v1, v2, lambda0, lambda1, lambda2);

                const double eps = 1e-6;
                if (lambda0 >= -eps && lambda1 >= -eps && lambda2 >= -eps) {
                    double pixel_val = lambda0 * c0 + lambda1 * c1 + lambda2 * c2;
                    output.at<uchar>(y, x) = (uchar)std::round(pixel_val);
                }
            }
        }
    }

    cv::imshow("debug_triangulation", debug_image);
    return output;
}


std::vector<Point2> select_mesh_points(const cv::Mat& img_pgm, unsigned long long amount, bool add_random_points) {
    try {
        std::vector<Point2> point_vec;
        
        // Function to check if a point already exists in the vector.
        auto point_exists = [&point_vec](const Point2& point) -> bool {
            return std::find(point_vec.begin(), point_vec.end(), point) != point_vec.end();
        };

        // Use goodFeaturesToTrack to detect corners.
        std::vector<cv::Point2f> corners;
        double qualityLevel = 0.005;     // minimal accepted quality
        double minDistance  = 1;        // minimum distance between corners
        cv::goodFeaturesToTrack(img_pgm, corners, amount, qualityLevel, minDistance);
        for (const auto &pt : corners) {
            Point2 p(pt.x, pt.y);
            if (!point_exists(p)) {
                point_vec.push_back(p);
            }
        }

        std::cout << "\nAmount of corner points: " << point_vec.size() << std::endl;

        //Optionally add random points if we still haven't reached the desired count.
        if (point_vec.size() < amount && add_random_points) {
            cv::RNG rng(cv::getTickCount());
            while (point_vec.size() < amount) {
                int x = rng.uniform(0, img_pgm.cols);
                int y = rng.uniform(0, img_pgm.rows);
                Point2 p(x, y);
                if (!point_exists(p)) {
                    point_vec.push_back(p);
                }
            }
        }

        std::cout << "\nAmount of points: " << point_vec.size() << std::endl;

        //Add border points to ensure the triangulation covers the image boundary.
        addBorderPoints(img_pgm, point_vec, 16);

        // Debug visualization.
        cv::Mat debug_img;
        cv::cvtColor(img_pgm, debug_img, cv::COLOR_GRAY2BGR);
        for (const auto &pt : point_vec) {
            cv::circle(debug_img, cv::Point(pt.x(), pt.y()), 1, cv::Scalar(0, 0, 255), -1);
        }
        // cv::imshow("Selected Points", debug_img);
        // cv::waitKey(1);

        return point_vec;
    }
    catch (const std::exception &e) {
        std::cerr << "Error in select_mesh_points: " << e.what() << std::endl;
        return std::vector<Point2>();
    }
}


void addBorderPoints(const cv::Mat& img_pgm, std::vector<Point2>& point_vec, int border_spacing) {
    // Function to check if a point exists
    auto point_exists = [&point_vec](const Point2& point) -> bool {
        return std::find(point_vec.begin(), point_vec.end(), point) != point_vec.end();
    };

    // Add top border points
    for (int x = 0; x < img_pgm.cols; x += border_spacing) {
        Point2 new_point(x, 0);
        if (!point_exists(new_point)) {
            point_vec.push_back(new_point);
        }
    }

    // Add bottom border points
    for (int x = 0; x < img_pgm.cols; x += border_spacing) {
        Point2 new_point(x, img_pgm.rows - 1);
        if (!point_exists(new_point)) {
            point_vec.push_back(new_point);
        }
    }

    // Add left border points (excluding corners already added)
    for (int y = border_spacing; y < img_pgm.rows - border_spacing; y += border_spacing) {
        Point2 new_point(0, y);
        if (!point_exists(new_point)) {
            point_vec.push_back(new_point);
        }
    }

    // Add right border points (excluding corners already added)
    for (int y = border_spacing; y < img_pgm.rows - border_spacing; y += border_spacing) {
        Point2 new_point(img_pgm.cols - 1, y);
        if (!point_exists(new_point)) {
            point_vec.push_back(new_point);
        }
    }

    // Add the four corner points
    Point2 corners[] = {
        Point2(0, 0),
        Point2(img_pgm.cols - 1, 0),
        Point2(0, img_pgm.rows - 1),
        Point2(img_pgm.cols - 1, img_pgm.rows - 1)
    };

    for (const auto& corner : corners) {
        if (!point_exists(corner)) {
            point_vec.push_back(corner);
        }
    }
}
