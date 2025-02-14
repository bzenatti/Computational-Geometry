#include <Types.h>
#include <CGAL/draw_polygon_2.h>
#include <Random.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <CGAL/draw_surface_mesh.h>


#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random.h>

#include <Partition.h>
#include <IO.h>

#define PI 3.14159265358979323846
#define MAXN 100

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
float getTriangleArea(std::vector<CGL::Point2> pts);

double dotProduct(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3,CGL::Point2 p4); // p1p2 and p3p4
double vectorMagnitude(CGL::Point2 p1, CGL::Point2 p2);
double angleBetween(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3,CGL::Point2 p4);

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

bool isEdgeShared(CGL::Mesh::Halfedge_index hedge, std::vector<CGL::Mesh::Face_index> triangle, CGL::Mesh mesh);

int main(int argc, char* argv[]) {

    /* Vector to store points. */
    std::vector<Point3> pts;
    std::vector<Vertex_index> indexes;

    /* Simple mesh with random vertices. */
    Mesh mesh_delaunay;

    int nv = 1000;    /* Number of vertices. */

    std::vector<Point2>   point_vec;
    CGAL::Random         rand;

    rand.get_seed();

    
    //generate for representing 512x512 pixels
    CGAL::copy_n_unique(Point_generator(512), nv, std::back_inserter(point_vec));

    Mesh mesh;
    mesh = delaunay_triangulation(point_vec);

    // CGL::print_mesh(mesh);
    CGAL::draw(mesh);
   
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
    // Making it 20% larger than necessary to ensure all points are inside
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

    for (auto p = vertices.begin(); p != vertices.end(); ++p) {
        insert_point(delaunay_mesh,*p);
    }

    std::cout << "Ended inserting points" << std::endl;

    // CGL::print_mesh(delaunay_mesh);

    remove_super_triangle(delaunay_mesh,super_indexes);

    for(unsigned i = 0;i < 3;i++)
        for(Halfedge_index h : delaunay_mesh.halfedges())
            legalize_edge(h,delaunay_mesh);

    if (!delaunay_mesh.is_valid(false))
        printf("Deu ruim");

    // CGL::print_mesh_vertices(delaunay_mesh);

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
            std::cout << "Triangulation completed successfully" << std::endl;
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
            std::cout << "Triangulation completed successfully on edge" << std::endl;
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
        if (inCircle(mesh.point(a), mesh.point(b), mesh.point(c), mesh.point(d))) {
            // std::cout << "Flipping edge between vertices " << a << " and " << b << std::endl;

            // CGAL::draw(mesh);
            // Flip the edge
            CGAL::Euler::flip_edge(hedge, mesh);
            // CGAL::draw(mesh);

            // Recursively legalize the affected edges
            legalize_edge(mesh.next(hedge), mesh);
            legalize_edge(mesh.prev(hedge), mesh);

            // Debug verification after flip
            if (mesh.is_valid(false)) {
                std::cout << "Edge flip completed successfully" << std::endl;
            } else {
                std::cout << "Mesh validation failed after flip" << std::endl;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error during edge legalization: " << e.what() << std::endl;
        throw;
    }
}

bool isEdgeShared(Halfedge_index hedge, std::vector<Face_index> triangle, CGL::Mesh mesh){
    for (Face_index f : triangle){
        if(mesh.face(mesh.opposite(hedge)).idx() == f)
            return true;
    }

    return false;
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

// Considering two vectors, p1p2==u and p3p4==v
double dotProduct(CGL::Point3 p1, CGL::Point3 p2,CGL::Point3 p3,CGL::Point3 p4) {
    int ux = p2.x() - p1.x(), uy = p2.y() - p1.y();
    int vx = p4.x() - p3.x(), vy = p4.y() - p3.y();

    return (ux*vx) + (uy*vy);
}

double vectorMagnitude(CGL::Point3 p1, CGL::Point3 p2){
    return std::sqrt(dotProduct(p1,p2,p1,p2));
}

//returns in radians the angle between two vectors, p1p2==u and p3p4==v
double angleBetween(CGL::Point3 p1, CGL::Point3 p2,CGL::Point3 p3,CGL::Point3 p4){
    return std::acos((dotProduct(p1,p2,p3,p4))/(vectorMagnitude(p1,p2) * vectorMagnitude(p3,p4)));
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