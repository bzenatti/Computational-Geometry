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

bool inCircle(CGL::Point3 a, CGL::Point3 b,CGL::Point3 c,CGL::Point3 d);

bool getOrientation(CGL::Polygon2 P);
bool getOrientationTriangle(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3);
float getTriangleArea(std::vector<CGL::Point2> pts);

double dotProduct(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3,CGL::Point2 p4); // p1p2 and p3p4
double vectorMagnitude(CGL::Point2 p1, CGL::Point2 p2);
double angleBetween(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3,CGL::Point2 p4);

bool isVerticeConvex(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3,bool polygonOrientation);
bool isVerticeEar(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3, CGL::Polygon2 P);

bool isInsideTriangle(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3, CGL::Point2 point);

double triangle_goodness (CGL::Point3 p0, CGL::Point3 p1, CGL::Point3 p2);
void edge_flip(CGL::Mesh::Halfedge_index& e, CGL::Mesh& polygon);

std::vector<CGL::Mesh::Vertex_index> get_face_vertices(CGL::Mesh::Face_index face, CGL::Mesh polygon);
std::vector<CGL::Mesh::Halfedge_index> get_face_halfedges(CGL::Mesh::Face_index face, CGL::Mesh polygon);

bool compare_goodness(CGL::Mesh::Halfedge_index first_triangle_hedge, CGL::Mesh::Halfedge_index adjacent_triangle_hedge, CGL::Mesh polygon);
bool can_flip(CGL::Mesh::Halfedge_index initial_hedge, CGL::Mesh polygon);

std::vector<CGL::Point2>  createSuperTriangle(const std::vector<CGL::Point2>& points);
CGL::Mesh delaunay_triangulation(std::vector<CGL::Point2> mesh_vertices);
void insert_point(CGL::Mesh& mesh, CGL::Point2 vi);
void insert_triangles(Halfedge_index halfedge,Vertex_index new_vertex, CGL::Mesh& mesh);
bool isEdgeShared(CGL::Mesh::Halfedge_index hedge, std::vector<CGL::Mesh::Face_index> triangle, CGL::Mesh mesh);

int main(int argc, char* argv[]) {

    /* Vector to store points. */
    std::vector<Point3> pts;
    std::vector<Vertex_index> indexes;

    /* Simple mesh with random vertices. */
    Mesh mesh_delaunay;

    int nv = 4;    /* Number of vertices. */

    std::vector<Point2>   point_vec;
    CGAL::Random         rand;

    rand.get_seed();

    
    //generate for representing 512x512 pixels
    CGAL::copy_n_unique(Point_generator(20), nv, std::back_inserter(point_vec));

    /*     
    for (const auto& point : point_set) {
        std::cout << point << std::endl;
    }
    */

    Mesh mesh;
    
    int i = 0 ;

    

    mesh = delaunay_triangulation(point_vec);

    CGL::print_mesh(mesh);

    //printf("teste2");
    //CGAL::draw(mesh);
   
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
        // std::cout << delaunay_mesh.opposite(delaunay_mesh.halfedge(vertex)) << 
        // delaunay_mesh.face(delaunay_mesh.opposite(delaunay_mesh.halfedge(vertex))) << std::endl;
    }

    insert_point(delaunay_mesh,vertices[0]);

    // for (auto p = vertices.begin(); p != vertices.end(); ++p) {
    //     insert_point(delaunay_mesh,*p);
    //     printf("Test");
    // }
    

    // //Remove super triangle
    // for(CGL::Mesh::Halfedge_index hedge : super_hedges){
    //     CGAL::Euler::remove_face(hedge,delaunay_mesh);
    // } 

    return delaunay_mesh;
}

void insert_point(CGL::Mesh& mesh, CGL::Point2 vi){

    Point3 p3(vi.x(), vi.y(), 0.0);  // Convert 2D point to 3D 
    Vertex_index new_vertice = mesh.add_vertex(p3);

    //Get face the point is inside
    Face_index boundary_triangle;
    for(Face_index f : mesh.faces()){
        Halfedge_index hedge = mesh.halfedge(f);
        Point2 t1(mesh.point(mesh.target(hedge)).x(), mesh.point(mesh.target(hedge)).y());
        Point2 t2(mesh.point(mesh.source(hedge)).x(), mesh.point(mesh.source(hedge)).y());
        Point2 t3(mesh.point(mesh.target(mesh.next(hedge))).x(), mesh.point(mesh.target(mesh.next(hedge))).y());

        if(isInsideTriangle(t1,t2,t3,vi)){
            boundary_triangle = f;
            break;
        }
    }

    Halfedge_index halfedge = mesh.halfedge(boundary_triangle);

    std::vector<Vertex_index> face_vertexes = get_face_vertices(boundary_triangle,mesh);
    //add_face do cgal n funciona pq ele recebe vertice como parametro, 
    //e os vertices aponta para as twins

    insert_triangles(halfedge,new_vertice,mesh);

    //Legalize Edges


    CGL::print_mesh(mesh);
    CGAL::draw(mesh);

}
void insert_triangles(Halfedge_index halfedge, Vertex_index new_vertex, CGL::Mesh& mesh) {
    // Debug output
    std::cout << "Starting insert_triangles operation" << std::endl;
    std::cout << "Initial halfedge: " << halfedge << ", New vertex: " << new_vertex << std::endl;

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

        // Debug output
        std::cout << "Created new halfedges: " << new_he1 << ", " << new_he2 << ", " << new_he3 << std::endl;
        std::cout << "Created new faces: " << new_triangle1 << ", " << new_triangle2 << std::endl;

        // Debug verification
        if (mesh.is_valid(false)) {
            std::cout << "Triangulation completed successfully" << std::endl;
        } else {
            std::cout << "Warning: Mesh validation failed after triangulation" << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error during triangulation: " << e.what() << std::endl;
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
// True for counterclockwise and false for clockwise
// Default: straight line is set for counterclockwise
bool getOrientation(CGL::Polygon2 P){

    int totalSum = 0;
    size_t P_size = P.size();
    for (size_t i = 0; i < P_size; i++) {
        // % operation for when get to i == P_size - 1, when adds 1 or 2, go to 0 and 1
        int p0 = i, p1 = ((i+1) % P_size);
        //If totalSum is positive, is cw
        totalSum += ((P[p1].x() - P[p0].x()) * (P[p1].y() + P[p0].y()) );
    }

    return !(totalSum >= 0);
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

bool isVerticeConvex(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3,bool polygonOrientation){
    if(getOrientationTriangle(p1,p2,p3) == polygonOrientation)
        return true;
    else 
        return false;
}

// Returns true if point is inside triangle
bool isInsideTriangle(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3, CGL::Point2 point){
    bool orientation1 = getOrientationTriangle(p1,p2,point);
    bool orientation2 = getOrientationTriangle(p2,p3,point);
    bool orientation3 = getOrientationTriangle(p3,p1,point);

    return (orientation1 == orientation2 && orientation2 == orientation3);
}

bool isVerticeEar(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3, CGL::Polygon2 P){
    size_t P_size = P.size();
    size_t i = P_size;

    // If vertice is concave, returns false
    if(!isVerticeConvex(p1,p2,p3,getOrientation(P)))
        return false;

    for (size_t i = 0; i < P_size; i++) {
        if(p1 == P[i] || p2 == P[i] || p3 == P[i]){
            continue;
        }
        if(isInsideTriangle(p1,p2,p3,P[i]))
            return false;
    }
    return true;
}


double triangle_goodness (CGL::Point3 p0, CGL::Point3 p1, CGL::Point3 p2){
    double scale = 0;
    double angle0 = angleBetween(p0,p1,p0,p2) * 180/PI;
    double angle1 = angleBetween(p1,p0,p1,p2) * 180/PI;
    double angle2 = angleBetween(p2,p0,p2,p1) * 180/PI;

    double std_deviation = std::sqrt((angle0 - 60)*(angle0 - 60) + 
                                     (angle1 - 60)*(angle1 - 60) +
                                     (angle2 - 60)*(angle2 - 60)) / 3;

    // Considering that for std deviation 0 is ideal and 60 is horrible, 0 returns 1 and 60 returns 0

    // Scale of 0-1 of goodness:
    scale = 1 - (std_deviation / 60);
    return scale;
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

void edge_flip(CGL::Mesh::Halfedge_index& e, CGL::Mesh& polygon){
    CGL::Mesh::Face_index face = polygon.face(e);
    CGL::Mesh::Halfedge_index first_triagle_hedge = polygon.next(e);
    CGL::Mesh::Halfedge_index adjacent_triangle_hedge = polygon.next(polygon.opposite(e));

    CGAL::Euler::join_face(e,polygon);
    e = CGAL::Euler::split_face(first_triagle_hedge,adjacent_triangle_hedge,polygon);
}

// Is adjacente triangle better?
bool compare_goodness(CGL::Mesh::Halfedge_index first_triangle_hedge, CGL::Mesh::Halfedge_index adjacent_triangle_hedge, CGL::Mesh polygon){
    std::vector<CGL::Mesh::Vertex_index> actual_triangle_vertices = get_face_vertices(polygon.face(first_triangle_hedge),polygon);
    std::vector<CGL::Mesh::Vertex_index> adjacent_triangle_vertices = get_face_vertices(polygon.face(adjacent_triangle_hedge),polygon);

    double x1 = triangle_goodness(polygon.point(actual_triangle_vertices[0]),polygon.point(actual_triangle_vertices[1]),polygon.point(actual_triangle_vertices[2]));
    double x2 = triangle_goodness(polygon.point(adjacent_triangle_vertices[0]),polygon.point(adjacent_triangle_vertices[1]),polygon.point(adjacent_triangle_vertices[2]));
    double total_goodness_1 = x1 + x2;

    CGL::Point3 actual_ext_point = polygon.point(polygon.target(first_triangle_hedge));
    CGL::Point3 adjacent_ext_point = polygon.point(polygon.target(adjacent_triangle_hedge));

    x1 = triangle_goodness(actual_ext_point,adjacent_ext_point,polygon.point(polygon.source(first_triangle_hedge)));
    x2 = triangle_goodness(actual_ext_point,adjacent_ext_point,polygon.point(polygon.source(adjacent_triangle_hedge)));
    double total_goodness_2 = x1 + x2;

    if(total_goodness_2 > total_goodness_1) 
        return true;
    else        
        return false;
}

bool can_flip(CGL::Mesh::Halfedge_index initial_hedge, CGL::Mesh polygon){
    CGL::Mesh::Halfedge_index first_triangle_hedge = polygon.next(initial_hedge);
    CGL::Mesh::Halfedge_index adjacent_triangle_hedge = polygon.next(polygon.opposite(initial_hedge));

    CGAL::Euler::join_face(initial_hedge,polygon);

    if(!CGL::is_diagonal(polygon,first_triangle_hedge,adjacent_triangle_hedge))
        return false;
    
    return true;
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