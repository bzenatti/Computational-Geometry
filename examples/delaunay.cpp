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

bool compare_goodness(CGL::Mesh::Halfedge_index first_triangle_hedge, CGL::Mesh::Halfedge_index adjacent_triangle_hedge, CGL::Mesh polygon);
bool can_flip(CGL::Mesh::Halfedge_index initial_hedge, CGL::Mesh polygon);

std::vector<CGL::Point2>  createSuperTriangle(const std::vector<CGL::Point2>& points);
CGL::Mesh delaunay_triangulation(std::vector<CGL::Point2> mesh_vertices);
void insert_point(CGL::Mesh& mesh, CGL::Point2 vi);

int main(int argc, char* argv[]) {

    /* Vector to store points. */
    std::vector<CGL::Point3> pts;
    std::vector<CGL::Mesh::Vertex_index> indexes;

    /* Simple mesh with random vertices. */
    CGL::Mesh mesh_delaunay;

    int nv = 10;    /* Number of vertices. */

    std::vector<CGL::Point2>   point_vec;
    CGAL::Random         rand;

    rand.get_seed();

    
    //generate for representing 512x512 pixels
    CGAL::copy_n_unique(Point_generator(20), nv, std::back_inserter(point_vec));

    /*     
    for (const auto& point : point_set) {
        std::cout << point << std::endl;
    }
    */

    CGL::Mesh mesh;
    
    int i = 0 ;

    // Add vertices to the mesh
    for (auto p = point_vec.begin(); p != point_vec.end(); ++p) {
        CGL::Point3 p3(p->x(), p->y(), 0.0);  // Convert 2D point to 3D 
        CGL::Mesh::Vertex_index vi = mesh.add_vertex(p3);
        // if(i < 3)
        //     indexes.push_back(vi);
        // i++;
    }

    // CGAL::Euler::add_face(indexes,mesh);
    std::vector<CGL::Point2> superTriangle = createSuperTriangle(point_vec);
    
    for(int i = 0;i<3;i++){
        CGL::Point2 p = superTriangle[i];
        CGL::Point3 p3(p.x(), p.y(), 0.0);  // Convert 2D point to 3D 
        CGL::Mesh::Vertex_index vi = mesh.add_vertex(p3);
    }

    //Verify if all vertexes are inside triangle
    for(int i = 0;i<nv;i++){
        if(!isInsideTriangle(superTriangle[0],superTriangle[1],superTriangle[2],point_vec[i]))
            return -1;
        
    }

    CGL::print_mesh(mesh);

    //CGAL::draw(mesh);
   
    return EXIT_SUCCESS;
}

std::vector<CGL::Point2> createSuperTriangle(const std::vector<CGL::Point2>& points) {
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

    CGL::Point2 p1(midx - 2 * dmax, midy - dmax);
    CGL::Point2 p2(midx + 2 * dmax, midy - dmax);
    CGL::Point2 p3(midx, midy + 2 * dmax);

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

    for (auto p = vertices.begin(); p != vertices.end(); ++p) {
        insert_point(delaunay_mesh,*p);
    }
    


}

void insert_point(CGL::Mesh& mesh, CGL::Point2 vi){

    std::vector<CGL::Mesh::Face_index> bad_triangles;

    // Triangles that are no longer valid
    for (CGL::Mesh::Face_index f : mesh.faces()) {
        std::vector<CGL::Mesh::Vertex_index> triangle_vertices = get_face_vertices(f,mesh);

        CGL::Point3 point(vi.x(), vi.y(), 0.0);  // Convert 2D point to 3D 
        CGL::Point3 t1 = mesh.point(triangle_vertices[0]);
        CGL::Point3 t2 = mesh.point(triangle_vertices[1]);
        CGL::Point3 t3 = mesh.point(triangle_vertices[2]);

        if(inCircle(t1,t2,t3, point)){
            bad_triangles.push_back(f);
        }
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

