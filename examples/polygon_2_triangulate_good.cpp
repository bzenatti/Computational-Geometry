#include <Types.h>
#include <CGAL/draw_polygon_2.h>
#include <Random.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <CGAL/draw_surface_mesh.h>

#include <Partition.h>
#include <IO.h>

#define PI 3.14159265358979323846

// Computational Geometry project - Bruno Emanuel Zenatti

float getDeterminant(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3);
float getDeterminant(std::vector<CGL::Point2> pts);

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
std::vector<CGL::Mesh::Vertex_index> get_face_vertices(CGL::Mesh::Face_index face, CGL::Mesh polygon);


int main(int argc, char* argv[]) {

    /* Check input. */
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " input.obj" << std::endl;
        return EXIT_FAILURE;
    }

    CGL::Mesh original_polygon = CGL::read_mesh(argv[1]);

    CGL::Mesh polygon = original_polygon;

    CGL::polygon_2_triangulate_naive(polygon);   //triangulates polygon received as input

    //CGL::print_mesh(polygon);

    /*
        The first step of this project is to evaluate the "goodness" of triangles formed by the naive triangulation.
        To do this, we identify all triangles within the polygon by traversing its faces, as each face represents a triangle.
    */
    

    std::vector<std::vector<CGL::Mesh::Vertex_index>> polygon_triangles;

    // This main loop looks into all faces
    for (CGL::Mesh::Face_index f : polygon.faces()){

        std::vector<CGL::Mesh::Vertex_index> triangle_vertices = get_face_vertices(f,polygon);

        polygon_triangles.push_back(triangle_vertices);

        /*
            Important to know that the outer face is -1
        */
    }

    for (CGL::Mesh::Face_index f : polygon.faces()) {
        std::vector<CGL::Mesh::Vertex_index> triangle = get_face_vertices(f,polygon);
        double goodness = triangle_goodness(polygon.point(triangle[0]),polygon.point(triangle[1]),polygon.point(triangle[2]));
        printf("Goodness of triangle of face_idx %u = %f\n", f.idx(), goodness);
    }

    for (CGL::Mesh::Face_index f : polygon.faces()) {
        // This loop is for verifying the three adjacent triangles
        for(int i = 0;i < 3;i++){
            CGL::Mesh::Halfedge_index first_hedge = polygon.halfedge(f);
            double total_goodness_1 = 0;
            double total_goodness_2 = 0;

            // If the face of the twin is the outer one, skip it
            if(polygon.face(polygon.opposite(first_hedge)).idx() == -1)
                continue;

            std::vector<CGL::Mesh::Vertex_index> actual_triangle = get_face_vertices(f,polygon);
            std::vector<CGL::Mesh::Vertex_index> adjacent_triangle = get_face_vertices(polygon.face(polygon.opposite(first_hedge)),polygon);

            double x1 = triangle_goodness(polygon.point(actual_triangle[0]),polygon.point(actual_triangle[1]),polygon.point(actual_triangle[2]));
            double x2 = triangle_goodness(polygon.point(adjacent_triangle[0]),polygon.point(adjacent_triangle[1]),polygon.point(adjacent_triangle[2]));
            total_goodness_1 = x1 + x2;

            CGL::Point3 actual_ext_point;
            CGL::Point3 adjacent_ext_point;

            // gets all points of quadrilateral (externals of both triangles)
            for(int j = 0;j < 3;j++) {
                if(polygon.source(first_hedge) != actual_triangle[j] && polygon.source(polygon.opposite(first_hedge)) != actual_triangle[j])
                    actual_ext_point = polygon.point(actual_triangle[j]);
                if(polygon.source(first_hedge) != adjacent_triangle[j] && polygon.source(polygon.opposite(first_hedge)) != adjacent_triangle[j])
                    adjacent_ext_point = polygon.point(adjacent_triangle[j]);
            }

            //test if a segment between both external points intersects polygon
            //to use is_diagonal function, need the half-edge for getting its target
            // (the vertice that the hedge points to)

        }
    }
    /*
    for (CGL::Mesh::Face_index f : polygon.faces()) {
        CGL::Mesh::Halfedge_index hedge = polygon.halfedge(f);
        while(polygon.face(polygon.opposite(hedge)).idx() == -1) {
            hedge = polygon.next(polygon.opposite(hedge));
        }
        CGAL::Euler::join_face(hedge, polygon);
        break;
    }   
    */ 
    CGL::print_mesh(polygon);

    CGAL::draw(polygon);
    return EXIT_SUCCESS;
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

    double std_deviation = std::sqrt((angle0 - 60)*(angle0 - 60) + (angle1 - 60)*(angle1 - 60)+(angle0 - 60)*(angle0 - 60)) / 3;
    // double mean = (angle0 + angle1 + angle2)/3;
    // mean *= 180/PI;
    //printf("This should be, ideally, 0. \t %f  ", std_deviation);

    // Considering that for std deviation 0 is ideal and 60 is horrible, 0 returns 1 and 60 returns 0
    // Scale of 0-1 of goodness:
    scale = 1 - (std_deviation / 60);

    //printf("Scale of goodness %f \n", scale);

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
            printf("v: %5d x: %+10.2f  y: %+10.2f\n", triangle_vertices[i].idx(), pt.x(), pt.y());
        }
        printf("\n\n");
    
    return triangle_vertices;
}