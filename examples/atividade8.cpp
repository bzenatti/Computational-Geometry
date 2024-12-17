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
#define MAXN 100

// Computational Geometry project - Bruno Emanuel Zenatti

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


int main(int argc, char* argv[]) {

    /* Vector to store points. */
    std::vector<CGL::Point3> pts;
    std::vector<CGL::Mesh::Vertex_index> indexes;

    /* Simple polygon with random vertices. */
    CGL::Mesh totalPolygon;

    unsigned char c = 'a';

    int predefined_points[4][2] = {
        {2, 2}, // Point A
        {4, 0}, // Point B
        {4, 4}, // Point C
        {6, 2}  // Point D
    };

    for(int i = 0;i < 4;i++){
        int x1 = predefined_points[i][0];
        int y1 = predefined_points[i][1];

        CGL::Point3 p(x1, y1,0);

        pts.push_back(p);
        totalPolygon;
        indexes.push_back(totalPolygon.add_vertex(p));
    }
    // //  Uncomment for reading from user
    // //  Receive from the user 4 points
    // for(int i = 0;i < 4;i++){
    //     printf("Insert the coordinates of point %c separated by space:\n", c);
    //     c++;
    //     int x1 = 0, y1=0;
    //     scanf("%d %d", &x1, &y1);

    //     CGL::Point3 p(x1, y1,0);
    //     pts.push_back(p);
    //     totalPolygon;
    //     indexes.push_back(totalPolygon.add_vertex(p));
    // }

    totalPolygon.add_face(indexes[0],indexes[1],indexes[2]);
    totalPolygon.add_face(indexes[2],indexes[1],indexes[3]);

    if(inCircle(pts[0],pts[1],pts[2],pts[3])){
        printf("\nPonto d está dentro do circulo formado pelo triangulo abc\n");
    }
    else{
        printf("\nPonto d nao esta dentro do circulo formado pelo triangulo abc\n");
    }

    if(inCircle(pts[2],pts[1],pts[3],pts[0])){
        printf("\nPonto a está dentro do circulo formado pelo triangulo cbd\n");
    }
    else{
        printf("\nPonto a nao esta dentro do circulo formado pelo triangulo cbd\n");
    }

    CGL::polygon_2_triangulate_naive(totalPolygon);

    CGAL::draw(totalPolygon);

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

bool inCircle(CGL::Point3 a, CGL::Point3 b,CGL::Point3 c,CGL::Point3 d){
    double matrix[4][4] = {{a.x(),a.y(), (a.x()*a.x()) + (a.y()*a.y()),1},
                           {b.x(),b.y(), (b.x()*b.x()) + (b.y()*b.y()),1},
                           {c.x(),c.y(), (c.x()*c.x()) + (c.y()*c.y()),1},
                           {d.x(),d.y(), (d.x()*d.x()) + (d.y()*d.y()),1}};

    //Includind if the point is in the circumference
    if(getDeterminant4x4(matrix, 4) >= 0)
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