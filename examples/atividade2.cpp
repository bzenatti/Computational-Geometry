#include <Types.h>
#include <CGAL/draw_polygon_2.h>
#include <Random.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>


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


int main() {
    /* Vector to store points. */
    std::vector<CGL::Point2> pts;
    
    /* Simple polygon with random vertices. */
    CGL::Polygon2 P;

    /* Create vertices. */
    int ns = 6; /* Number of vertices. */
    CGL::Random rng(-10,10); /* "bounding box" for points. "*/

    int i = 0;

    // Generates a random polygon until is met the condition
    do 
    {
        P.clear();
        pts.clear();
        for (i = 0; i < ns; i++){
            CGL::Point2 p(rng.get_int(),rng.get_int());
            pts.push_back(p);
            P.push_back(p);
        }
    }
    while (!P.is_simple());

    bool orientation = getOrientation(P);

    if(orientation) std::cout << "The orientation is "<< "Counterclockwise" << std::endl;
    else std::cout << "The orientation is " << "Clockwise" << std::endl;

    size_t P_size = P.size();

    for (size_t i = P_size; i < 2*P.size(); i++){
        int actualVertice = i % P_size;
        int prevVertice = (i-1) % P_size;
        int nextVertice = (i+1) % P_size;

        if(isVerticeConvex(P[prevVertice],P[actualVertice],P[nextVertice],orientation))
            std::cout << "The vertice ("<< P[actualVertice] << ") "<< "is convex"<<std::endl;
        else 
            std::cout << "The vertice ("<< P[actualVertice] << ") "<< "is concave"<<std::endl;
    }

    for (size_t i = P_size; i < 2*P_size; i++){
        int actualVertice = i % P_size;
        int prevVertice = (i-1) % P_size;
        int nextVertice = (i+1) % P_size;
        bool isEar = true;
        int count = 0;

        isEar = isVerticeEar(P[prevVertice],P[actualVertice],P[nextVertice],P);

        if(isEar)
            std::cout << "The vertice ("<< P[actualVertice] << ") "<< "is ear"<<std::endl;
        else 
            std::cout << "The vertice ("<< P[actualVertice] << ") "<< "is not ear"<<std::endl;
    }


    /* Draw points. */
    CGAL::draw(P);

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
double dotProduct(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3,CGL::Point2 p4) {
    int ux = p2.x() - p1.x(), uy = p2.y() - p1.y();
    int vx = p4.x() - p3.x(), vy = p4.y() - p3.y();

    return (ux*vx) + (uy+vy);
}

double vectorMagnitude(CGL::Point2 p1, CGL::Point2 p2){
    return std::sqrt(dotProduct(p1,p2,p1,p2));
}

//returns in radians the angle between two vectors, p1p2==u and p3p4==v
double angleBetween(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3,CGL::Point2 p4){
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