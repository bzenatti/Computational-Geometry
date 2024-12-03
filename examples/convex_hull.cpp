/**
 * @file draw_points_2.cpp
 * @ingroup cglexamples
 * @brief Example of drawing 2D points.
 *
 * This example creates and draws a set of random 2D points.
 *
 * @author Ricardo Dutra da Silva
 */


#include <Types.h>
#include <DrawPoints2.h>
#include <CGAL/draw_polygon_2.h>
#include <Random.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

CGL::Polygon2 gift_wrapping(std::vector<CGL::Point2> pts);
bool left(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3);
float getDeterminant(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3);

/**
 * Code example.
 * @return Success.
 */
int main()
{

    /* Vector to store points. */
    std::vector<CGL::Point2> pts;
    CGL::Polygon2 P;
    /* Create pointss with random endpoints. */
    int ns = 7; /* Number of points. */
    CGL::Random rng(-10,10); /* "bounding box" for points. "*/
    for (int i = 1; i <= ns; i++)
    {
        /* Point. */
        CGL::Point2 p(rng.get_int(), rng.get_int());
        std::cout << "Point " << i << " (" << p << ")" << std::endl;

        /* Insert to list. */
        pts.push_back(p);
    }

    P = gift_wrapping(pts);

    /* Draw points. */
    CGL::draw(pts);
    CGAL::draw(P);

    return EXIT_SUCCESS;
}

CGL::Polygon2 gift_wrapping(std::vector<CGL::Point2> pts){
    CGL::Polygon2 hull;
    int leftmostIndex = 0;
    //In case of tie, lowest y is the leftmost
    for (int i = 1; i < pts.size(); i++)
    {
        if(pts[i].x() < pts[leftmostIndex].x()) {
            leftmostIndex = i;
        }
        else if(pts[i].x() == pts[leftmostIndex].x()) {
            if(pts[i].y() < pts[leftmostIndex].y())
                leftmostIndex = i;
        }
    }

    int p = leftmostIndex, q;
    do {
        hull.push_back(pts[p]);  // Add to the hull
        q = (p + 1) % pts.size(); // Initialize q to the next point

        for (int i = 0; i < pts.size(); i++)
            //If i is left pq, so q = i, because I want all points on the right
            if (left(pts[p], pts[q], pts[i])) 
                q = i;  
        p = q; 
        //CGAL::draw(hull);
    } while (p != leftmostIndex);  

    return hull;
}

//This function receives as parameter three points
bool left(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3){
    return (getDeterminant(p1,p2,p3) > 0);
}

//This function receives as parameter three points
float getDeterminant(CGL::Point2 p1, CGL::Point2 p2,CGL::Point2 p3){
    // Formula to get determinant with the points
    return ((p2.x() - p1.x())*(p3.y() - p1.y()) - (p2.y() - p1.y())*(p3.x() - p1.x()));
}


