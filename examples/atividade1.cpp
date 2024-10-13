#include <Types.h>
#include <DrawPoints2.h>
#include <CGAL/draw_polygon_2.h>
#include <Random.h>
#include <iostream>
#include <stdlib.h>
#include <vector>


float getTriangleArea(std::vector<CGL::Point2> pts);

int main(){

    /* Vector to store points. */
    std::vector<CGL::Point2> pts;

    /* Simple polygon with random vertices. */
    CGL::Polygon2 P;

    //  Receive from the user 3 points
    for(int i = 0;i < 3;i++){
        printf("Insert the coordinates of three points separated by space: \n");
        int x1 = 0, y1=0;
        scanf("%d %d", &x1, &y1);

        CGL::Point2 p(x1, y1);
        pts.push_back(p);
        P.push_back(p);
    }

    // int x = pts[0].x();
    // printf("%d",x);

    /* Draw triangle. */
    CGAL::draw(P);

    float Parea = getTriangleArea(pts);
    printf("\nArea of the triangle : %.2f \n",Parea);

    if(Parea >= 0)
        printf("The orientation of the three points is counterclockwise\n");
    else
        printf("The orientation of the three points is clockwise\n");

    return 0;
}


float getTriangleArea(std::vector<CGL::Point2> pts){
    int x[3] = {0}, y[3] = {0};

    for(int i = 0;i < 3;i++){
        x[i] = pts[i].x();
        y[i] = pts[i].y();
    }
    // Formula to get triangle area with the points
    return 0.5*((x[1] - x[0])*(y[2] - y[0]) - (y[1] - y[0])*(x[2] - x[0]));
}