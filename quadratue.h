/*
 *
 *  Copyright (C) 2016, Yimin Zhong, University of Texas at Austin.
 *
 *  Quadrature rules for singular integral
 *
 */

#ifndef MOLECULE_QUADRATUE_H
#define MOLECULE_QUADRATUE_H

#include "mapping.h"

inline scalar_t singularAngle(point& a, point &b, int id, vector<point>& centers) {
    return angle(a, b, centers[id]);
}

inline void gauss(point &singularPoint, point &left, point &right, int id, vector<scalar_t> &radii,
                  vector<point> &centers, vector<point> &quadraturePoint, vector<scalar_t> &quadratureWeight) {

    vector<scalar_t> gaussPoint = {(1 - std::sqrt(3.) / sqrt(5.)) / 2., 0.5, (1 + std::sqrt(3.) / sqrt(5.)) / 2.};
    vector<scalar_t> gaussWeight = {5. / 18., 4. / 9., 5. / 18.};

    int N = (int) gaussPoint.size();

    point pre_projection;

    for (int i = 0; i < N; ++i) {
        double curWeight = gaussWeight[i];
        double curPoint1 = gaussPoint[i];
        for (int j = 0; j < N; ++j) {
            double combinedWeight = gaussWeight[j] * curWeight;
            /*
             * rescale due to singularity
             */
            double curPoint2 = gaussPoint[j] * (1.0 - curPoint1);
            quadratureWeight.push_back(2.0 * (1.0 - curPoint1) * combinedWeight);

            /*
             * Cartesian coordinate of quadrature point
             */

            pre_projection.x =
                    singularPoint.x * curPoint1 + right.x * curPoint2 + left.x * (1.0 - curPoint1 - curPoint2);
            pre_projection.y =
                    singularPoint.y * curPoint1 + right.y * curPoint2 + left.y * (1.0 - curPoint1 - curPoint2);
            pre_projection.z =
                    singularPoint.z * curPoint1 + right.z * curPoint2 + left.z * (1.0 - curPoint1 - curPoint2);

            /*
             * project to sphere
             */
            double scale = radii[id] / norm(pre_projection, centers[id]);
            quadraturePoint.push_back(point(
                    (pre_projection.x - centers[id].x) * scale + centers[id].x,
                    (pre_projection.y - centers[id].y) * scale + centers[id].y,
                    (pre_projection.z - centers[id].z) * scale + centers[id].z));
        }
    }
}


inline double
singularIntegral(point &singular, vector<point> &vertices, int id, vector<scalar_t> &radii, vector<point> &centers,
                 std::function<double(point &, point &)> eval) noexcept {
    /*
     * get 3 vertices of triangle enclosing singularity
     */
    point &pointA = vertices[3 * singular.triangleId];
    point &pointB = vertices[3 * singular.triangleId + 1];
    point &pointC = vertices[3 * singular.triangleId + 2];

    double angleA = singularAngle(pointB, pointC, id, centers);
    double angleB = singularAngle(pointC, pointA, id, centers);
    double angleC = singularAngle(pointA, pointB, id, centers);
    double triangleArea = area(angleA, angleB, angleC, radii[id]);

    double angleSubA, angleSubB, angleSubC, triangleSubArea;

    angleSubA = singularAngle(pointA, singular, id, centers);
    angleSubB = singularAngle(pointB, singular, id, centers);
    angleSubC = singularAngle(pointC, singular, id, centers);

    vector<point> quadraturePoint;
    vector<double> quadratureWeight;

    double ret = 0.;

    /*
     * AB Part
     */
    gauss(singular, pointA, pointB, id, radii, centers, quadraturePoint, quadratureWeight);
    triangleSubArea = area(angleC, angleSubA, angleSubB, radii[id]);

    for (int i = 0; i < quadraturePoint.size(); ++i) {
        ret += quadratureWeight[i] * triangleSubArea / triangleArea * eval(quadraturePoint[i], singular);
    }

    quadraturePoint.clear();
    quadratureWeight.clear();

    /*
     * BC part
     */
    gauss(singular, pointB, pointC, id, radii, centers, quadraturePoint, quadratureWeight);
    triangleSubArea = area(angleA, angleSubB, angleSubC, radii[id]);

    for (int i = 0; i < quadraturePoint.size(); ++i) {
        ret += quadratureWeight[i] * triangleSubArea / triangleArea * eval(quadraturePoint[i], singular);
    }
    quadraturePoint.clear();
    quadratureWeight.clear();


    /*
     * CA part
     */
    gauss(singular, pointC, pointA, id, radii, centers, quadraturePoint, quadratureWeight);
    triangleSubArea = area(angleB, angleSubC, angleSubA, radii[id]);

    for (int i = 0; i < quadraturePoint.size(); ++i) {
        ret += quadratureWeight[i] * triangleSubArea / triangleArea * eval(quadraturePoint[i], singular);
    }
    quadraturePoint.clear();
    quadratureWeight.clear();

    return ret;

}


inline double nearFieldIntegral(point& source, point& target, vector<point>& vertices, int id, vector<scalar_t>
&radii, vector<point> &centers, std::function<double(point &, point &)> eval) {

    point &pointA = vertices[3 * source.triangleId];
    point &pointB = vertices[3 * source.triangleId + 1];
    point &pointC = vertices[3 * source.triangleId + 2];

    // all at leaf and not same triangles.
    assert(source.triangleId >= 0);

    scalar_t dist;
    scalar_t scale;

    point pc((pointC.x  -centers[id].x) * 2.0/3.0 +(pointA.x - centers[id].x) / 6.0 + (pointB.x  - centers[id].x)/ 6.0,
             (pointC.y  -centers[id].y) * 2.0/3.0 +(pointA.y - centers[id].y) / 6.0 + (pointB.y  - centers[id].y)/ 6.0,
             (pointC.z  -centers[id].z) * 2.0/3.0 +(pointA.z - centers[id].z) / 6.0 + (pointB.z  - centers[id].z)/ 6.0);

    dist = norm(pc, centers[id]);
    scale = radii[id]/dist;

    pc.x = pc.x * scale + centers[id].x;
    pc.y = pc.y * scale + centers[id].y;
    pc.z = pc.z * scale + centers[id].z;

    point pb((pointB.x  -centers[id].x) * 2.0/3.0 +(pointA.x - centers[id].x) / 6.0 + (pointC.x  - centers[id].x)/ 6.0,
             (pointB.y  -centers[id].y) * 2.0/3.0 +(pointA.y - centers[id].y) / 6.0 + (pointC.y  - centers[id].y)/ 6.0,
             (pointB.z  -centers[id].z) * 2.0/3.0 +(pointA.z - centers[id].z) / 6.0 + (pointC.z  - centers[id].z)/ 6.0);

    dist = norm(pb, centers[id]);
    scale = radii[id]/dist;

    pb.x = pb.x * scale + centers[id].x;
    pb.y = pb.y * scale + centers[id].y;
    pb.z = pb.z * scale + centers[id].z;

    point pa((pointA.x  -centers[id].x) * 2.0/3.0 +(pointC.x - centers[id].x) / 6.0 + (pointB.x  - centers[id].x)/ 6.0,
             (pointA.y  -centers[id].y) * 2.0/3.0 +(pointC.y - centers[id].y) / 6.0 + (pointB.y  - centers[id].y)/ 6.0,
             (pointA.z  -centers[id].z) * 2.0/3.0 +(pointC.z - centers[id].z) / 6.0 + (pointB.z  - centers[id].z)/ 6.0);

    dist = norm(pa, centers[id]);
    scale = radii[id]/dist;

    pa.x = pa.x * scale + centers[id].x;
    pa.y = pa.y * scale + centers[id].y;
    pa.z = pa.z * scale + centers[id].z;

    return 0.5 * (eval(pa, target) + eval(pb, target) + eval(pc, target) - eval(source, target));
}
#endif //MOLECULE_QUADRATUE_H
