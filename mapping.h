/*
 *
 *  Copyright (C) 2016, Yimin Zhong, University of Texas at Austin.
 *
 *  Mapping ico-triangles onto spherical triangles.
 *
 *
 */
#ifndef MOLECULE_MAPPING_H
#define MOLECULE_MAPPING_H

#include "utils.h"
#include "bbfmm.h"

using namespace bbfmm;


using std::sqrt;
using std::acos;
using std::atan;


inline scalar_t inner(point &a, point &b, point &center) {
    return (a.x - center.x) * (b.x - center.x) + (a.y - center.y) * (b.y - center.y) +
           (a.z - center.z) * (b.z - center.z);
}

inline scalar_t norm(point &a, point &center) {
    return sqrt(inner(a, a, center));
}

inline scalar_t angle(point &a, point &b, point &center) {
    auto ra = norm(a, center);
    auto rb = norm(b, center);
    return acos(inner(a, b, center) / ra / rb);
}

inline scalar_t area(scalar_t a, scalar_t b, scalar_t c, scalar_t sRadius = 1.0) {
    scalar_t s = (a + b + c) * 0.5;
    return 4.0 * sRadius * sRadius *
           atan(sqrt(tan(0.5 * s) * tan(0.5 * (s - a)) * tan(0.5 * (s - b)) * tan(0.5 * (s - c))));
}

inline void quadrature(vector<point> &source, vector<scalar_t> weight, scalar_t sRadius, scalar_t triangleArea) {

}

inline void stdIcosahedron(vector<point> &vertices, vector<vector<int>> &faces, vector<vector<int>> &edges) {
    /*
     * define standard icosahedron
     */
    if (vertices.size() != 12) { vertices.resize(12); }
    if (faces.size() != 20) {
        faces.resize(20);
        for (int i = 0; i < 20; ++i) { faces[i].resize(3); }
    }
    if (edges.size() != 30) {
        edges.resize(30);
        for (int i= 0; i < 30; ++i) { edges[i].resize(2);}
    }

    scalar_t phi = 0.5 * (1.0 + sqrt(5.0));

    for (int i = 0; i < 4; ++i) {
        vertices[i].x = 0.;
        vertices[i].y = (i & 1) * 2 - 1;
        vertices[i].z = phi * (((i >> 1) & 1) * 2 - 1);
    }

    for (int i = 4; i < 8; ++i) {
        vertices[i].y = 0.;
        vertices[i].z = (i & 1) * 2 - 1;
        vertices[i].x = phi * (((i >> 1) & 1) * 2 - 1);
    }

    for (int i = 8; i < 12; ++i) {
        vertices[i].z = 0.;
        vertices[i].x = (i & 1) * 2 - 1;
        vertices[i].y = phi * (((i >> 1) & 1) * 2 - 1);
    }

    faces[0] = {3, 10, 11};
    faces[1] = {3, 7, 11};
    faces[2] = {3, 10, 5};
    faces[3] = {3, 2, 5};
    faces[4] = {3, 2, 7};
    faces[5] = {2, 5, 8};
    faces[6] = {2, 8, 9};
    faces[7] = {2, 9, 7};
    faces[8] = {7, 6, 9};
    faces[9] = {7, 6, 11};
    faces[10] = {5, 4, 10};
    faces[11] = {5, 4, 8};
    faces[12] = {0, 8, 9};
    faces[13] = {0, 9, 6};
    faces[14] = {0, 8, 4};
    faces[15] = {0, 1, 4};
    faces[16] = {0, 1, 6};
    faces[17] = {1, 4, 10};
    faces[18] = {1, 10, 11};
    faces[19] = {1, 11, 6};

    unordered_set<int> edge_set;
    for (int i = 0; i < faces.size(); ++i) {
        if (faces[i][0] < faces[i][1]) {
            edge_set.insert(faces[i][0] * 12 + faces[i][1]);
        }
        else {
            edge_set.insert(faces[i][1] * 12 + faces[i][0]);
        }
        if (faces[i][1] < faces[i][2]) {
            edge_set.insert(faces[i][1] * 12 + faces[i][2]);
        }
        else {
            edge_set.insert(faces[i][2] * 12 + faces[i][1]);
        }
        if (faces[i][2] < faces[i][0]) {
            edge_set.insert(faces[i][2] * 12 + faces[i][0]);
        }
        else {
            edge_set.insert(faces[i][0] * 12 + faces[i][2]);
        }
    }

    assert(edge_set.size() == 30);

    int i = 0;
    for (auto e : edge_set) {
        edges[i][0] = e / 12;
        edges[i][1] =  e % 12;
        ++i;
    }
}

/*
 * for a single unit ball
 */
inline void stdIcosahedronMapping(int N, vector<point> &sources, vector<point>& targets, vector<scalar_t> &weight,
                                  vector<point> &triangle) {
    vector<point> vertices;
    vector<vector<int>> faces;
    vector<vector<int>> edges;
    stdIcosahedron(vertices, faces, edges);
    /*
     * subdivide each triangle into small equilateral triangles.
     */
    point center(0., 0., 0.);
    scalar_t delta = 1.0 / scalar_t(N);
    scalar_t alpha, beta, gamma;
    scalar_t lambda = 1.0 / 3.0, mu = 1.0 / 3.0;
    point sub_a, sub_b, sub_c, sub_center, projection;
    point sub_center_a, sub_center_b, sub_center_c;
    point projection_a, projection_b, projection_c;
    scalar_t arc_a, arc_b, arc_c, _area;
    scalar_t _area_a, _area_b, _area_c;

    int curId = 0;
    /*
     * iterate over each faces.
     */
    for (int i = 0; i < 20; ++i) {
        point &a = vertices[faces[i][0]];
        point &b = vertices[faces[i][1]];
        point &c = vertices[faces[i][2]];

        for (int j = 0; j < N; ++j) {
            beta = j * delta;
            for (int k = 0; k < N - j; ++k) {
                alpha = k * delta;
                gamma = 1.0 - beta - alpha;

                sub_c.x = alpha * a.x + beta * b.x + gamma * c.x;
                sub_c.y = alpha * a.y + beta * b.y + gamma * c.y;
                sub_c.z = alpha * a.z + beta * b.z + gamma * c.z;

                sub_a.x = sub_c.x + delta * (a.x - c.x);
                sub_a.y = sub_c.y + delta * (a.y - c.y);
                sub_a.z = sub_c.z + delta * (a.z - c.z);

                sub_b.x = sub_c.x + delta * (b.x - c.x);
                sub_b.y = sub_c.y + delta * (b.y - c.y);
                sub_b.z = sub_c.z + delta * (b.z - c.z);

                sub_center.x = lambda * sub_a.x + mu * sub_b.x + (1 - lambda - mu) * sub_c.x;
                sub_center.y = lambda * sub_a.y + mu * sub_b.y + (1 - lambda - mu) * sub_c.y;
                sub_center.z = lambda * sub_a.z + mu * sub_b.z + (1 - lambda - mu) * sub_c.z;


                sub_center_a.x = 2.0/3.0 * sub_a.x + 1.0/6.0 * sub_b.x + 1.0/6.0 * sub_c.x;
                sub_center_a.y = 2.0/3.0 * sub_a.y + 1.0/6.0 * sub_b.y + 1.0/6.0 * sub_c.y;
                sub_center_a.z = 2.0/3.0 * sub_a.z + 1.0/6.0 * sub_b.z + 1.0/6.0 * sub_c.z;

                sub_center_b.x = 2.0/3.0 * sub_b.x + 1.0/6.0 * sub_a.x + 1.0/6.0 * sub_c.x;
                sub_center_b.y = 2.0/3.0 * sub_b.y + 1.0/6.0 * sub_a.y + 1.0/6.0 * sub_c.y;
                sub_center_b.z = 2.0/3.0 * sub_b.z + 1.0/6.0 * sub_a.z + 1.0/6.0 * sub_c.z;

                sub_center_c.x = 2.0/3.0 * sub_c.x + 1.0/6.0 * sub_b.x + 1.0/6.0 * sub_a.x;
                sub_center_c.y = 2.0/3.0 * sub_c.y + 1.0/6.0 * sub_b.y + 1.0/6.0 * sub_a.y;
                sub_center_c.z = 2.0/3.0 * sub_c.z + 1.0/6.0 * sub_b.z + 1.0/6.0 * sub_a.z;

                /*
                 * projection to unit sphere.
                 */
                scalar_t dist = norm(sub_center, center);
                projection = {sub_center.x / dist, sub_center.y / dist, sub_center.z / dist};

                targets.push_back(projection);
                targets[targets.size() - 1].triangleId = curId;

                sources.push_back(projection);
                sources[sources.size() - 1].triangleId = curId;

                dist = norm(sub_center_a, center);
                projection_a = {sub_center_a.x / dist, sub_center_a.y / dist, sub_center_a.z / dist};
                sources.push_back(projection_a);
                sources[sources.size() - 1].triangleId = curId;


                dist = norm(sub_center_b, center);
                projection_b = {sub_center_b.x / dist, sub_center_b.y / dist, sub_center_b.z / dist};
                sources.push_back(projection_b);
                sources[sources.size() - 1].triangleId = curId;

                dist = norm(sub_center_c, center);
                projection_c = {sub_center_c.x / dist, sub_center_c.y / dist, sub_center_c.z / dist};
                sources.push_back(projection_c);
                sources[sources.size() - 1].triangleId = curId;



                curId++;
                triangle.push_back(sub_a);
                triangle.push_back(sub_b);
                triangle.push_back(sub_c);


                arc_a = angle(sub_b, sub_c, center);
                arc_b = angle(sub_c, sub_a, center);
                arc_c = angle(sub_a, sub_b, center);
                _area = area(arc_a, arc_b, arc_c);

                point mid_ab(0.5 * (sub_a.x + sub_b.x), 0.5 * (sub_a.y + sub_b.y), 0.5 * (sub_a.z + sub_b.z));
                point mid_bc(0.5 * (sub_c.x + sub_b.x), 0.5 * (sub_c.y + sub_b.y), 0.5 * (sub_c.z + sub_b.z));
                point mid_ca(0.5 * (sub_a.x + sub_c.x), 0.5 * (sub_a.y + sub_c.y), 0.5 * (sub_a.z + sub_c.z));

                arc_a = angle(sub_a, mid_ab, center);
                arc_b = angle(sub_a, mid_ca, center);
                arc_c = angle(mid_ab, mid_ca, center);

                _area_a = area(arc_a, arc_b, arc_c);

                arc_a = angle(sub_b, mid_ab, center);
                arc_b = angle(sub_b, mid_bc, center);
                arc_c = angle(mid_ab, mid_bc, center);

                _area_b = area(arc_a, arc_b, arc_c);


                arc_a = angle(sub_c, mid_bc, center);
                arc_b = angle(sub_c, mid_ca, center);
                arc_c = angle(mid_ca, mid_bc, center);

                _area_c = area(arc_a, arc_b, arc_c);

                /*
                 * todo: update weight
                 */
                weight.push_back(_area - 2 * _area_a - 2 * _area_b - 2 * _area_c);
                weight.push_back(2 * _area_a);
                weight.push_back(2 * _area_b);
                weight.push_back(2 * _area_c);
            }
        }

        for (int j = 1; j < N; ++j) {
            beta = j * delta;
            for (int k = 1; k <= N - j; ++k) {
                alpha = k * delta;
                gamma = 1.0 - beta - alpha;

                sub_c.x = alpha * a.x + beta * b.x + gamma * c.x;
                sub_c.y = alpha * a.y + beta * b.y + gamma * c.y;
                sub_c.z = alpha * a.z + beta * b.z + gamma * c.z;

                sub_a.x = sub_c.x - delta * (a.x - c.x);
                sub_a.y = sub_c.y - delta * (a.y - c.y);
                sub_a.z = sub_c.z - delta * (a.z - c.z);

                sub_b.x = sub_c.x - delta * (b.x - c.x);
                sub_b.y = sub_c.y - delta * (b.y - c.y);
                sub_b.z = sub_c.z - delta * (b.z - c.z);

                sub_center.x = lambda * sub_a.x + mu * sub_b.x + (1 - lambda - mu) * sub_c.x;
                sub_center.y = lambda * sub_a.y + mu * sub_b.y + (1 - lambda - mu) * sub_c.y;
                sub_center.z = lambda * sub_a.z + mu * sub_b.z + (1 - lambda - mu) * sub_c.z;


                sub_center_a.x = 2.0/3.0 * sub_a.x + 1.0/6.0 * sub_b.x + 1.0/6.0 * sub_c.x;
                sub_center_a.y = 2.0/3.0 * sub_a.y + 1.0/6.0 * sub_b.y + 1.0/6.0 * sub_c.y;
                sub_center_a.z = 2.0/3.0 * sub_a.z + 1.0/6.0 * sub_b.z + 1.0/6.0 * sub_c.z;

                sub_center_b.x = 2.0/3.0 * sub_b.x + 1.0/6.0 * sub_a.x + 1.0/6.0 * sub_c.x;
                sub_center_b.y = 2.0/3.0 * sub_b.y + 1.0/6.0 * sub_a.y + 1.0/6.0 * sub_c.y;
                sub_center_b.z = 2.0/3.0 * sub_b.z + 1.0/6.0 * sub_a.z + 1.0/6.0 * sub_c.z;

                sub_center_c.x = 2.0/3.0 * sub_c.x + 1.0/6.0 * sub_b.x + 1.0/6.0 * sub_a.x;
                sub_center_c.y = 2.0/3.0 * sub_c.y + 1.0/6.0 * sub_b.y + 1.0/6.0 * sub_a.y;
                sub_center_c.z = 2.0/3.0 * sub_c.z + 1.0/6.0 * sub_b.z + 1.0/6.0 * sub_a.z;

                /*
                 * projection to unit sphere.
                 */
                scalar_t dist = norm(sub_center, center);
                projection = {sub_center.x / dist, sub_center.y / dist, sub_center.z / dist};

                targets.push_back(projection);
                targets[targets.size() - 1].triangleId = curId;

                sources.push_back(projection);
                sources[sources.size() - 1].triangleId = curId;

                dist = norm(sub_center_a, center);
                projection_a = {sub_center_a.x / dist, sub_center_a.y / dist, sub_center_a.z / dist};
                sources.push_back(projection_a);
                sources[sources.size() - 1].triangleId = curId;


                dist = norm(sub_center_b, center);
                projection_b = {sub_center_b.x / dist, sub_center_b.y / dist, sub_center_b.z / dist};
                sources.push_back(projection_b);
                sources[sources.size() - 1].triangleId = curId;

                dist = norm(sub_center_c, center);
                projection_c = {sub_center_c.x / dist, sub_center_c.y / dist, sub_center_c.z / dist};
                sources.push_back(projection_c);
                sources[sources.size() - 1].triangleId = curId;

                curId++;
                triangle.push_back(sub_a);
                triangle.push_back(sub_b);
                triangle.push_back(sub_c);


                arc_a = angle(sub_b, sub_c, center);
                arc_b = angle(sub_c, sub_a, center);
                arc_c = angle(sub_a, sub_b, center);
                _area = area(arc_a, arc_b, arc_c);



                point mid_ab(0.5 * (sub_a.x + sub_b.x), 0.5 * (sub_a.y + sub_b.y), 0.5 * (sub_a.z + sub_b.z));
                point mid_bc(0.5 * (sub_c.x + sub_b.x), 0.5 * (sub_c.y + sub_b.y), 0.5 * (sub_c.z + sub_b.z));
                point mid_ca(0.5 * (sub_a.x + sub_c.x), 0.5 * (sub_a.y + sub_c.y), 0.5 * (sub_a.z + sub_c.z));

                arc_a = angle(sub_a, mid_ab, center);
                arc_b = angle(sub_a, mid_ca, center);
                arc_c = angle(mid_ab, mid_ca, center);

                _area_a = area(arc_a, arc_b, arc_c);

                arc_a = angle(sub_b, mid_ab, center);
                arc_b = angle(sub_b, mid_bc, center);
                arc_c = angle(mid_ab, mid_bc, center);

                _area_b = area(arc_a, arc_b, arc_c);


                arc_a = angle(sub_c, mid_bc, center);
                arc_b = angle(sub_c, mid_ca, center);
                arc_c = angle(mid_ca, mid_bc, center);

                _area_c = area(arc_a, arc_b, arc_c);



                /*
                 * todo: update weight
                 */
                weight.push_back(_area - 2 * _area_a - 2 * _area_b - 2 * _area_c);
                weight.push_back(2 * _area_a);
                weight.push_back(2 * _area_b);
                weight.push_back(2 * _area_c);
            }
        }
    }
}


#endif //MOLECULE_MAPPING_H
