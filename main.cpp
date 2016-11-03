#include "linalg.h"
#include "bbfmm.h"
#include "mapping.h"
#include "quadratue.h"
#include "gmres.h"

using namespace bbfmm;

int main() {
#ifdef RUN_OMP
    omp_set_num_threads(omp_get_max_threads());
#endif

    vector<point> source, target;
    vector<scalar_t> weight;
    vector<point> triangle;

    stdIcosahedronMapping(16, source, target, weight, triangle);

    assert(source.size() == weight.size());

    std::cout << "source : " << source.size() << std::endl;
    std::cout << "target : " << target.size() << std::endl;

    vector<point> centers;
    centers.push_back(point(0., 0., 0.));
    vector<scalar_t> radii;
    radii.push_back(1.0);

    kernel k0, k1;

    scalar_t  dI = 2.0;
    scalar_t  dE = 80.0;

    Vector p(target.size());
    for (int i = 0; i < target.size(); ++i) {
        p(i) = weight[i];
    }

    k0.initialize(8, source, target, p, source.size(), target.size(), 640, 10);
    k1.initialize(8, source, target, p, source.size(), target.size(), 640, 10);

    auto G0 = [&](point &a, point &b) {
        scalar_t r = norm(a, b);
        return 1.0 / 4 / M_PI / r;
    };
    auto G1 = [&](point &a, point &b) {
        scalar_t r = norm(a, b);
        return -(a.x * (a.x - b.x) + a.y * (a.y - b.y) + a.z * (a.z - b.z)) / 4 / M_PI / r / r / r;
    };

    k0.eval = [&](point &a, point &b) {
        if (a == b) {
            return 0.;
        }
        else {
            return G0(a, b);
        }

    };

    k1.eval = [&](point &a, point &b) {
        if (a == b) {
            return 0.;
        }
        else {
            return G1(a, b);
        }

    };


    auto A = [&](Vector &f) {
        Vector q;
        Vector v;
        int half = f.row()/ 2;

        Vector f0(half);
        Vector f1(half);

        Vector ret(f.row());

        Vector ff(half * 4);
        Vector gg(half * 4);
        for (int i = 0; i < half; ++i) {
            ff(4 * i) = f(i) * weight[4 * i];
            ff(4 * i + 1) = f(i) * weight[4 * i + 1];
            ff(4 * i + 2) = f(i) * weight[4 * i + 2];
            ff(4 * i + 3) = f(i) * weight[4 * i + 3];
            gg(4 * i) = f(i + half) * weight[4 * i];
            gg(4 * i + 1) = f(i + half) * weight[4 * i + 1];
            gg(4 * i + 2) = f(i + half) * weight[4 * i + 2];
            gg(4 * i + 3) = f(i + half) * weight[4 * i + 3];
            f0(i) = f(i);
            f1(i) = f(i + half);
        }
        k0.chargeTree = ff; k1.chargeTree = gg;
        k0.reset(); k1.reset();
        k0.run(q); k1.run(v);

        for (int i = 0; i < half; ++i) {
            ret(i) = 0.5 * f(i + half) + v(i) - q(i);
            ret(i + half) = 0.5 * f(i + half) - v(i) + dI/dE * q(i);
        }

        return ret;
    };


    Vector x(target.size() * 2);
    Vector b(target.size() * 2);
    Vector y(target.size() * 2);
    for (int i = 0; i < target.size(); ++i) {
        b(i) = 1.0/4.0/M_PI/dI;
        y(i) = -1.0/4.0/M_PI/dI;
        y(i + target.size()) = 1.0/4.0/M_PI/dE;
    }
//    setValue(b, 1.0/ 4.0/M_PI / dI);

    GMRES(A, x, b, 10, 40, 1e-5);

    scalar_t ret = 0., ret_ = 0.;
    for (int i = 0; i < x.row(); ++i) {
        ret += (x(i) - y(i)) * (x(i) - y(i));
        ret_ += y(i) * y(i);
    }
    std::cout << sqrt(ret) / sqrt(ret_) << std::endl;
}
