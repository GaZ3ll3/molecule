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

    kernel k;

    Vector p(target.size());
    for (int i = 0; i < target.size(); ++i) {
        p(i) = weight[i];
    }

    k.initialize(6, source, target, p, source.size(), target.size(), 320, 10);

    auto G0 = [&](point &a, point &b) {
        scalar_t r = norm(a, b);
        return (a.x * (a.x - b.x) + a.y * (a.y - b.y) + a.z * (a.z - b.z)) / 4 / M_PI / r / r / r;
    };

    k.eval = [&](point &a, point &b) {
        if (a == b) {
            return 0.;
        }
        else {
            return G0(a, b);
        }

    };


    auto A = [&](Vector &f) {
        Vector q;
        Vector ff(f.row() * 4);
        for (int i = 0; i < f.row(); ++i) {
            ff(4 * i) = f(i) * weight[4 * i];
            ff(4 * i + 1) = f(i) * weight[4 * i + 1];
            ff(4 * i + 2) = f(i) * weight[4 * i + 2];
            ff(4 * i + 3) = f(i) * weight[4 * i + 3];
        }
        k.chargeTree = ff;
        k.reset();
        k.run(q);
        daxpy(0.5, f, q);
        return q;
    };


    Vector x(target.size());
    Vector b(target.size());
    setValue(b, 1.0);

    GMRES(A, x, b, 10, 40, 1e-9);

    scalar_t ret = 0.;
    for (int i = 0; i < x.row(); ++i) {
        ret += (x(i) - 1.0) * (x(i) - 1.0);
    }
    std::cout << sqrt(ret) / sqrt(x.row()) << std::endl;
}
