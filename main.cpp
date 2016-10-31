#include "linalg.h"
#include "bbfmm.h"
#include "mapping.h"

using namespace bbfmm;

int main() {
#ifdef RUN_OMP
    omp_set_num_threads(omp_get_max_threads());
#endif

    vector<point> source;
    vector<scalar_t> weight;

    stdIcosahedronMapping(32, source, weight);


    kernel k;

    Vector p(source.size());
    for (int i = 0; i < source.size(); ++i) {
        p(i) = weight[i];
    }

    k.initialize(3, source, source, p, source.size(), source.size(), 40, 10);

    k.eval = [&](point &a, point &b) {
        if (a == b) return 0.;
        else {
            scalar_t r = norm(a, b);
            return (a.x * (a.x - b.x) + a.y * (a.y - b.y) + a.z * (a.z - b.z)) * 1.0 / 4 / M_PI / r / r / r;
        }

    };
    Vector q(source.size());
    k.run(q);

    std::cout << q << std::endl;


}