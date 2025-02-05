#include "attribute.h"

#include <assert.h>

#include "common.h"

bool get_xn(attribute x, int n) {
    // printf("x, n: %d,%d\n", x, n);
    assert(n > 0);
        // printf("%d\n", PARAMS.K);
    assert(n <= PARAMS.Att_num);

    // printf("x: %lld get x_n· n:%d: %d\n", x, n, (x >> (n - 1)) & 1);
    return (x >> (n - 1)) & 1ULL;
}
