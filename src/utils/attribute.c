#include "attribute.h"

#include <assert.h>

#include "common.h"

bool get_xn(attribute x, int n) {
    // printf("x, n: %d,%d\n", x, n);
    assert(n > 0);
    assert(n <= PARAMS.K);
    return (x >> (n - 1)) & 1;
}
