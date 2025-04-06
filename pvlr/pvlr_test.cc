#include <cassert>
#include <random>
#include <vector>

import pvlr;

void TestMatrix(void) {
    Matrix mat = {nullptr};
    pvlr_matrix_create(&mat, 10, 10);
    pvlr_matrix_destroy(&mat);
}

void TestIntVector(void) {
    std::vector<int> ix_arange(10);
    for (auto it = 0u; it != ix_arange.size(); ++it) {
        ix_arange[it] = static_cast<int>(it + 1);
    }
    IntVector ix = {nullptr};
    pvlr_int_vector_create(&ix, 10);
    pvlr_int_vector_from_buffer(&ix, ix_arange.data());
    pvlr_int_vector_destroy(&ix);
}

void TestMaxvol(void) {
    constexpr int n = 1'000;
    constexpr int k = 10;
    constexpr int maxsteps = 2;
    constexpr int maxswaps = 4 * k;

    std::mt19937_64 rng(42);
    std::uniform_real_distribution<> uniform(-1, 1);
    std::vector<double> elems;
    elems.reserve(n * n);
    for (auto it = 0; it != n * n; ++it) {
        elems.push_back(uniform(rng));
    }
    auto A = pvlr::Matrix::FromBuffer(elems, n, n);
    {
        auto [rows, cols] = A.GetShape();
        assert(rows == n && "Wrong number of rows.");
        assert(cols == n && "Wrong number of columns.");
    }

    pvlr::Matrix C, CA, AR;
    auto per1 = pvlr::IntVector::Iota(1, n + 1);
    auto per2 = pvlr::IntVector::Iota(1, n + 1);

    pvlr_maxvol(A.get(), n, n, k, per1.get(), per2.get(), C.get(), AR.get(),
                maxsteps, maxswaps, CA.get(), true);

    {
        auto [rows, cols] = C.GetShape();
        assert(rows == n && "Wrong number of rows.");
        assert(cols == n && "Wrong number of columns.");
    }
    {
        auto [rows, cols] = AR.GetShape();
        assert(rows == k && "Wrong number of rows.");
        assert(cols == k && "Wrong number of columns.");
    }
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
    TestIntVector();
    TestMatrix();
    TestMaxvol();
    return 0;
}
