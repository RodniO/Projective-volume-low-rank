#include <stdexcept>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

import pvlr;

namespace nb = nanobind;

namespace {

NB_MODULE(extension, m) {
    using pvlr::IntVector;
    nb::class_<IntVector>(m, "IntVector")
        .def("from_buffer",
             [](nb::ndarray<int, nb::ro> const &arr) {
                 if (arr.ndim() != 1) {
                     throw std::runtime_error("Wrong number of dimensions.");
                 }
                 auto size = static_cast<int>(arr.shape(0));
                 return IntVector::FromBuffer(arr.data(), size);
             })
        .def("to_buffer", [](IntVector const &vec) {
            auto buffer = vec.ToBuffer();
            auto owner = new std::vector<int>(std::move(buffer));
            nb::capsule deleter(owner, [](void *ptr) noexcept {
                delete reinterpret_cast<std::vector<int> *>(ptr);
            });
            return nb::ndarray<int, nb::numpy, nb::c_contig>(
                owner->data(), {owner->size()}, deleter);
        });

    using pvlr::Matrix;
    nb::class_<Matrix>(m, "Matrix")
        .def(nb::init<>())
        .def("from_buffer",
             [](nb::ndarray<double, nb::ro, nb::c_contig> const &arr) {
                 if (arr.ndim() != 2) {
                     throw std::runtime_error("Wrong number of dimensions.");
                 }
                 auto rows = static_cast<int>(arr.shape(0));
                 auto cols = static_cast<int>(arr.shape(1));
                 return Matrix::FromBuffer(arr.data(), rows, cols);
             })
        .def("to_buffer", [](Matrix const &mat) {
            auto shape = mat.GetShape();
            size_t rows = static_cast<size_t>(shape.first);
            size_t cols = static_cast<size_t>(shape.second);

            auto buffer = mat.ToBuffer();
            auto owner = new std::vector<double>(std::move(buffer));
            nb::capsule deleter(owner, [](void *ptr) noexcept {
                delete reinterpret_cast<std::vector<double> *>(ptr);
            });

            return nb::ndarray<double, nb::numpy, nb::c_contig>(
                owner->data(), {rows, cols}, deleter);
        });

    m.def("aca", &pvlr::ACA);
    m.def("maxvol", &pvlr::Maxvol);
    m.def("maxvol_proj", &pvlr::MaxvolProj);
}

} // namespace
