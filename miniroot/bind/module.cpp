#include "Miniroot.hpp"
#include "nanobind/nanobind.h"
#include "nanobind/stl/string.h"
#include "nanobind/stl/vector.h"

namespace nb = nanobind;
using namespace nb::literals;

/**
 * @brief Construct the pyminiroot python module 
 */

NB_MODULE(pyminiroot, m)
{
    nb::class_<Miniroot>(m, "MinirootBase")
        .def(nb::init<const std::string&>())
        .def("branches", &Miniroot::GetBranches)
        .def("get_float", &Miniroot::Get<float>, "branch_name"_a)
        .def("get_double", &Miniroot::Get<double>, "branch_name"_a)
        .def("get_byte", &Miniroot::GetBytes, "branch_name"_a)
        .def("get_short", &Miniroot::Get<short>, "branch_name"_a)
        .def("get_integer", &Miniroot::Get<int>, "branch_name"_a)
        .def("get_long", &Miniroot::Get<long>, "branch_name"_a);
}
