#include "miniroot.hpp"
#include "typeutils.hpp"
#include "nanobind/nanobind.h"
#include "nanobind/stl/string.h"
#include "nanobind/stl/vector.h"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(pyminiroot, m)
{
    nb::class_<Miniroot>(m, "Miniroot")
        .def(nb::init<const std::string&>())
        .def("branches", &Miniroot::GetBranches)
        .def("get", &Miniroot::Get, "branch_name"_a);
}
