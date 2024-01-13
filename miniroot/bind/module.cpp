#include "miniroot.hpp"
#include "nanobind/nanobind.h"
#include "nanobind/operators.h"
#include "nanobind/stl/string.h"
#include "nanobind/stl/vector.h"
#include "typeutils.hpp"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(pyminiroot, m)
{
    nb::enum_<Type>(m, "DataType")
        .value("Float", Type::Float)
        .value("Double", Type::Double)
        .value("Byte", Type::Byte)
        .value("Short", Type::Short)
        .value("Integer", Type::Integer)
        .value("Long", Type::Long);

    nb::class_<Miniroot>(m, "MinirootBase")
        .def(nb::init<const std::string&>())
        .def("branches", &Miniroot::GetBranches)
        .def("get_float", &Miniroot::Get<float>, "branch_name"_a)
        .def("get_double", &Miniroot::Get<float>, "branch_name"_a)
        .def("get_byte", &Miniroot::GetBytes, "branch_name"_a)
        .def("get_short", &Miniroot::Get<float>, "branch_name"_a)
        .def("get_integer", &Miniroot::Get<float>, "branch_name"_a)
        .def("get_long", &Miniroot::Get<float>, "branch_name"_a);
}
