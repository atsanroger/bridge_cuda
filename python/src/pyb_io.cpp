/*!
        @file    pyb_io.cpp

        @brief   python binding of file I/O classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#include "pyb_defs.h"
#include "pyb_field.h"
#include "pyb_io.h"
#include "bridge.h"
#include "Field/field.h"
#include "IO/gaugeConfig.h"
#include "IO/dataIO.h"
#include "IO/dataIO_Text.h"
#include "IO/fieldIO.h"
#include "IO/io_format.h"
#include "Parameters/parameters.h"

enum IO_Format_type {
  UNKNOWN,
  TRIVIAL,
  GAUGE_ILDG,
  GAUGE_JLQCD
};

FieldIO *create_FieldIO(const std::string& subtype, IO_Format_type format_type)
{
  const IO_Format::Format *format = nullptr;

  switch (format_type)
  {
  case TRIVIAL: format = IO_Format::Trivial; break;
  case GAUGE_ILDG: format = IO_Format::Gauge::ILDG; break;
  case GAUGE_JLQCD: format = IO_Format::Gauge::JLQCD; break;
  default:
    warn("unknown io_format_type %d", format_type);
  }

  if (!format) return nullptr;


  FieldIO *obj = nullptr;

  if (subtype == "None") {
    return new FieldIO_None(format);

  } else if (subtype == "Text") {
    return new FieldIO_Text(format);

  } else if (subtype == "Binary") {
    return new FieldIO_Binary(format);

  } else if (subtype == "Fortran") {
    return new FieldIO_Fortran(format);

#ifdef USE_LIME
  } else if (subtype == "LIME") {
    return new FieldIO_LIME(format);
#endif

#ifdef USE_MPI
  } else if (subtype == "Binary_Parallel") {
    return new FieldIO_Binary_Parallel(format);

  } else if (subtype == "Binary_Distributed") {
    return new FieldIO_Binary_Distributed(format);

  } else if (subtype == "LIME_Parallel") {
    return new FieldIO_LIME_Parallel(format);
#endif

  } else {
    warn("unknown fieldio subtype \"%s\"", subtype.c_str());
  }

  return obj;
}

void pybind_io(py::module m)
{
  typedef void (GaugeConfig::*method_signature)(Field_G*, const std::string&);

  //! GaugeConfig
  //
  py::class_<GaugeConfig>(m, "GaugeConfig", R"pybdoc(
File I/O of gauge configurations
)pybdoc")
    .def(py::init<std::string>(), REDIRECT)

    // .def("read",  (method_signature) &GaugeConfig::read)
    // .def("write", (method_signature) &GaugeConfig::write)

    .def("read_file",  (method_signature) &GaugeConfig::read_file, REDIRECT)
    .def("write_file", (method_signature) &GaugeConfig::write_file, REDIRECT)
    ;

  //! DataIO
  //
  // example
  //   io = DataIO(type, [options])
  //
  // arguments
  //   type: "Text", "Binary" (todo)
  //   options:
  //     format_precision (Text) -- number of digits (default=14)
  //     endian (Binary) -- "native"(default), "big_endian", "little_endian"  
  //
  py::class_<DataIO>(m, "DataIO", R"pybdoc(
File I/O of general collection of data
)pybdoc")
    .def(py::init([](const std::string& type,
                     int format_prec,
                     const std::string& endian
                  ) -> DataIO*
                  {
                    if (type == "Text") {
                      auto obj = new DataIO_Text;
                      obj->set_parameter(format_prec);
                      return obj;
                    }
                    return nullptr;
                  }),
         py::arg("type"),
         py::arg("format_precision") = 14, //default_format_precision,
         py::arg("endian") = "native",
         REDIRECT
    )

//    .def("read_file", [](DataIO& self, std::vector<double>& data, const std::string& filename){
//        return self.read_file(data, filename);
//      }, REDIRECT)

    .def("read_file", [](DataIO& self, const std::string& filename, int length){
        std::vector<double> buf(length);
        self.read_file(buf, filename);
        return buf;
      }, py::arg("file"), py::arg("length"),
      REDIRECT)

    .def("write_file", [](DataIO& self, std::vector<double>& data, const std::string& filename, bool append){
        return self.write_file(data, filename, append);
      }, py::arg("data"), py::arg("file"), py::arg("append")=true,
      REDIRECT)
    ;
          
  //! FieldIO
  //
  py::enum_<IO_Format_type>(m, "IO_Format_type")
    .value("IO_Format_Unknown", IO_Format_type::UNKNOWN)
    .value("IO_Format_Trivial", IO_Format_type::TRIVIAL)
    .value("IO_Format_Gauge_ILDG", IO_Format_type::GAUGE_ILDG)
    .value("IO_Format_Gauge_JLQCD", IO_Format_type::GAUGE_JLQCD)
    .export_values();


  // typedef void (FieldIO::*func_single_t)(Field*, const std::string&);
  // typedef void (FieldIO::*func_array_t)(std::vector<Field*>&, const std::string&);

  py::class_<FieldIO>(m, "FieldIO", R"pybdoc(
File I/O of space-time distributed data
)pybdoc")
    .def(py::init([](const std::string& type, IO_Format_type format) -> FieldIO*
                  {
                    return create_FieldIO(type, format);
                  }),
         REDIRECT)

    .def("read_file", [](FieldIO& self, Field* f, const std::string& filename)
         {
           self.read_file(f, filename);
         }, REDIRECT)
    .def("write_file", [](FieldIO& self, Field* f, const std::string& filename)
         {
           self.write_file(f, filename);
         }, REDIRECT)
    
    .def("read_file", [](FieldIO& self, std::vector<Field*>& ff, const std::string& filename)
         {
           self.read_file(ff, filename);
         }, REDIRECT)
    .def("write_file", [](FieldIO& self, std::vector<Field*>& ff, const std::string& filename)
         {
           self.write_file(ff, filename);
         }, REDIRECT)

    // .def("read_file",  (func_single_t) &FieldIO::read_file,  REDIRECT)
    // .def("write_file", (func_single_t) &FieldIO::write_file, REDIRECT)

    // .def("read_file",  (func_array_t)  &FieldIO::read_file,  REDIRECT)
    // .def("write_file", (func_array_t)  &FieldIO::write_file, REDIRECT)

    .def("deliver", &FieldIO::deliver, REDIRECT)
    .def("gather",  &FieldIO::gather,  REDIRECT)
    ;
  

}
