/*!
        @file    gaugeConfig.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "gaugeConfig.h"
#include "Tools/randomNumberManager.h"
#include "io_format_gauge.h"

const std::string GaugeConfig::class_name = "GaugeConfig";

//====================================================================
class GaugeConfig::DataSource
{
 public:
  virtual ~DataSource() {}
  virtual void set(Field_G& U) = 0;
};

//--------------------------------------------------------------------
class GaugeConfig::DataSource_Unit : public DataSource
{
 public:
  virtual void set(Field_G& U)
  {
    U.set_unit();
  }
};

//--------------------------------------------------------------------
class GaugeConfig::DataSource_Random : public DataSource
{
 public:
  DataSource_Random()
  {
    m_rand = RandomNumberManager::getInstance();
  }

  ~DataSource_Random()
  {
  }

  virtual void set(Field_G& U)
  {
    U.set_random(m_rand);
  }

 private:
  RandomNumbers *m_rand;
};

//====================================================================
GaugeConfig::GaugeConfig(const string& type)
  : m_vl(CommonParameters::Vlevel())
  , m_fieldio(NULL)
  , m_format(NULL)
  , m_datasource(NULL)
{
  if (type == "Text") {
    m_format  = new IO_Format::Gauge::ILDG_Format;
    m_fieldio = new FieldIO_Text(m_format);
  } else if (type == "Text_4x4x4x8") {
    m_format  = new IO_Format::Gauge::ILDG_Format;
    m_fieldio = new FieldIO_Text_4x4x4x8(m_format);
  } else if (type == "Binary") {
    m_format  = new IO_Format::Gauge::ILDG_Format;
    m_fieldio = new FieldIO_Binary(m_format);
  } else if (type == "Fortran_JLQCD") {
    m_format  = new IO_Format::Gauge::JLQCD_Format;
    m_fieldio = new FieldIO_Fortran(m_format);
  } else if (type == "Fortran_ILDG") {
    m_format  = new IO_Format::Gauge::ILDG_Format;
    m_fieldio = new FieldIO_Fortran(m_format);
  } else if (type == "ILDG") {
    m_format  = new IO_Format::Gauge::ILDG_Format;
    m_fieldio = new FieldIO_LIME(m_format);
  } else if (type == "Binary_Parallel") {
    m_format  = new IO_Format::Gauge::ILDG_Format;
    m_fieldio = new FieldIO_Binary_Parallel(m_format);
  } else if (type == "Binary_Distributed") {
    m_format  = new IO_Format::Gauge::ILDG_Format;
    m_fieldio = new FieldIO_Binary_Distributed(m_format);
  } else if (type == "ILDG_Parallel") {
    m_format  = new IO_Format::Gauge::ILDG_Format;
    m_fieldio = new FieldIO_LIME_Parallel(m_format);
  } else if (type == "NERSC") {
    m_format  = new IO_Format::Gauge::ILDG_Format;
    m_fieldio = new FieldIO_NERSC(m_format);
  } else if ((type == "None") || (type == "Null")) {  // for backward compatibility
    m_format  = new IO_Format::Trivial_Format;
    m_fieldio = new FieldIO_None(m_format);
  } else if (type == "Unit") {
    m_datasource = new DataSource_Unit;
  } else if (type == "Random") {
    m_datasource = new DataSource_Random;
  } else {
    vout.crucial("Error at %s: unsupported type \"%s\".\n", class_name.c_str(), type.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
GaugeConfig::~GaugeConfig()
{
  if (m_fieldio) delete m_fieldio;
  if (m_format) delete m_format;
  if (m_datasource) delete m_datasource;
}


//====================================================================
void GaugeConfig::read(Field_G& U, const string& filename)
{
  if (m_fieldio) {
    return m_fieldio->read_file(U, filename);
  }

  if (m_datasource) {
    return m_datasource->set(U);
  }

  vout.crucial("Error at %s::read_file(): FieldIO or DataSource not set.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
void GaugeConfig::write(Field_G& U, const string& filename)
{
  if (m_fieldio) {
    if (filename != "NO_OUTPUT") {
      return m_fieldio->write_file(U, filename);
    } else {
      return;  // do nothig
    }
  }

  vout.crucial("Error at %s::write_file(): FieldIO not set.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
//============================================================END=====
