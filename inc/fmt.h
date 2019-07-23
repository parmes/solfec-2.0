/* https://stackoverflow.com/questions/11989374/floating-point-format-for-stdostream */
#include <iostream>
#include <iomanip>
#include <stdio.h>

class FMT
{
public:
  explicit FMT(const char* fmt): m_fmt(fmt) {}

private:
  class fmter /* actual worker class */
  {
  public:
    explicit fmter(std::ostream& strm, const FMT& fmt): m_strm(strm), m_fmt(fmt.m_fmt) {}

    /* output next object (any type) to stream */
    template<typename TYPE>
    std::ostream& operator<<(const TYPE& value)
    {
      /* return m_strm << "FMT(" << m_fmt << "," << value << ")"; */
      char buf[128]; /* enlarge as needed */
      snprintf(buf, sizeof(buf), m_fmt, value);
      return m_strm << buf;
    }
  private:
    std::ostream& m_strm;
    const char* m_fmt;
  };

  const char* m_fmt; /* save fmt string for inner class */

  /* return derived stream to allow operator overloading */
  friend FMT::fmter operator<<(std::ostream& strm, const FMT& fmt)
  {
    return FMT::fmter(strm, fmt);
  }
};
