#ifndef __Solvers_options_hpp
#define __Solvers_options_hpp

#include  <iostream>
#include  <Alat/map.hpp>
#include  <Mesh/enums.hpp>

/*--------------------------------------------------------------------------*/
namespace solver_options
{
  /*--------------------------------------------------------------------------*/
  enum output_manager_data
  {
    datadir,
    datafilename,
    solutionfilename,
    meshfilename
  };
  struct output_manager
  {
    std::string path;
    alat::Map<output_manager_data, std::string> data;
    inline explicit output_manager(std::string datadir_in="Data")
    {
      data[datadir]=datadir_in;
      data[solutionfilename]="solution";
      data[datafilename]="data";
      data[meshfilename]="mesh";
    }
    void init();
    std::string getSolutionFileName(meshEnums::meshunittype type=meshEnums::None, std::string name="u", int it=-1) const;
    std::string getDataFileName(meshEnums::meshunittype type=meshEnums::None, int it=-1) const;
    std::string getMeshFileName() const;
  };
  std::string output_manager_data_to_string(output_manager_data data);
  std::ostream& operator<<(std::ostream& os, const output_manager& om);

  /*--------------------------------------------------------------------------*/
  typedef unsigned int flag_type;
  struct opts
  {
    flag_type flags;
    inline explicit opts(const flag_type in_flags) : flags(in_flags) {}
    inline const opts operator+(const opts& rhs) const
    {
      const opts result( flags | rhs.flags );
      return result;
    }
    inline opts& operator+=(const opts& rhs)
    {
      flags = ( flags | rhs.flags );
      return *this;
    }
    inline bool operator<(const opts& rhs) const
    {
      return this->flags<rhs.flags;
    }
  };
  std::ostream& operator<<(std::ostream& os, const opts& op);
  bool has_armamat(const opts& op);
  bool has_dynamic(const opts& op);
  bool has_erasedatadir(const opts& op);

  static const flag_type flag_none    = flag_type(0      );
  static const flag_type flag_armamat = flag_type(1u << 0);
  static const flag_type flag_dynamic = flag_type(1u << 1);
  static const flag_type flag_erasedatadir = flag_type(1u << 2);

  struct opts_none    : public opts { inline opts_none()    : opts(flag_none     ) {} };
  struct opts_armamat : public opts { inline opts_armamat() : opts(flag_armamat  ) {} };
  struct opts_dynamic  : public opts { inline opts_dynamic()  : opts(flag_dynamic) {} };
  struct opts_erasedatadir  : public opts { inline opts_erasedatadir()  : opts(flag_erasedatadir) {} };

  static const opts_none      none;
  static const opts_armamat   armamat;
  static const opts_dynamic   dynamic;
  static const opts_erasedatadir   erasedatadir;

  namespace errors
  {
    static const flag_type flag_none = flag_type(0      );
    static const flag_type flag_L2   = flag_type(1u << 0);
    static const flag_type flag_H1   = flag_type(1u << 1);
    static const flag_type flag_L1   = flag_type(1u << 2);
    static const flag_type flag_Linf   = flag_type(1u << 3);
    static const flag_type flag_E   = flag_type(1u << 4);

    struct opts : solver_options::opts
    {
      inline explicit opts(const flag_type in_flags) : solver_options::opts(in_flags) {}
    };
    std::ostream& operator<<(std::ostream& os, const opts& op);

    struct opts_none : public opts { inline opts_none() : opts(flag_none) {} };
    struct opts_L2   : public opts { inline opts_L2()   : opts(flag_L2  ) {} };
    struct opts_H1   : public opts { inline opts_H1()   : opts(flag_H1  ) {} };
    struct opts_L1   : public opts { inline opts_L1()   : opts(flag_L1  ) {} };
    struct opts_Linf   : public opts { inline opts_Linf()   : opts(flag_Linf  ) {} };
    struct opts_E   : public opts { inline opts_E()   : opts(flag_E  ) {} };

    bool has_L2(const opts& op);
    bool has_H1(const opts& op);
    bool has_L1(const opts& op);
    bool has_Linf(const opts& op);
    bool has_E(const opts& op);

    static const opts_none none;
    static const opts_L2   L2;
    static const opts_H1   H1;
    static const opts_L1   L1;
    static const opts_Linf Linf;
    static const opts_E E;
  }
  namespace pdepart
  {
    static const flag_type flag_none    = flag_type(0      );
    static const flag_type flag_cell    = flag_type(1u << 0);
    static const flag_type flag_bdry    = flag_type(1u << 1);
    static const flag_type flag_iside    = flag_type(1u << 2);

    struct opts : solver_options::opts
    {
      inline explicit opts(const flag_type in_flags) : solver_options::opts(in_flags) {}
      inline const opts operator+(const opts& rhs) const
      {
        const opts result( flags | rhs.flags );
        return result;
      }
    };
    std::ostream& operator<<(std::ostream& os, const opts& op);

    struct opts_none    : public opts { inline opts_none()    : opts(flag_none    ) {} };
    struct opts_cell    : public opts { inline opts_cell()    : opts(flag_cell    ) {} };
    struct opts_bdry    : public opts { inline opts_bdry()    : opts(flag_bdry    ) {} };
    struct opts_iside    : public opts { inline opts_iside()    : opts(flag_iside    ) {} };

    bool has_cell(const opts& op);
    bool has_bdry(const opts& op);
    bool has_iside(const opts& op);

    static const opts_none none;
    static const opts_cell     cell;
    static const opts_bdry     bdry;
    static const opts_iside    iside;
  }
}

/*--------------------------------------------------------------------------*/
#endif
