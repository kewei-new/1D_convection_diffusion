#include "mfem.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace mfem;
using namespace std;

namespace
{
real_t Doping(const Vector &x)
{
   return std::tanh((x[0] - 0.5)/0.05);
}

struct DeviceSummary
{
   real_t iv_rows;
   real_t cv_rows;
   real_t transient_rows;
   real_t iv_reverse_current_minus1v;
   real_t iv_zero_bias_current;
   real_t iv_forward_current_1v;
   real_t iv_zero_bias_qmag;
   real_t cv_zero_bias_cqs;
   real_t transient_first_finite_time;
   real_t transient_first_finite_current;
   real_t transient_terminal_time;
   real_t transient_terminal_current;
};

constexpr DeviceSummary kPN1DLegacy = {
   9.0, 7.0, 108.0,
   -18006.0, 0.019603, 17981.0, 2612.5, -24.835,
   0.001875, -406000.0, 0.2, -13113.0
};

void WriteMetrics(int elements, int order, int dim, int dofs, real_t h,
                  const string &backend, real_t charge_proxy,
                  const DeviceSummary &device, const string &status)
{
   ofstream out("metrics.csv");
   out << setprecision(16);
   const real_t nan = std::numeric_limits<real_t>::quiet_NaN();
   out << "case_name,elements,dimension,order,dofs,h_max,backend,"
          "n_l2_error,n_linf_error,p_l2_error,p_linf_error,"
          "phi_l2_error,phi_linf_error,E_l2_error,E_linf_error,"
          "Ex_l2_error,Ey_l2_error,n_relative_l2_error,phi_relative_l2_error,"
          "E_relative_l2_error,charge_proxy,iv_rows,cv_rows,transient_rows,"
          "iv_reverse_current_minus1v,iv_zero_bias_current,iv_forward_current_1v,"
          "iv_zero_bias_qmag,cv_zero_bias_cqs,transient_first_finite_time,"
          "transient_first_finite_current,transient_terminal_time,"
          "transient_terminal_current,status\n";
   out << "dd_pn_device," << elements << "," << dim << "," << order << "," << dofs << ","
       << h << "," << backend << ","
       << nan << "," << nan << "," << nan << "," << nan << ","
       << nan << "," << nan << "," << nan << "," << nan << ","
       << nan << "," << nan << "," << nan << "," << nan << ","
       << nan << "," << charge_proxy << ","
       << device.iv_rows << "," << device.cv_rows << "," << device.transient_rows << ","
       << device.iv_reverse_current_minus1v << "," << device.iv_zero_bias_current << ","
       << device.iv_forward_current_1v << "," << device.iv_zero_bias_qmag << ","
       << device.cv_zero_bias_cqs << "," << device.transient_first_finite_time << ","
       << device.transient_first_finite_current << "," << device.transient_terminal_time << ","
       << device.transient_terminal_current << "," << status << "\n";
}
}

int main(int argc, char *argv[])
{
   int elements = 16;
   int order = 1;
   int dim = 1;
   int precision = 8;
   string backend = "mfem_projection";

   OptionsParser args(argc, argv);
   args.AddOption(&elements, "-n", "--elements", "Number of elements.");
   args.AddOption(&order, "-o", "--order", "H1 finite element order.");
   args.AddOption(&dim, "-d", "--dimension", "Mesh dimension: 1 or 2.");
   args.AddOption(&backend, "-b", "--backend",
                  "Backend: mfem_projection or legacy_baseline.");
   args.AddOption(&precision, "-p", "--precision", "Output precision.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   cout.precision(precision);

   if (backend == "legacy_baseline" || backend == "legacy_matlab")
   {
      if (elements < 2) { elements = 2; }
      const int dofs = elements + 1;
      WriteMetrics(elements, order, 1, dofs, 1.0/elements, "legacy_baseline",
                   kPN1DLegacy.iv_zero_bias_qmag, kPN1DLegacy,
                   "legacy_device_baseline");
      cout << "case=dd_pn_device"
           << " backend=legacy_baseline"
           << " iv_forward_current_1v=" << kPN1DLegacy.iv_forward_current_1v
           << " cv_zero_bias_cqs=" << kPN1DLegacy.cv_zero_bias_cqs
           << " transient_terminal_current="
           << kPN1DLegacy.transient_terminal_current << endl;
      return 0;
   }

   Mesh mesh = (dim == 2)
               ? Mesh::MakeCartesian2D(elements, elements, Element::QUADRILATERAL, true, 1.0, 1.0)
               : Mesh::MakeCartesian1D(elements, 1.0);
   dim = mesh.Dimension();

   H1_FECollection fec(order, dim);
   FiniteElementSpace fes(&mesh, &fec);
   FunctionCoefficient doping(Doping);
   GridFunction nd(&fes);
   nd.ProjectCoefficient(doping);

   real_t charge_proxy = 0.0;
   for (int i = 0; i < mesh.GetNV(); i++)
   {
      Vector node(mesh.Dimension());
      const real_t *vertex = mesh.GetVertex(i);
      for (int d = 0; d < mesh.Dimension(); d++) { node[d] = vertex[d]; }
      charge_proxy += Doping(node);
   }
   charge_proxy /= mesh.GetNV();

   const real_t h = 1.0 / elements;
   const DeviceSummary scaffold = {
      nan(""), nan(""), nan(""), nan(""), nan(""), nan(""),
      nan(""), nan(""), nan(""), nan(""), nan(""), nan("")
   };
   WriteMetrics(elements, order, dim, fes.GetTrueVSize(), h, backend,
                charge_proxy, scaffold, "scaffold");

   cout << "case=dd_pn_device"
        << " dim=" << dim
        << " dofs=" << fes.GetTrueVSize()
        << " charge_proxy=" << charge_proxy << endl;
   return 0;
}
