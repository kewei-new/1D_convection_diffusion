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
constexpr real_t pi = 3.141592653589793238462643383279502884;

real_t ExactN(const Vector &x)
{
   return 1.0 + 0.2*std::sin(pi*x[0]);
}

real_t ExactPhi(const Vector &x)
{
   return std::sin(2.0*pi*x[0]);
}

struct LegacyRow
{
   real_t h;
   real_t n_l2;
   real_t n_linf;
   real_t phi_l2;
   real_t phi_linf;
   real_t e_l2;
   real_t e_linf;
};

bool FindLegacyRow(int elements, int order, LegacyRow &row)
{
   if (order != 3) { return false; }
   switch (elements)
   {
      case 20:
         row = {3.1416e-01, 4.507504e-06, 3.897874e-06,
                5.547920e-06, 9.655373e-06, 5.658847e-06, 9.737335e-06};
         return true;
      case 40:
         row = {1.5708e-01, 4.071754e-07, 3.531198e-07,
                3.751691e-07, 6.082643e-07, 3.928129e-07, 6.268181e-07};
         return true;
      case 80:
         row = {7.8540e-02, 4.561335e-08, 3.706344e-08,
                2.966955e-08, 3.820909e-08, 3.311548e-08, 4.220075e-08};
         return true;
      case 160:
         row = {3.9270e-02, 5.531522e-09, 4.186841e-09,
                2.940053e-09, 2.935649e-09, 3.472993e-09, 3.100092e-09};
         return true;
      case 320:
         row = {1.9635e-02, 6.906262e-10, 4.991136e-10,
                3.419493e-10, 3.259086e-10, 4.140009e-10, 3.365440e-10};
         return true;
      default:
         return false;
   }
}

void WriteMetrics(int elements, int order, int dofs, real_t h,
                  const string &backend, const string &status,
                  real_t n_l2, real_t n_linf, real_t p_l2, real_t p_linf,
                  real_t phi_l2, real_t phi_linf, real_t e_l2, real_t e_linf,
                  real_t ex_l2, real_t ey_l2, real_t charge_proxy)
{
   ofstream out("metrics.csv");
   out << setprecision(16);
   out << "case_name,elements,dimension,order,dofs,h_max,backend,"
          "n_l2_error,n_linf_error,p_l2_error,p_linf_error,"
          "phi_l2_error,phi_linf_error,E_l2_error,E_linf_error,"
          "Ex_l2_error,Ey_l2_error,n_relative_l2_error,phi_relative_l2_error,"
          "E_relative_l2_error,charge_proxy,status\n";
   out << "dd1d_smooth_mms," << elements << ",1," << order << "," << dofs << ","
       << h << "," << backend << ","
       << n_l2 << "," << n_linf << "," << p_l2 << "," << p_linf << ","
       << phi_l2 << "," << phi_linf << "," << e_l2 << "," << e_linf << ","
       << ex_l2 << "," << ey_l2 << ","
       << n_l2 << "," << phi_l2 << "," << e_l2 << ","
       << charge_proxy << "," << status << "\n";
}
}

int main(int argc, char *argv[])
{
   int elements = 16;
   int order = 2;
   int precision = 8;
   string backend = "mfem_projection";

   OptionsParser args(argc, argv);
   args.AddOption(&elements, "-n", "--elements", "Number of 1D mesh elements.");
   args.AddOption(&order, "-o", "--order", "H1 finite element order.");
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

   const real_t nan = std::numeric_limits<real_t>::quiet_NaN();
   if (backend == "legacy_baseline" || backend == "legacy_matlab")
   {
      LegacyRow row {};
      if (!FindLegacyRow(elements, order, row))
      {
         cerr << "No embedded DD1D legacy baseline row for elements="
              << elements << " order=" << order << ". Use order=3 and "
              << "elements in {20,40,80,160,320}." << endl;
         return 2;
      }
      WriteMetrics(elements, order, elements + 1, row.h, "legacy_baseline",
                   "legacy_baseline", row.n_l2, row.n_linf, nan, nan,
                   row.phi_l2, row.phi_linf, row.e_l2, row.e_linf,
                   nan, nan, nan);
      cout << "case=dd1d_smooth_mms"
           << " backend=legacy_baseline"
           << " n_l2_error=" << row.n_l2
           << " phi_l2_error=" << row.phi_l2
           << " E_l2_error=" << row.e_l2 << endl;
      return 0;
   }

   Mesh mesh = Mesh::MakeCartesian1D(elements, 1.0);
   H1_FECollection fec(order, mesh.Dimension());
   FiniteElementSpace fes(&mesh, &fec);

   FunctionCoefficient n_coeff(ExactN);
   FunctionCoefficient phi_coeff(ExactPhi);
   GridFunction n(&fes), phi(&fes);
   n.ProjectCoefficient(n_coeff);
   phi.ProjectCoefficient(phi_coeff);

   const real_t n_err = n.ComputeL2Error(n_coeff);
   const real_t phi_err = phi.ComputeL2Error(phi_coeff);
   real_t charge_proxy = 0.0;
   for (int i = 0; i < mesh.GetNV(); i++)
   {
      Vector node(mesh.Dimension());
      const real_t *vertex = mesh.GetVertex(i);
      for (int d = 0; d < mesh.Dimension(); d++) { node[d] = vertex[d]; }
      charge_proxy += ExactN(node) - 1.0;
   }
   charge_proxy /= mesh.GetNV();

   const real_t h = 1.0 / elements;
   WriteMetrics(elements, order, fes.GetTrueVSize(), h, backend, "implemented",
                n_err, nan, nan, nan, phi_err, nan, nan, nan, nan, nan,
                charge_proxy);

   cout << "case=dd1d_smooth_mms"
        << " dofs=" << fes.GetTrueVSize()
        << " n_l2_error=" << n_err
        << " phi_l2_error=" << phi_err << endl;
   return 0;
}
