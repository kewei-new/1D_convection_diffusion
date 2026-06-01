#include "mfem.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace mfem;
using namespace std;

namespace
{
constexpr real_t pi = 3.141592653589793238462643383279502884;

real_t ExactN(const Vector &x)
{
   return 1.0 + 0.1*std::sin(pi*x[0])*std::sin(pi*x[1]);
}

real_t ExactPhi(const Vector &x)
{
   return std::sin(2.0*pi*x[0])*std::sin(pi*x[1]);
}

void WriteMetrics(int elements, int order, int dofs, real_t h, real_t n_err,
                  real_t phi_err, real_t charge_proxy)
{
   ofstream out("metrics.csv");
   out << setprecision(16);
   out << "case_name,elements,dimension,order,dofs,h_max,n_l2_error,phi_l2_error,"
          "n_relative_l2_error,phi_relative_l2_error,charge_proxy,status\n";
   out << "dd2d_smooth_mms," << elements << ",2," << order << "," << dofs << ","
       << h << "," << n_err << "," << phi_err << "," << n_err << ","
       << phi_err << "," << charge_proxy << ",implemented\n";
}
}

int main(int argc, char *argv[])
{
   int elements = 8;
   int order = 2;
   int precision = 8;

   OptionsParser args(argc, argv);
   args.AddOption(&elements, "-n", "--elements", "Number of elements per direction.");
   args.AddOption(&order, "-o", "--order", "H1 finite element order.");
   args.AddOption(&precision, "-p", "--precision", "Output precision.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   cout.precision(precision);

   Mesh mesh = Mesh::MakeCartesian2D(elements, elements, Element::QUADRILATERAL, true, 1.0, 1.0);
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
   WriteMetrics(elements, order, fes.GetTrueVSize(), h, n_err, phi_err, charge_proxy);

   cout << "case=dd2d_smooth_mms"
        << " dofs=" << fes.GetTrueVSize()
        << " n_l2_error=" << n_err
        << " phi_l2_error=" << phi_err << endl;
   return 0;
}
