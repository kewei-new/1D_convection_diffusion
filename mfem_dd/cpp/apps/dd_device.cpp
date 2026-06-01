#include "mfem.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace mfem;
using namespace std;

namespace
{
real_t Doping(const Vector &x)
{
   return std::tanh((x[0] - 0.5)/0.05);
}

void WriteMetrics(int elements, int order, int dim, int dofs, real_t h,
                  real_t charge_proxy)
{
   ofstream out("metrics.csv");
   out << setprecision(16);
   out << "case_name,elements,dimension,order,dofs,h_max,n_l2_error,phi_l2_error,"
          "n_relative_l2_error,phi_relative_l2_error,charge_proxy,status\n";
   out << "dd_pn_device," << elements << "," << dim << "," << order << "," << dofs << ","
       << h << ",0,0,0,0," << charge_proxy << ",scaffold\n";
}
}

int main(int argc, char *argv[])
{
   int elements = 16;
   int order = 1;
   int dim = 1;
   int precision = 8;

   OptionsParser args(argc, argv);
   args.AddOption(&elements, "-n", "--elements", "Number of elements.");
   args.AddOption(&order, "-o", "--order", "H1 finite element order.");
   args.AddOption(&dim, "-d", "--dimension", "Mesh dimension: 1 or 2.");
   args.AddOption(&precision, "-p", "--precision", "Output precision.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   cout.precision(precision);

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
   WriteMetrics(elements, order, dim, fes.GetTrueVSize(), h, charge_proxy);

   cout << "case=dd_pn_device"
        << " dim=" << dim
        << " dofs=" << fes.GetTrueVSize()
        << " charge_proxy=" << charge_proxy << endl;
   return 0;
}
