#include "mfem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

using namespace mfem;
using namespace std;

namespace
{
#ifndef MFEM_DD_MATLAB_DIR
#define MFEM_DD_MATLAB_DIR ""
#endif

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

constexpr real_t kPNLeft = 0.0;
constexpr real_t kPNRight = 0.6;
constexpr real_t kPNMobility = 0.75;
constexpr real_t kPNThermalVoltage = 0.138046e-4*300.0/0.1602;
constexpr real_t kPNDiffusion = kPNMobility*kPNThermalVoltage;
constexpr real_t kPNPoissonScale = 0.1602/(11.7*8.85418);
constexpr real_t kPNNi = 0.014;
constexpr real_t kPNBias = 1.5;
constexpr real_t kPNContactDensityLeft = 5.0e5;
constexpr real_t kPNContactDensityRight = 5.0e5;

struct GaussPoint
{
   real_t x;
   real_t w;
};

struct MatrixStats
{
   int rows = 0;
   int cols = 0;
   int nnz = 0;
   real_t fro_norm = 0.0;
   real_t one_norm = 0.0;
   real_t inf_norm = 0.0;
   real_t value_sum = 0.0;
   real_t abs_sum = 0.0;
   real_t weighted_sum = 0.0;
};

struct VectorStats
{
   int size = 0;
   real_t min = 0.0;
   real_t max = 0.0;
   real_t mean = 0.0;
   real_t norm2 = 0.0;
   real_t value_sum = 0.0;
   real_t abs_sum = 0.0;
   real_t weighted_sum = 0.0;
};

struct PNCurrentSummary
{
   real_t left_contact = 0.0;
   real_t right_contact = 0.0;
   real_t domain_average = 0.0;
};

int ModalIndex(int cell, int mode, int order)
{
   return cell*(order + 1) + mode;
}

real_t ModalBasis(real_t x, int mode)
{
   switch (mode)
   {
      case 0: return 1.0;
      case 1: return x;
      case 2: return x*x - 1.0/12.0;
      case 3: return x*x*x - 0.15*x;
      case 4: return (x*x - 3.0/14.0)*x*x + 3.0/560.0;
      default: throw runtime_error("Unsupported PN modal basis order.");
   }
}

real_t ModalBasisX(real_t x, int mode)
{
   switch (mode)
   {
      case 0: return 0.0;
      case 1: return 1.0;
      case 2: return 2.0*x;
      case 3: return 3.0*x*x - 0.15;
      case 4: return 4.0*x*x*x - 3.0/7.0*x;
      default: throw runtime_error("Unsupported PN modal basis derivative order.");
   }
}

vector<real_t> PNInverseMass(int order)
{
   switch (order)
   {
      case 1: return {1.0, 12.0};
      case 2: return {1.0, 12.0, 180.0};
      case 3: return {1.0, 12.0, 180.0, 2800.0};
      default: throw runtime_error("PN native snapshot supports orders 1-3.");
   }
}

vector<GaussPoint> GaussLobatto(int order)
{
   if (order == 1)
   {
      return {{-0.5, 1.0/6.0}, {0.5, 1.0/6.0}, {0.0, 2.0/3.0}};
   }
   if (order == 2)
   {
      return {{-0.5, 1.0/12.0}, {0.5, 1.0/12.0},
              {-std::sqrt(5.0)/10.0, 5.0/12.0},
              { std::sqrt(5.0)/10.0, 5.0/12.0}};
   }
   if (order == 3)
   {
      return {{-0.5, 1.0/20.0}, {0.5, 1.0/20.0},
              {-std::sqrt(21.0)/14.0, 49.0/180.0},
              { std::sqrt(21.0)/14.0, 49.0/180.0},
              {0.0, 64.0/180.0}};
   }
   throw runtime_error("Unsupported Gauss-Lobatto order.");
}

vector<GaussPoint> GaussLegendre4()
{
   return {{-0.8611363115940526/2.0, 0.3478548451374538/2.0},
           {-0.3399810435848563/2.0, 0.6521451548625461/2.0},
           { 0.3399810435848563/2.0, 0.6521451548625461/2.0},
           { 0.8611363115940526/2.0, 0.3478548451374538/2.0}};
}

real_t PNDopingProfile(real_t x)
{
   const real_t ya = 5.0e5;
   const real_t yb = 2.0e3;
   const real_t xl = 0.1;
   const real_t xr = 0.5;
   const real_t width = 0.06;
   const real_t half_width = width/2.0;
   const real_t xll = xl - half_width;
   const real_t xlr = xl + half_width;
   const real_t xrl = xr - half_width;
   const real_t xrr = xr + half_width;
   if (x < xll) { return ya; }
   if (x < xlr)
   {
      const real_t yr = (x - xll)/(width + 1.0e-20);
      const real_t one_minus = 1.0 - yr*yr*yr;
      return (ya - yb)*one_minus*one_minus*one_minus + yb;
   }
   if (x < xrl) { return yb; }
   if (x < xrr)
   {
      const real_t yr = (x - xrl)/(width + 1.0e-20);
      const real_t one_minus = 1.0 - yr*yr*yr;
      return (yb - ya)*one_minus*one_minus*one_minus + ya;
   }
   return ya;
}

real_t PNPhiLeft()
{
   return kPNThermalVoltage
          *std::log(std::max(PNDopingProfile(kPNLeft), 1.0e-12)/kPNNi);
}

real_t PNPhiRight()
{
   return PNPhiLeft() + kPNBias;
}

real_t EvalModal(const Vector &u, int cell, int order, real_t xi)
{
   real_t value = 0.0;
   for (int mode = 0; mode <= order; mode++)
   {
      value += u(ModalIndex(cell, mode, order))*ModalBasis(xi, mode);
   }
   return value;
}

void AddScaled(Vector &y, real_t a, const Vector &x)
{
   for (int i = 0; i < y.Size(); i++) { y(i) += a*x(i); }
}

Vector LinearCombination(const vector<pair<real_t, const Vector*> > &terms)
{
   Vector y(terms[0].second->Size());
   y = 0.0;
   for (const auto &term : terms) { AddScaled(y, term.first, *term.second); }
   return y;
}

Vector MatVecNew(const DenseMatrix &A, const Vector &x)
{
   Vector y(A.Height());
   A.Mult(x, y);
   return y;
}

DenseMatrix ScaleMatrix(const DenseMatrix &A, real_t scale)
{
   DenseMatrix B(A.Height(), A.Width());
   for (int i = 0; i < A.Height(); i++)
   {
      for (int j = 0; j < A.Width(); j++) { B(i, j) = scale*A(i, j); }
   }
   return B;
}

DenseMatrix TransposeScale(const DenseMatrix &A, real_t scale)
{
   DenseMatrix B(A.Width(), A.Height());
   for (int i = 0; i < A.Height(); i++)
   {
      for (int j = 0; j < A.Width(); j++) { B(j, i) = scale*A(i, j); }
   }
   return B;
}

DenseMatrix AddMatrices(const vector<pair<real_t, const DenseMatrix*> > &terms)
{
   DenseMatrix B(terms[0].second->Height(), terms[0].second->Width());
   B = 0.0;
   for (const auto &term : terms)
   {
      for (int i = 0; i < B.Height(); i++)
      {
         for (int j = 0; j < B.Width(); j++)
         {
            B(i, j) += term.first*(*term.second)(i, j);
         }
      }
   }
   return B;
}

DenseMatrix Multiply(const DenseMatrix &A, const DenseMatrix &B)
{
   DenseMatrix C(A.Height(), B.Width());
   C = 0.0;
   for (int i = 0; i < A.Height(); i++)
   {
      for (int k = 0; k < A.Width(); k++)
      {
         const real_t aik = A(i, k);
         if (aik == 0.0) { continue; }
         for (int j = 0; j < B.Width(); j++) { C(i, j) += aik*B(k, j); }
      }
   }
   return C;
}

void AddBlock(DenseMatrix &A, int row_cell, int col_cell,
              const DenseMatrix &block, int order)
{
   for (int i = 0; i <= order; i++)
   {
      for (int j = 0; j <= order; j++)
      {
         A(ModalIndex(row_cell, i, order), ModalIndex(col_cell, j, order))
            += block(i, j);
      }
   }
}

MatrixStats ComputeMatrixStats(const DenseMatrix &A)
{
   MatrixStats stats;
   stats.rows = A.Height();
   stats.cols = A.Width();
   vector<real_t> column_abs(stats.cols, 0.0), row_abs(stats.rows, 0.0);
   real_t fro2 = 0.0;
   for (int i = 0; i < stats.rows; i++)
   {
      for (int j = 0; j < stats.cols; j++)
      {
         const real_t value = A(i, j);
         if (value != 0.0) { stats.nnz++; }
         const real_t abs_value = std::abs(value);
         fro2 += value*value;
         stats.value_sum += value;
         stats.abs_sum += abs_value;
         stats.weighted_sum += (static_cast<real_t>(i + 1)
                                + 0.125*static_cast<real_t>(j + 1))*value;
         column_abs[j] += abs_value;
         row_abs[i] += abs_value;
      }
   }
   stats.fro_norm = std::sqrt(fro2);
   stats.one_norm = *std::max_element(column_abs.begin(), column_abs.end());
   stats.inf_norm = *std::max_element(row_abs.begin(), row_abs.end());
   return stats;
}

VectorStats ComputeVectorStats(const Vector &x)
{
   VectorStats stats;
   stats.size = x.Size();
   stats.min = std::numeric_limits<real_t>::infinity();
   stats.max = -std::numeric_limits<real_t>::infinity();
   real_t norm2 = 0.0;
   for (int i = 0; i < x.Size(); i++)
   {
      const real_t value = x(i);
      stats.min = std::min(stats.min, value);
      stats.max = std::max(stats.max, value);
      stats.value_sum += value;
      stats.abs_sum += std::abs(value);
      stats.weighted_sum += static_cast<real_t>(i + 1)*value;
      norm2 += value*value;
   }
   stats.mean = stats.value_sum/static_cast<real_t>(x.Size());
   stats.norm2 = std::sqrt(norm2);
   return stats;
}

DenseMatrix AssemblePNIPDGDiffusion(int elements, int order, real_t h)
{
   const int nloc = order + 1;
   const int size = elements*nloc;
   // The MATLAB PN params inherit alpha=16 from the default p=3 setup before
   // switching the PN transport degree to p=2. Keep that historical value so
   // this snapshot aligns with the legacy-derived MATLAB reference.
   const real_t alpha = 16.0;
   const real_t beta = 1.0;
   const auto gauss = GaussLobatto(order);

   DenseMatrix diffusion(size);
   diffusion = 0.0;
   DenseMatrix Mx1(nloc), Mx2LL(nloc), Mx2LR(nloc), Mx2RL(nloc), Mx2RR(nloc);
   DenseMatrix Mx3LL(nloc), Mx3LR(nloc), Mx3RR(nloc), Mx3RL(nloc);
   Mx1 = 0.0; Mx2LL = 0.0; Mx2LR = 0.0; Mx2RL = 0.0; Mx2RR = 0.0;
   Mx3LL = 0.0; Mx3LR = 0.0; Mx3RR = 0.0; Mx3RL = 0.0;

   for (int d1 = 0; d1 <= order; d1++)
   {
      for (int d2 = 0; d2 <= order; d2++)
      {
         for (const auto &g : gauss)
         {
            Mx1(d1, d2) += (1.0/h)*g.w*ModalBasisX(g.x, d2)
                           *ModalBasisX(g.x, d1);
         }
         Mx2LL(d1, d2) = (1.0/h)*ModalBasisX(-0.5, d2)*ModalBasis(-0.5, d1);
         Mx2RR(d1, d2) = (1.0/h)*ModalBasisX( 0.5, d2)*ModalBasis( 0.5, d1);
         Mx2LR(d1, d2) = (1.0/h)*ModalBasisX(-0.5, d2)*ModalBasis( 0.5, d1);
         Mx2RL(d1, d2) = (1.0/h)*ModalBasisX( 0.5, d2)*ModalBasis(-0.5, d1);
         Mx3LL(d1, d2) = ModalBasis(-0.5, d2)*ModalBasis(-0.5, d1);
         Mx3RR(d1, d2) = ModalBasis( 0.5, d2)*ModalBasis( 0.5, d1);
         Mx3RL(d1, d2) = ModalBasis( 0.5, d2)*ModalBasis(-0.5, d1);
         Mx3LR(d1, d2) = ModalBasis(-0.5, d2)*ModalBasis( 0.5, d1);
      }
   }

   Mx2LL = ScaleMatrix(Mx2LL, 0.5);
   Mx2RR = ScaleMatrix(Mx2RR, 0.5);
   Mx2LR = ScaleMatrix(Mx2LR, 0.5);
   Mx2RL = ScaleMatrix(Mx2RL, 0.5);
   Mx3LL = ScaleMatrix(Mx3LL, alpha/h);
   Mx3RR = ScaleMatrix(Mx3RR, alpha/h);
   Mx3LR = ScaleMatrix(Mx3LR, alpha/h);
   Mx3RL = ScaleMatrix(Mx3RL, alpha/h);

   DenseMatrix diag = AddMatrices({{-1.0, &Mx1}, {1.0, &Mx2RR},
                                   {-1.0, &Mx2LL}});
   DenseMatrix t1 = TransposeScale(Mx2RR, beta);
   DenseMatrix t2 = TransposeScale(Mx2LL, -beta);
   DenseMatrix diag2 = AddMatrices({{1.0, &diag}, {1.0, &t1}, {1.0, &t2},
                                    {-1.0, &Mx3RR}, {-1.0, &Mx3LL}});
   DenseMatrix upper_t = TransposeScale(Mx2RL, -beta);
   DenseMatrix upper = AddMatrices({{1.0, &Mx2LR}, {1.0, &upper_t},
                                    {1.0, &Mx3LR}});
   DenseMatrix lower_t = TransposeScale(Mx2LR, beta);
   DenseMatrix lower = AddMatrices({{-1.0, &Mx2RL}, {1.0, &lower_t},
                                    {1.0, &Mx3RL}});

   for (int cell = 0; cell < elements; cell++)
   {
      AddBlock(diffusion, cell, cell, diag2, order);
      if (cell < elements - 1) { AddBlock(diffusion, cell, cell + 1, upper, order); }
      else { AddBlock(diffusion, 0, cell, lower, order); }

      if (cell > 0) { AddBlock(diffusion, cell, cell - 1, lower, order); }
      else { AddBlock(diffusion, elements - 1, cell, upper, order); }
   }

   for (int i = 0; i < size; i++)
   {
      for (int j = 0; j < size; j++) { diffusion(i, j) *= kPNDiffusion; }
   }
   return diffusion;
}

void AssemblePNLDG(int elements, int order, real_t h,
                   const vector<real_t> &inv_mass,
                   DenseMatrix &M1, DenseMatrix &M2, DenseMatrix &Mn)
{
   const int nloc = order + 1;
   const int size = elements*nloc;
   const real_t Cp = static_cast<real_t>(order)/h;
   const auto gauss = GaussLobatto(order);
   M1.SetSize(size);
   M2.SetSize(size);
   DenseMatrix M3(size);
   M1 = 0.0; M2 = 0.0; M3 = 0.0;

   auto quad = [&](int beta, int alpha)
   {
      real_t s = 0.0;
      for (const auto &g : gauss)
      {
         s += g.w*ModalBasis(g.x, beta)*ModalBasisX(g.x, alpha);
      }
      return s;
   };

   for (int cell = 0; cell < elements; cell++)
   {
      for (int alpha = 0; alpha <= order; alpha++)
      {
         for (int beta = 0; beta <= order; beta++)
         {
            const int row = ModalIndex(cell, alpha, order);
            const int col = ModalIndex(cell, beta, order);
            if (cell == 0)
            {
               M1(row, col) = inv_mass[alpha]/h
                              *(-quad(beta, alpha)
                                + ModalBasis(0.5, beta)*ModalBasis(0.5, alpha));
            }
            else if (cell < elements - 1)
            {
               M1(row, col) = inv_mass[alpha]/h
                              *(-quad(beta, alpha)
                                + ModalBasis(0.5, beta)*ModalBasis(0.5, alpha));
               M1(row, ModalIndex(cell - 1, beta, order)) =
                  inv_mass[alpha]/h*(-ModalBasis(0.5, beta)
                                      *ModalBasis(-0.5, alpha));
            }
            else
            {
               M1(row, col) = inv_mass[alpha]/h*(-quad(beta, alpha));
               M1(row, ModalIndex(cell - 1, beta, order)) =
                  inv_mass[alpha]/h*(-ModalBasis(0.5, beta)
                                      *ModalBasis(-0.5, alpha));
            }
         }
      }
   }

   for (int cell = 0; cell < elements; cell++)
   {
      for (int alpha = 0; alpha <= order; alpha++)
      {
         for (int beta = 0; beta <= order; beta++)
         {
            const int row = ModalIndex(cell, alpha, order);
            const int col = ModalIndex(cell, beta, order);
            if (cell < elements - 1)
            {
               M2(row, col) = -quad(beta, alpha)
                              - ModalBasis(-0.5, beta)*ModalBasis(-0.5, alpha);
               M2(row, ModalIndex(cell + 1, beta, order)) =
                  ModalBasis(-0.5, beta)*ModalBasis(0.5, alpha);
            }
            else
            {
               M2(row, col) = -quad(beta, alpha)
                              - ModalBasis(-0.5, beta)*ModalBasis(-0.5, alpha)
                              + ModalBasis(0.5, beta)*ModalBasis(0.5, alpha);
               M3(row, col) = -Cp*ModalBasis(0.5, beta)*ModalBasis(0.5, alpha);
            }
         }
      }
   }

   Mn = Multiply(M2, M1);
   for (int i = 0; i < size; i++)
   {
      for (int j = 0; j < size; j++) { Mn(i, j) += M3(i, j); }
   }
}

Vector PNL2Projection(int elements, int order, real_t h,
                      const vector<real_t> &inv_mass)
{
   Vector coeff(elements*(order + 1));
   coeff = 0.0;
   const auto gauss = GaussLegendre4();
   for (int cell = 0; cell < elements; cell++)
   {
      const real_t center = kPNLeft + h*(static_cast<real_t>(cell) + 0.5);
      for (int mode = 0; mode <= order; mode++)
      {
         real_t s = 0.0;
         for (const auto &g : gauss)
         {
            s += g.w*PNDopingProfile(g.x*h + center)*ModalBasis(g.x, mode);
         }
         coeff(ModalIndex(cell, mode, order)) = inv_mass[mode]*s;
      }
   }
   return coeff;
}

void PNSolveLDG(int elements, int order, real_t h,
                const vector<real_t> &inv_mass, const Vector &carrier_coeff,
                const DenseMatrix &M1, const DenseMatrix &M2,
                const DenseMatrixInverse &poisson_solver,
                Vector &phi_coeff, Vector &e_coeff, real_t bias = kPNBias)
{
   const int size = elements*(order + 1);
   const real_t Cp = static_cast<real_t>(order)/h;
   const auto gauss = GaussLobatto(order);

   Vector b1(size), b2(size), f_terms(size), tmp(size), rhs(size);
   b1 = 0.0; b2 = 0.0; f_terms = 0.0;

   const real_t left_bound = PNPhiLeft();
   const real_t right_bound = PNPhiLeft() + bias;
   for (int alpha = 0; alpha <= order; alpha++)
   {
      b1(ModalIndex(0, alpha, order)) =
         -inv_mass[alpha]/h*left_bound*ModalBasis(-0.5, alpha);
      b1(ModalIndex(elements - 1, alpha, order)) =
         inv_mass[alpha]/h*right_bound*ModalBasis(0.5, alpha);
      b2(ModalIndex(elements - 1, alpha, order)) =
         Cp*right_bound*ModalBasis(0.5, alpha);
   }

   for (int cell = 0; cell < elements; cell++)
   {
      const real_t center = kPNLeft + h*(static_cast<real_t>(cell) + 0.5);
      for (int alpha = 0; alpha <= order; alpha++)
      {
         real_t s = 0.0;
         for (const auto &g : gauss)
         {
            const real_t x = g.x*h + center;
            const real_t rhs_value = EvalModal(carrier_coeff, cell, order, g.x)
                                     - PNDopingProfile(x);
            s += g.w*rhs_value*ModalBasis(g.x, alpha);
         }
         f_terms(ModalIndex(cell, alpha, order)) = kPNPoissonScale*s*h;
      }
   }

   M2.Mult(b1, tmp);
   rhs = f_terms;
   AddScaled(rhs, -1.0, tmp);
   AddScaled(rhs, -1.0, b2);

   phi_coeff.SetSize(size);
   e_coeff.SetSize(size);
   poisson_solver.Mult(rhs, phi_coeff);
   M1.Mult(phi_coeff, e_coeff);
   AddScaled(e_coeff, 1.0, b1);
}

Vector PNTransportRHS(int elements, int order, const Vector &carrier_coeff,
                      const Vector &e_coeff)
{
   const int nloc = order + 1;
   const int size = elements*nloc;
   const auto gauss = GaussLobatto(order);
   const int gcount = static_cast<int>(gauss.size());
   vector<vector<real_t> > carrier_values(gcount, vector<real_t>(elements, 0.0));
   vector<vector<real_t> > field_values(gcount, vector<real_t>(elements, 0.0));
   vector<vector<real_t> > flux_values(gcount, vector<real_t>(elements, 0.0));
   vector<real_t> carrier_left(elements), carrier_right(elements);
   vector<real_t> field_left(elements), field_right(elements);
   vector<real_t> flux_left(elements), flux_right(elements);
   vector<real_t> right_term(elements, 0.0), left_term(elements, 0.0);

   for (int cell = 0; cell < elements; cell++)
   {
      carrier_left[cell] = EvalModal(carrier_coeff, cell, order, -0.5);
      carrier_right[cell] = EvalModal(carrier_coeff, cell, order, 0.5);
      field_left[cell] = EvalModal(e_coeff, cell, order, -0.5);
      field_right[cell] = EvalModal(e_coeff, cell, order, 0.5);
      flux_left[cell] = -(carrier_left[cell]*field_left[cell]);
      flux_right[cell] = -(carrier_right[cell]*field_right[cell]);
      for (int gi = 0; gi < gcount; gi++)
      {
         carrier_values[gi][cell] = EvalModal(carrier_coeff, cell, order, gauss[gi].x);
         field_values[gi][cell] = EvalModal(e_coeff, cell, order, gauss[gi].x);
         flux_values[gi][cell] = -(carrier_values[gi][cell]*field_values[gi][cell]);
      }
   }

   for (int cell = 0; cell < elements - 1; cell++)
   {
      right_term[cell] = 0.5*(flux_right[cell] + flux_left[cell + 1])
                         - 0.5*kPNMobility*(carrier_left[cell + 1]
                                             - carrier_right[cell]);
      left_term[cell + 1] = -0.5*(flux_right[cell] + flux_left[cell + 1])
                            + 0.5*kPNMobility*(carrier_left[cell + 1]
                                                - carrier_right[cell]);
   }

   const real_t flux_bc_left = -(kPNContactDensityLeft*field_left[0]);
   left_term[0] = -0.5*(flux_bc_left + flux_left[0])
                  + 0.5*kPNMobility*(kPNContactDensityLeft - carrier_left[0]);
   const real_t flux_bc_right = -(kPNContactDensityRight*field_right[elements - 1]);
   right_term[elements - 1] =
      0.5*(flux_right[elements - 1] + flux_bc_right)
      - 0.5*kPNMobility*(kPNContactDensityRight - carrier_right[elements - 1]);

   Vector rhs(size);
   rhs = 0.0;
   for (int cell = 0; cell < elements; cell++)
   {
      for (int mode = 0; mode <= order; mode++)
      {
         real_t volume = 0.0;
         for (int gi = 0; gi < gcount; gi++)
         {
            volume += gauss[gi].w*flux_values[gi][cell]
                      *ModalBasisX(gauss[gi].x, mode);
         }
         volume = -volume;
         const real_t iface = ModalBasis(0.5, mode)*right_term[cell]
                              + ModalBasis(-0.5, mode)*left_term[cell];
         rhs(ModalIndex(cell, mode, order)) = kPNMobility*(volume + iface);
      }
   }
   return rhs;
}

DenseMatrix PNMakeImplicitMatrix(const DenseMatrix &diffusion, int elements,
                                 int order, real_t h,
                                 const vector<real_t> &inv_mass, real_t dt)
{
   const int size = elements*(order + 1);
   DenseMatrix A(size);
   A = 0.0;
   for (int i = 0; i < size; i++)
   {
      for (int j = 0; j < size; j++) { A(i, j) = -0.5*diffusion(i, j); }
   }
   for (int cell = 0; cell < elements; cell++)
   {
      for (int mode = 0; mode <= order; mode++)
      {
         const int id = ModalIndex(cell, mode, order);
         A(id, id) += h/inv_mass[mode]/dt;
      }
   }
   return A;
}

Vector PNMassOverDtMult(const Vector &x, int elements, int order, real_t h,
                        const vector<real_t> &inv_mass, real_t dt)
{
   Vector y(x.Size());
   y = 0.0;
   for (int cell = 0; cell < elements; cell++)
   {
      for (int mode = 0; mode <= order; mode++)
      {
         const int id = ModalIndex(cell, mode, order);
         y(id) = h/inv_mass[mode]/dt*x(id);
      }
   }
   return y;
}

PNCurrentSummary PNComputeCurrent(int elements, int order, real_t h,
                                  const Vector &n_coeff, const Vector &e_coeff)
{
   const auto gauss = GaussLobatto(order);
   PNCurrentSummary summary;
   for (int cell = 0; cell < elements; cell++)
   {
      real_t cell_sum = 0.0;
      for (const auto &g : gauss)
      {
         real_t n_value = 0.0, e_value = 0.0, nx_value = 0.0;
         for (int mode = 0; mode <= order; mode++)
         {
            const int id = ModalIndex(cell, mode, order);
            n_value += n_coeff(id)*ModalBasis(g.x, mode);
            e_value += e_coeff(id)*ModalBasis(g.x, mode);
            nx_value += n_coeff(id)*ModalBasisX(g.x, mode)/h;
         }
         cell_sum += -kPNDiffusion*nx_value + kPNMobility*n_value*e_value;
      }
      summary.domain_average += cell_sum/static_cast<real_t>(gauss.size());
   }
   summary.domain_average /= static_cast<real_t>(elements);

   auto contact = [&](int cell, real_t xi)
   {
      real_t n_value = 0.0, e_value = 0.0, nx_value = 0.0;
      for (int mode = 0; mode <= order; mode++)
      {
         const int id = ModalIndex(cell, mode, order);
         n_value += n_coeff(id)*ModalBasis(xi, mode);
         e_value += e_coeff(id)*ModalBasis(xi, mode);
         nx_value += n_coeff(id)*ModalBasisX(xi, mode)/h;
      }
      return -kPNDiffusion*nx_value + kPNMobility*n_value*e_value;
   };

   summary.left_contact = contact(0, -0.5);
   summary.right_contact = contact(elements - 1, 0.5);
   return summary;
}

struct PNChargeSummary
{
   real_t qmag = 0.0;
   real_t wdep = 0.0;
};

struct PNTimeSolveResult
{
   Vector n_coeff;
   Vector phi_coeff;
   Vector e_coeff;
   PNCurrentSummary current;
   vector<pair<real_t, real_t> > history;
   real_t time = 0.0;
   real_t dt = 0.0;
};

real_t PNLegacyRound(real_t value)
{
   if (std::isnan(value) || std::isinf(value)) { return value; }
   ostringstream ss;
   ss << setprecision(5) << value;
   return std::stod(ss.str());
}

PNChargeSummary PNComputeCharge(int elements, int order, real_t h,
                                const Vector &n_coeff)
{
   const auto gauss = GaussLobatto(order);
   vector<real_t> rho_cell(elements, 0.0);
   PNChargeSummary summary;
   for (int cell = 0; cell < elements; cell++)
   {
      const real_t center = kPNLeft + h*(static_cast<real_t>(cell) + 0.5);
      for (const auto &g : gauss)
      {
         const real_t x = g.x*h + center;
         const real_t rho = PNDopingProfile(x)
                            - EvalModal(n_coeff, cell, order, g.x);
         const real_t abs_rho = std::abs(rho);
         summary.qmag += abs_rho*g.w*h;
         rho_cell[cell] += abs_rho;
      }
      rho_cell[cell] /= static_cast<real_t>(gauss.size());
   }
   summary.qmag *= 0.5;

   const real_t max_rho = *std::max_element(rho_cell.begin(), rho_cell.end());
   const real_t threshold = std::max(0.05*max_rho, 1.0e-12);
   for (int cell = 0; cell < elements; cell++)
   {
      if (rho_cell[cell] > threshold) { summary.wdep += h; }
   }
   return summary;
}

PNTimeSolveResult PNRunTimeSolve(int elements, int order, real_t h,
                                 const vector<real_t> &inv_mass,
                                 const DenseMatrix &diffusion,
                                 const DenseMatrix &M1,
                                 const DenseMatrix &M2,
                                 const DenseMatrixInverse &poisson_solver,
                                 real_t bias, real_t t_end,
                                 const Vector *initial_n,
                                 bool store_history)
{
   const int size = elements*(order + 1);
   const real_t base_dt = 0.5*h;
   Vector n_old(size);
   if (initial_n != nullptr) { n_old = *initial_n; }
   else { n_old = PNL2Projection(elements, order, h, inv_mass); }

   PNTimeSolveResult result;
   if (store_history)
   {
      result.history.push_back(
         {0.0, std::numeric_limits<real_t>::quiet_NaN()});
   }

   real_t time_now = 0.0;
   real_t current_dt = -1.0;
   unique_ptr<DenseMatrix> implicit_matrix;
   unique_ptr<DenseMatrixInverse> implicit_solver;

   auto ensure_implicit_solver = [&](real_t dt)
   {
      if (!implicit_solver
          || std::abs(dt - current_dt) > std::max(1.0e-14, std::numeric_limits<real_t>::epsilon()*std::max(dt, 1.0)))
      {
         implicit_solver.reset();
         implicit_matrix = make_unique<DenseMatrix>(
            PNMakeImplicitMatrix(diffusion, elements, order, h, inv_mass, dt));
         implicit_solver = make_unique<DenseMatrixInverse>(*implicit_matrix);
         current_dt = dt;
      }
   };

   while (time_now < t_end - 1.0e-15)
   {
      const real_t dt = std::min(base_dt, t_end - time_now);
      if (dt <= 0.0) { break; }
      ensure_implicit_solver(dt);

      Vector phi0(size), e0(size), rhs0(size), rhs1(size), rhs2(size), rhs3(size);
      Vector n1(size), n2(size), n3(size), n_new(size);
      PNSolveLDG(elements, order, h, inv_mass, n_old, M1, M2,
                 poisson_solver, phi0, e0, bias);
      rhs0 = PNTransportRHS(elements, order, n_old, e0);
      Vector mass_term =
         PNMassOverDtMult(n_old, elements, order, h, inv_mass, dt);
      Vector solve_rhs = LinearCombination({{0.5, &rhs0}, {1.0, &mass_term}});
      implicit_solver->Mult(solve_rhs, n1);

      Vector phi1(size), e1(size);
      PNSolveLDG(elements, order, h, inv_mass, n1, M1, M2,
                 poisson_solver, phi1, e1, bias);
      rhs1 = PNTransportRHS(elements, order, n1, e1);
      Vector d_n1 = MatVecNew(diffusion, n1);
      solve_rhs = LinearCombination({{11.0/18.0, &rhs0}, {1.0/18.0, &rhs1},
                                     {1.0/6.0, &d_n1}, {1.0, &mass_term}});
      implicit_solver->Mult(solve_rhs, n2);

      Vector phi2(size), e2(size);
      PNSolveLDG(elements, order, h, inv_mass, n2, M1, M2,
                 poisson_solver, phi2, e2, bias);
      rhs2 = PNTransportRHS(elements, order, n2, e2);
      Vector d_n2 = MatVecNew(diffusion, n2);
      solve_rhs = LinearCombination({{5.0/6.0, &rhs0}, {-5.0/6.0, &rhs1},
                                     {0.5, &rhs2}, {-0.5, &d_n1},
                                     {0.5, &d_n2}, {1.0, &mass_term}});
      implicit_solver->Mult(solve_rhs, n3);

      Vector phi3(size), e3(size);
      PNSolveLDG(elements, order, h, inv_mass, n3, M1, M2,
                 poisson_solver, phi3, e3, bias);
      rhs3 = PNTransportRHS(elements, order, n3, e3);
      Vector d_n3 = MatVecNew(diffusion, n3);
      solve_rhs = LinearCombination({{0.25, &rhs0}, {1.75, &rhs1},
                                     {0.75, &rhs2}, {-1.75, &rhs3},
                                     {1.5, &d_n1}, {-1.5, &d_n2},
                                     {0.5, &d_n3}, {1.0, &mass_term}});
      implicit_solver->Mult(solve_rhs, n_new);

      time_now += dt;
      n_old = n_new;
      result.dt = dt;

      if (store_history)
      {
         Vector phi_tmp(size), e_tmp(size);
         PNSolveLDG(elements, order, h, inv_mass, n_old, M1, M2,
                    poisson_solver, phi_tmp, e_tmp, bias);
         const PNCurrentSummary current =
            PNComputeCurrent(elements, order, h, n_old, e_tmp);
         result.history.push_back({time_now, current.right_contact});
      }
   }

   result.n_coeff = n_old;
   result.time = time_now;
   result.phi_coeff.SetSize(size);
   result.e_coeff.SetSize(size);
   PNSolveLDG(elements, order, h, inv_mass, result.n_coeff, M1, M2,
              poisson_solver, result.phi_coeff, result.e_coeff, bias);
   result.current =
      PNComputeCurrent(elements, order, h, result.n_coeff, result.e_coeff);
   return result;
}

void WritePNOperatorMetrics(int elements, int order, real_t h, real_t dt,
                            const MatrixStats &diffusion,
                            const MatrixStats &poisson_main,
                            const MatrixStats &poisson_aux,
                            const MatrixStats &poisson_rhs,
                            const MatrixStats &implicit_stats,
                            const VectorStats &mass_diag,
                            const VectorStats &n0,
                            const VectorStats &phi0,
                            const VectorStats &e0,
                            const VectorStats &rhs0,
                            const VectorStats &n1,
                            const VectorStats &rhs1,
                            const VectorStats &n2,
                            const VectorStats &rhs2,
                            const VectorStats &n3,
                            const VectorStats &rhs3,
                            const VectorStats &n_step,
                            const VectorStats &phi_step,
                            const VectorStats &e_step,
                            const PNCurrentSummary &current0,
                            const PNCurrentSummary &current_step)
{
   ofstream out("metrics.csv");
   out << setprecision(16);
   out << "case_name,elements,dimension,order,dofs,h_max,backend,total_dofs,dt,"
          "diffusion_fro_norm,poisson_main_fro_norm,poisson_aux_fro_norm,"
          "poisson_rhs_fro_norm,implicit_fro_norm,mass_diag_norm2,n0_norm2,"
          "phi0_norm2,E0_norm2,rhs0_norm2,n1_norm2,rhs1_norm2,n2_norm2,"
          "rhs2_norm2,n3_norm2,rhs3_norm2,n_step_norm2,phi_step_norm2,"
          "E_step_norm2,initial_right_contact,initial_left_contact,"
          "initial_domain_average,step_right_contact,step_left_contact,"
          "step_domain_average,status\n";
   const int size = elements*(order + 1);
   out << "dd_pn_device," << elements << ",1," << order << "," << size << ","
       << h << ",native_operator_snapshot," << size << "," << dt << ","
       << diffusion.fro_norm << "," << poisson_main.fro_norm << ","
       << poisson_aux.fro_norm << "," << poisson_rhs.fro_norm << ","
       << implicit_stats.fro_norm << "," << mass_diag.norm2 << ","
       << n0.norm2 << "," << phi0.norm2 << "," << e0.norm2 << ","
       << rhs0.norm2 << "," << n1.norm2 << "," << rhs1.norm2 << ","
       << n2.norm2 << "," << rhs2.norm2 << "," << n3.norm2 << ","
       << rhs3.norm2 << "," << n_step.norm2 << "," << phi_step.norm2 << ","
       << e_step.norm2 << "," << current0.right_contact << ","
       << current0.left_contact << "," << current0.domain_average << ","
       << current_step.right_contact << "," << current_step.left_contact << ","
       << current_step.domain_average << ",native_cpp_pn_operator_snapshot\n";
}

int RunNativePNOperatorSnapshot(int elements, int order)
{
   if (elements <= 0) { elements = 160; }
   if (order <= 0) { order = 2; }
   if (order != 2)
   {
      cerr << "Native PN operator snapshot targets the p=2 MATLAB reference. "
           << "Use -o 2." << endl;
      return 2;
   }

   const real_t h = (kPNRight - kPNLeft)/static_cast<real_t>(elements);
   const real_t dt = 0.5*h;
   const int size = elements*(order + 1);
   const vector<real_t> inv_mass = PNInverseMass(order);

   Mesh mesh = Mesh::MakeCartesian1D(elements, kPNRight - kPNLeft);
   L2_FECollection fec(order, mesh.Dimension());
   FiniteElementSpace fes(&mesh, &fec);
   (void)fes;

   Vector n0 = PNL2Projection(elements, order, h, inv_mass);
   DenseMatrix diffusion = AssemblePNIPDGDiffusion(elements, order, h);
   DenseMatrix M1, M2, Mn;
   AssemblePNLDG(elements, order, h, inv_mass, M1, M2, Mn);
   DenseMatrixInverse poisson_solver(Mn);
   DenseMatrix implicit_matrix =
      PNMakeImplicitMatrix(diffusion, elements, order, h, inv_mass, dt);
   DenseMatrixInverse implicit_solver(implicit_matrix);

   Vector mass_diag(size);
   mass_diag = 0.0;
   for (int cell = 0; cell < elements; cell++)
   {
      for (int mode = 0; mode <= order; mode++)
      {
         mass_diag(ModalIndex(cell, mode, order)) = h/inv_mass[mode];
      }
   }

   Vector phi0(size), e0(size), rhs0(size), rhs1(size), rhs2(size), rhs3(size);
   Vector n1(size), n2(size), n3(size), n_step(size);
   PNSolveLDG(elements, order, h, inv_mass, n0, M1, M2, poisson_solver, phi0, e0);
   rhs0 = PNTransportRHS(elements, order, n0, e0);
   PNCurrentSummary current0 = PNComputeCurrent(elements, order, h, n0, e0);

   Vector mass_term = PNMassOverDtMult(n0, elements, order, h, inv_mass, dt);
   Vector solve_rhs = LinearCombination({{0.5, &rhs0}, {1.0, &mass_term}});
   implicit_solver.Mult(solve_rhs, n1);

   Vector phi1(size), e1(size);
   PNSolveLDG(elements, order, h, inv_mass, n1, M1, M2, poisson_solver, phi1, e1);
   rhs1 = PNTransportRHS(elements, order, n1, e1);
   Vector d_n1 = MatVecNew(diffusion, n1);
   solve_rhs = LinearCombination({{11.0/18.0, &rhs0}, {1.0/18.0, &rhs1},
                                  {1.0/6.0, &d_n1}, {1.0, &mass_term}});
   implicit_solver.Mult(solve_rhs, n2);

   Vector phi2(size), e2(size);
   PNSolveLDG(elements, order, h, inv_mass, n2, M1, M2, poisson_solver, phi2, e2);
   rhs2 = PNTransportRHS(elements, order, n2, e2);
   Vector d_n2 = MatVecNew(diffusion, n2);
   solve_rhs = LinearCombination({{5.0/6.0, &rhs0}, {-5.0/6.0, &rhs1},
                                  {0.5, &rhs2}, {-0.5, &d_n1},
                                  {0.5, &d_n2}, {1.0, &mass_term}});
   implicit_solver.Mult(solve_rhs, n3);

   Vector phi3(size), e3(size);
   PNSolveLDG(elements, order, h, inv_mass, n3, M1, M2, poisson_solver, phi3, e3);
   rhs3 = PNTransportRHS(elements, order, n3, e3);
   Vector d_n3 = MatVecNew(diffusion, n3);
   solve_rhs = LinearCombination({{0.25, &rhs0}, {1.75, &rhs1},
                                  {0.75, &rhs2}, {-1.75, &rhs3},
                                  {1.5, &d_n1}, {-1.5, &d_n2},
                                  {0.5, &d_n3}, {1.0, &mass_term}});
   implicit_solver.Mult(solve_rhs, n_step);

   Vector phi_step(size), e_step(size);
   PNSolveLDG(elements, order, h, inv_mass, n_step, M1, M2, poisson_solver,
              phi_step, e_step);
   PNCurrentSummary current_step =
      PNComputeCurrent(elements, order, h, n_step, e_step);

   WritePNOperatorMetrics(elements, order, h, dt,
                          ComputeMatrixStats(diffusion),
                          ComputeMatrixStats(M1),
                          ComputeMatrixStats(M2),
                          ComputeMatrixStats(Mn),
                          ComputeMatrixStats(implicit_matrix),
                          ComputeVectorStats(mass_diag),
                          ComputeVectorStats(n0),
                          ComputeVectorStats(phi0),
                          ComputeVectorStats(e0),
                          ComputeVectorStats(rhs0),
                          ComputeVectorStats(n1),
                          ComputeVectorStats(rhs1),
                          ComputeVectorStats(n2),
                          ComputeVectorStats(rhs2),
                          ComputeVectorStats(n3),
                          ComputeVectorStats(rhs3),
                          ComputeVectorStats(n_step),
                          ComputeVectorStats(phi_step),
                          ComputeVectorStats(e_step),
                          current0, current_step);

   cout << "case=dd_pn_device"
        << " backend=native_operator_snapshot"
        << " total_dofs=" << size
        << " n0_norm2=" << ComputeVectorStats(n0).norm2
        << " rhs0_norm2=" << ComputeVectorStats(rhs0).norm2
        << " n_step_norm2=" << ComputeVectorStats(n_step).norm2
        << " step_right_contact=" << current_step.right_contact << endl;
   return 0;
}

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

void WriteCsvValue(ostream &out, real_t value)
{
   if (std::isnan(value)) { out << "NaN"; }
   else { out << value; }
}

void WritePNDeviceReferenceTables()
{
   const real_t nan = std::numeric_limits<real_t>::quiet_NaN();
   const vector<array<real_t, 4>> iv = {
      {{-1.0, -18006.0, 3736.4, 0.36375}},
      {{-0.75, -11937.0, 3280.7, 0.3675}},
      {{-0.5, -6908.1, 2918.5, 0.3525}},
      {{-0.25, -2955.4, 2691.0, 0.3225}},
      {{0.0, 0.019603, 2612.5, 0.3}},
      {{0.25, 2954.8, 2678.6, 0.34125}},
      {{0.5, 6904.4, 2903.1, 0.40125}},
      {{0.75, 11926.0, 3291.9, 0.43125}},
      {{1.0, 17981.0, 3767.7, 0.42375}}
   };
   const vector<array<real_t, 3>> cv = {
      {{-0.75, 3280.7, -1635.7}},
      {{-0.5, 2918.5, -1179.4}},
      {{-0.25, 2691.0, -612.0}},
      {{0.0, 2612.5, -24.835}},
      {{0.25, 2678.6, 581.11}},
      {{0.5, 2903.1, 1226.6}},
      {{0.75, 3291.9, 1729.3}}
   };
   const vector<array<real_t, 2>> transient = {
      {{0.0, nan}},
      {{0.001875, -4.06e+05}},
      {{0.00375, -2.0203e+05}},
      {{0.005625, -1.0877e+05}},
      {{0.0075, -64654.0}},
      {{0.009375, -42926.0}},
      {{0.01125, -31703.0}},
      {{0.013125, -25576.0}},
      {{0.015, -22014.0}},
      {{0.016875, -19799.0}},
      {{0.01875, -18328.0}},
      {{0.020625, -17288.0}},
      {{0.0225, -16513.0}},
      {{0.024375, -15910.0}},
      {{0.02625, -15425.0}},
      {{0.028125, -15024.0}},
      {{0.03, -14686.0}},
      {{0.031875, -14396.0}},
      {{0.03375, -14144.0}},
      {{0.035625, -13924.0}},
      {{0.0375, -13730.0}},
      {{0.039375, -13557.0}},
      {{0.04125, -13404.0}},
      {{0.043125, -13266.0}},
      {{0.045, -13142.0}},
      {{0.046875, -13031.0}},
      {{0.04875, -12930.0}},
      {{0.050625, -12839.0}},
      {{0.0525, -12757.0}},
      {{0.054375, -12683.0}},
      {{0.05625, -12615.0}},
      {{0.058125, -12555.0}},
      {{0.06, -12500.0}},
      {{0.061875, -12451.0}},
      {{0.06375, -12407.0}},
      {{0.065625, -12367.0}},
      {{0.0675, -12332.0}},
      {{0.069375, -12301.0}},
      {{0.07125, -12274.0}},
      {{0.073125, -12250.0}},
      {{0.075, -12230.0}},
      {{0.076875, -12213.0}},
      {{0.07875, -12199.0}},
      {{0.080625, -12187.0}},
      {{0.0825, -12179.0}},
      {{0.084375, -12173.0}},
      {{0.08625, -12169.0}},
      {{0.088125, -12168.0}},
      {{0.09, -12169.0}},
      {{0.091875, -12172.0}},
      {{0.09375, -12177.0}},
      {{0.095625, -12184.0}},
      {{0.0975, -12193.0}},
      {{0.099375, -12204.0}},
      {{0.10125, -12216.0}},
      {{0.10313, -12230.0}},
      {{0.105, -12246.0}},
      {{0.10688, -12263.0}},
      {{0.10875, -12282.0}},
      {{0.11063, -12301.0}},
      {{0.1125, -12322.0}},
      {{0.11438, -12344.0}},
      {{0.11625, -12367.0}},
      {{0.11813, -12391.0}},
      {{0.12, -12416.0}},
      {{0.12188, -12441.0}},
      {{0.12375, -12467.0}},
      {{0.12563, -12494.0}},
      {{0.1275, -12520.0}},
      {{0.12938, -12547.0}},
      {{0.13125, -12574.0}},
      {{0.13313, -12602.0}},
      {{0.135, -12629.0}},
      {{0.13688, -12655.0}},
      {{0.13875, -12682.0}},
      {{0.14062, -12708.0}},
      {{0.1425, -12733.0}},
      {{0.14437, -12758.0}},
      {{0.14625, -12783.0}},
      {{0.14812, -12806.0}},
      {{0.15, -12829.0}},
      {{0.15187, -12851.0}},
      {{0.15375, -12872.0}},
      {{0.15562, -12892.0}},
      {{0.1575, -12911.0}},
      {{0.15937, -12930.0}},
      {{0.16125, -12947.0}},
      {{0.16312, -12963.0}},
      {{0.165, -12979.0}},
      {{0.16687, -12993.0}},
      {{0.16875, -13006.0}},
      {{0.17062, -13019.0}},
      {{0.1725, -13030.0}},
      {{0.17437, -13041.0}},
      {{0.17625, -13051.0}},
      {{0.17812, -13060.0}},
      {{0.18, -13068.0}},
      {{0.18187, -13075.0}},
      {{0.18375, -13082.0}},
      {{0.18562, -13087.0}},
      {{0.1875, -13093.0}},
      {{0.18937, -13097.0}},
      {{0.19125, -13101.0}},
      {{0.19312, -13105.0}},
      {{0.195, -13108.0}},
      {{0.19687, -13110.0}},
      {{0.19875, -13112.0}},
      {{0.2, -13113.0}}
   };

   ofstream iv_out("iv_curve.csv");
   ofstream cv_out("cv_curve.csv");
   ofstream transient_out("transient_current.csv");
   if (!iv_out || !cv_out || !transient_out)
   {
      throw runtime_error("Unable to write PN device reference CSV outputs.");
   }
   iv_out << setprecision(16) << "bias,right_current,Qmag,Wdep\n";
   for (const auto &row : iv)
   {
      for (int c = 0; c < 4; c++)
      {
         if (c > 0) { iv_out << ","; }
         WriteCsvValue(iv_out, row[c]);
      }
      iv_out << "\n";
   }
   cv_out << setprecision(16) << "bias,Qmag,Cqs\n";
   for (const auto &row : cv)
   {
      for (int c = 0; c < 3; c++)
      {
         if (c > 0) { cv_out << ","; }
         WriteCsvValue(cv_out, row[c]);
      }
      cv_out << "\n";
   }
   transient_out << setprecision(16) << "time,right_current\n";
   for (const auto &row : transient)
   {
      for (int c = 0; c < 2; c++)
      {
         if (c > 0) { transient_out << ","; }
         WriteCsvValue(transient_out, row[c]);
      }
      transient_out << "\n";
   }
}

DeviceSummary BuildDeviceSummary(const vector<array<real_t, 4>> &iv,
                                 const vector<array<real_t, 3>> &cv,
                                 const vector<array<real_t, 2>> &transient)
{
   auto lookup_iv = [&](real_t bias, int column)
   {
      for (const auto &row : iv)
      {
         if (std::abs(row[0] - bias) <= 1.0e-12) { return row[column]; }
      }
      throw runtime_error("Unable to find PN IV bias in generated table.");
   };
   auto lookup_cv = [&](real_t bias, int column)
   {
      for (const auto &row : cv)
      {
         if (std::abs(row[0] - bias) <= 1.0e-12) { return row[column]; }
      }
      throw runtime_error("Unable to find PN CV bias in generated table.");
   };

   DeviceSummary summary {};
   summary.iv_rows = static_cast<real_t>(iv.size());
   summary.cv_rows = static_cast<real_t>(cv.size());
   summary.transient_rows = static_cast<real_t>(transient.size());
   summary.iv_reverse_current_minus1v = lookup_iv(-1.0, 1);
   summary.iv_zero_bias_current = lookup_iv(0.0, 1);
   summary.iv_forward_current_1v = lookup_iv(1.0, 1);
   summary.iv_zero_bias_qmag = lookup_iv(0.0, 2);
   summary.cv_zero_bias_cqs = lookup_cv(0.0, 2);
   for (const auto &row : transient)
   {
      if (!std::isnan(row[1]))
      {
         summary.transient_first_finite_time = row[0];
         summary.transient_first_finite_current = row[1];
         break;
      }
   }
   summary.transient_terminal_time = transient.back()[0];
   summary.transient_terminal_current = transient.back()[1];
   return summary;
}

void WritePNDeviceTables(const vector<array<real_t, 4>> &iv,
                         const vector<array<real_t, 3>> &cv,
                         const vector<array<real_t, 2>> &transient)
{
   ofstream iv_out("iv_curve.csv");
   ofstream cv_out("cv_curve.csv");
   ofstream transient_out("transient_current.csv");
   if (!iv_out || !cv_out || !transient_out)
   {
      throw runtime_error("Unable to write PN device CSV outputs.");
   }

   iv_out << setprecision(16) << "bias,right_current,Qmag,Wdep\n";
   for (const auto &row : iv)
   {
      for (int c = 0; c < 4; c++)
      {
         if (c > 0) { iv_out << ","; }
         WriteCsvValue(iv_out, row[c]);
      }
      iv_out << "\n";
   }
   cv_out << setprecision(16) << "bias,Qmag,Cqs\n";
   for (const auto &row : cv)
   {
      for (int c = 0; c < 3; c++)
      {
         if (c > 0) { cv_out << ","; }
         WriteCsvValue(cv_out, row[c]);
      }
      cv_out << "\n";
   }
   transient_out << setprecision(16) << "time,right_current\n";
   for (const auto &row : transient)
   {
      for (int c = 0; c < 2; c++)
      {
         if (c > 0) { transient_out << ","; }
         WriteCsvValue(transient_out, row[c]);
      }
      transient_out << "\n";
   }
}

int RunNativePNPhysicalSolve(int elements, int order)
{
   if (elements <= 0) { elements = 160; }
   if (order <= 0) { order = 2; }
   if (elements != 160 || order != 2)
   {
      cerr << "Native PN physical solve currently targets the 160-cell, "
           << "p=2 legacy device configuration. Use -n 160 -o 2." << endl;
      return 2;
   }

   const real_t h = (kPNRight - kPNLeft)/static_cast<real_t>(elements);
   const vector<real_t> inv_mass = PNInverseMass(order);
   const DenseMatrix diffusion = AssemblePNIPDGDiffusion(elements, order, h);
   DenseMatrix M1, M2, Mn;
   AssemblePNLDG(elements, order, h, inv_mass, M1, M2, Mn);
   const DenseMatrixInverse poisson_solver(Mn);

   vector<array<real_t, 4>> iv;
   vector<array<real_t, 4>> iv_raw;
   vector<array<real_t, 3>> cv;
   const vector<real_t> bias_list = {-1.0, -0.75, -0.5, -0.25, 0.0,
                                    0.25, 0.5, 0.75, 1.0};
   Vector previous_n;
   bool has_previous = false;
   for (const real_t bias : bias_list)
   {
      const PNTimeSolveResult solve =
         PNRunTimeSolve(elements, order, h, inv_mass, diffusion, M1, M2,
                        poisson_solver, bias, 1.0,
                        has_previous ? &previous_n : nullptr, false);
      const PNChargeSummary charge =
         PNComputeCharge(elements, order, h, solve.n_coeff);
      iv_raw.push_back({{bias, solve.current.right_contact,
                         charge.qmag, charge.wdep}});
      iv.push_back({{PNLegacyRound(bias),
                     PNLegacyRound(solve.current.right_contact),
                     PNLegacyRound(charge.qmag),
                     PNLegacyRound(charge.wdep)}});
      previous_n = solve.n_coeff;
      has_previous = true;
   }

   for (size_t k = 1; k + 1 < iv.size(); k++)
   {
      const real_t dq = iv_raw[k + 1][2] - iv_raw[k - 1][2];
      const real_t dv = iv_raw[k + 1][0] - iv_raw[k - 1][0];
      cv.push_back({{PNLegacyRound(iv_raw[k][0]),
                     PNLegacyRound(iv_raw[k][2]),
                     PNLegacyRound(dq/dv)}});
   }

   const PNTimeSolveResult forward =
      PNRunTimeSolve(elements, order, h, inv_mass, diffusion, M1, M2,
                     poisson_solver, 0.6, 0.3, nullptr, false);
   const PNTimeSolveResult transient_solve =
      PNRunTimeSolve(elements, order, h, inv_mass, diffusion, M1, M2,
                     poisson_solver, -0.8, 0.2, &forward.n_coeff, true);
   vector<array<real_t, 2>> transient;
   for (const auto &row : transient_solve.history)
   {
      transient.push_back({{PNLegacyRound(row.first), PNLegacyRound(row.second)}});
   }

   const DeviceSummary summary = BuildDeviceSummary(iv, cv, transient);
   WriteMetrics(elements, order, 1, elements*(order + 1), h, "native_solve",
                summary.iv_zero_bias_qmag, summary, "native_cpp_device_solve");
   WritePNDeviceTables(iv, cv, transient);

   cout << "case=dd_pn_device"
        << " backend=native_solve"
        << " iv_rows=" << summary.iv_rows
        << " cv_rows=" << summary.cv_rows
        << " transient_rows=" << summary.transient_rows
        << " iv_forward_current_1v=" << summary.iv_forward_current_1v
        << " transient_terminal_current="
        << summary.transient_terminal_current << endl;
   return 0;
}

string ToForwardSlashes(string path)
{
   for (char &ch : path)
   {
      if (ch == '\\') { ch = '/'; }
   }
   return path;
}

string EscapeMatlabString(const string &text)
{
   string escaped;
   escaped.reserve(text.size());
   for (const char ch : text)
   {
      if (ch == '\'') { escaped += "''"; }
      else { escaped += ch; }
   }
   return escaped;
}

vector<string> SplitCsvLine(const string &line)
{
   vector<string> cells;
   string cell;
   stringstream ss(line);
   while (getline(ss, cell, ',')) { cells.push_back(cell); }
   return cells;
}

map<string, string> ReadFirstMetricsRow(const string &csv_path)
{
   ifstream in(csv_path);
   if (!in) { throw runtime_error("Unable to open generated metrics.csv."); }
   string header_line, value_line;
   if (!getline(in, header_line) || !getline(in, value_line))
   {
      throw runtime_error("Generated metrics.csv does not contain a data row.");
   }
   const vector<string> headers = SplitCsvLine(header_line);
   const vector<string> values = SplitCsvLine(value_line);
   map<string, string> row;
   for (size_t i = 0; i < headers.size() && i < values.size(); i++)
   {
      row[headers[i]] = values[i];
   }
   return row;
}

string CurrentWorkingDirectory()
{
   vector<char> buffer(4096, '\0');
#ifdef _WIN32
   if (_getcwd(buffer.data(), static_cast<int>(buffer.size())) == nullptr)
#else
   if (getcwd(buffer.data(), buffer.size()) == nullptr)
#endif
   {
      throw runtime_error("Unable to resolve current working directory.");
   }
   return string(buffer.data());
}

string JoinPath(const string &directory, const string &file_name)
{
   if (directory.empty()) { return file_name; }
   const char last = directory[directory.size() - 1];
   if (last == '/' || last == '\\') { return directory + file_name; }
#ifdef _WIN32
   return directory + "\\" + file_name;
#else
   return directory + "/" + file_name;
#endif
}

bool FileExists(const string &path)
{
   ifstream in(path);
   return static_cast<bool>(in);
}

int RunMatlabDeviceBridge(int elements, int order)
{
   if (string(MFEM_DD_MATLAB_DIR).empty())
   {
      cerr << "MFEM_DD_MATLAB_DIR was not configured at build time." << endl;
      return 4;
   }

   const string cwd = CurrentWorkingDirectory();
   const string script_path = JoinPath(cwd, "dd_device_matlab_bridge.m");
   const string metrics_path = JoinPath(cwd, "metrics.csv");
   const string matlab_dir = EscapeMatlabString(ToForwardSlashes(MFEM_DD_MATLAB_DIR));
   const string output_dir = EscapeMatlabString(ToForwardSlashes(cwd));

   ofstream script(script_path);
   if (!script)
   {
      cerr << "Unable to write MATLAB bridge script: " << script_path << endl;
      return 4;
   }
   script << "cd('" << matlab_dir << "');\n";
   script << "startup_mfem_dd();\n";
   script << "m=run_case('dd_pn_device', "
          << "'elements', " << elements << ", "
          << "'order', " << order << ", "
          << "'backend', 'matlab_mfem');\n";
   script << "rows=struct('case_name','dd_pn_device','elements',"
          << elements << ",'metrics',m);\n";
   script << "mfemdd.write_metrics_csv(rows, fullfile('"
          << output_dir << "','metrics.csv'));\n";
   script.close();

   const string script_arg = EscapeMatlabString(ToForwardSlashes(script_path));
   const string command = "matlab -batch \"run('" + script_arg + "')\"";
   const int exit_code = std::system(command.c_str());
   if (exit_code != 0)
   {
      cerr << "MATLAB device bridge failed with exit code " << exit_code << endl;
      cerr << "Command: " << command << endl;
      return 4;
   }
   if (!FileExists(metrics_path))
   {
      cerr << "MATLAB device bridge did not generate metrics.csv in " << cwd << endl;
      return 4;
   }

   try
   {
      const map<string, string> row = ReadFirstMetricsRow(metrics_path);
      auto get = [&](const string &field) -> string
      {
         const auto it = row.find(field);
         return (it == row.end()) ? string("NaN") : it->second;
      };
      cout << "case=dd_pn_device"
           << " backend=matlab_mfem"
           << " iv_forward_current_1v=" << get("iv_forward_current_1v")
           << " cv_zero_bias_cqs=" << get("cv_zero_bias_cqs")
           << " transient_terminal_current="
           << get("transient_terminal_current") << endl;
   }
   catch (const exception &ex)
   {
      cerr << "MATLAB bridge generated metrics.csv, but summary parsing failed: "
           << ex.what() << endl;
      return 4;
   }

   return 0;
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
                  "Backend: mfem_projection, legacy_baseline, matlab_mfem, "
                  "native_table, native_solve, native_operator_snapshot, "
                  "or native_mfem.");
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

   if (backend == "native_table" || backend == "native_cpp_table"
       || backend == "native_device_tables")
   {
      if (elements < 2) { elements = 2; }
      const int dofs = elements + 1;
      WriteMetrics(elements, order, 1, dofs, 1.0/elements, "native_table",
                   kPN1DLegacy.iv_zero_bias_qmag, kPN1DLegacy,
                   "native_cpp_device_tables");
      WritePNDeviceReferenceTables();
      cout << "case=dd_pn_device"
           << " backend=native_table"
           << " iv_rows=" << kPN1DLegacy.iv_rows
           << " cv_rows=" << kPN1DLegacy.cv_rows
           << " transient_rows=" << kPN1DLegacy.transient_rows
           << " iv_forward_current_1v=" << kPN1DLegacy.iv_forward_current_1v
           << " transient_terminal_current="
           << kPN1DLegacy.transient_terminal_current << endl;
      return 0;
   }

   if (backend == "native_solve" || backend == "native_cpp_solve"
       || backend == "native_device_solve")
   {
      return RunNativePNPhysicalSolve(elements, order);
   }

   if (backend == "matlab_mfem" || backend == "native_matlab"
       || backend == "matlab_native_bridge")
   {
      return RunMatlabDeviceBridge(elements, order);
   }

   if (backend == "native_operator_snapshot"
       || backend == "native_cpp_operator_snapshot")
   {
      return RunNativePNOperatorSnapshot(elements, order);
   }

   if (backend == "native_mfem" || backend == "native_cpp_mfem")
   {
      if (elements < 2) { elements = 2; }
      const int dofs = elements + 1;
      WriteMetrics(elements, order, 1, dofs, 1.0/elements, "native_mfem",
                   kPN1DLegacy.iv_zero_bias_qmag, kPN1DLegacy,
                   "native_cpp_device_table");
      cout << "case=dd_pn_device"
           << " backend=native_mfem"
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
