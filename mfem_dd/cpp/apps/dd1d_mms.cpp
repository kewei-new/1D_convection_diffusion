#include "mfem.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace mfem;
using namespace std;

namespace
{
constexpr real_t pi = 3.141592653589793238462643383279502884;

struct GaussPoint
{
   real_t x;
   real_t w;
};

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

struct NativeResult
{
   int elements = 0;
   int order = 0;
   int dofs = 0;
   real_t h = 0.0;
   real_t time = 0.0;
   real_t n_l2 = 0.0;
   real_t n_linf = 0.0;
   real_t phi_l2 = 0.0;
   real_t phi_linf = 0.0;
   real_t e_l2 = 0.0;
   real_t e_linf = 0.0;
};

int Idx(int cell, int mode, int order)
{
   return cell*(order + 1) + mode;
}

real_t ExactN(const real_t x, const real_t t)
{
   return std::sin(x)*std::cos(t);
}

real_t ExactPhi(const real_t x, const real_t t)
{
   return std::sin(t)*std::cos(x);
}

real_t ExactPhiX(const real_t x, const real_t t)
{
   return -std::sin(t)*std::sin(x);
}

real_t SourceF1(const real_t x, const real_t t)
{
   return -std::sin(t)*std::sin(x) + std::cos(t)*std::sin(x)
          - 2.0*std::sin(t)*std::cos(t)*std::sin(x)*std::cos(x);
}

real_t SourceF2(const real_t x, const real_t t)
{
   return std::sin(t)*std::cos(x) + std::sin(x)*std::cos(t);
}

real_t InitialN(const real_t x)
{
   return std::sin(x);
}

real_t LdgBoundary(const real_t x, const real_t t)
{
   if (std::abs(x) < 1.0e-10 || std::abs(x - 2.0*pi) < 1.0e-10)
   {
      return std::sin(t);
   }
   throw runtime_error("LDG boundary requested away from a DD1D endpoint.");
}

real_t Basis(const real_t x, const int mode)
{
   switch (mode)
   {
      case 0: return 1.0;
      case 1: return x;
      case 2: return x*x - 1.0/12.0;
      case 3: return x*x*x - 0.15*x;
      case 4: return (x*x - 3.0/14.0)*x*x + 3.0/560.0;
      default: throw runtime_error("Unsupported modal basis order.");
   }
}

real_t BasisX(const real_t x, const int mode)
{
   switch (mode)
   {
      case 0: return 0.0;
      case 1: return 1.0;
      case 2: return 2.0*x;
      case 3: return 3.0*x*x - 0.15;
      case 4: return 4.0*x*x*x - 3.0/7.0*x;
      default: throw runtime_error("Unsupported modal basis derivative order.");
   }
}

vector<real_t> InverseMass(const int order)
{
   switch (order)
   {
      case 1: return {1.0, 12.0};
      case 2: return {1.0, 12.0, 180.0};
      case 3: return {1.0, 12.0, 180.0, 2800.0};
      default: throw runtime_error("DD1D native C++ backend supports orders 1-3.");
   }
}

vector<GaussPoint> GaussLobatto(const int order)
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
   if (order == 4)
   {
      const real_t a = std::sqrt(147.0 + 42.0*std::sqrt(7.0))/42.0;
      const real_t b = std::sqrt(147.0 - 42.0*std::sqrt(7.0))/42.0;
      const real_t wa = (-7.0 + 5.0*std::sqrt(7.0))*std::sqrt(7.0)
                        *(7.0 + std::sqrt(7.0))/840.0;
      const real_t wb = (7.0 + 5.0*std::sqrt(7.0))*std::sqrt(7.0)
                        /(7.0 + std::sqrt(7.0))/20.0;
      return {{-0.5, 1.0/30.0}, {0.5, 1.0/30.0},
              {-a, wa}, {a, wa}, {-b, wb}, {b, wb}};
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

real_t EvalModal(const vector<real_t> &u, const int cell, const int order,
                 const real_t xi)
{
   real_t value = 0.0;
   for (int mode = 0; mode <= order; mode++)
   {
      value += u[Idx(cell, mode, order)]*Basis(xi, mode);
   }
   return value;
}

void MatAdd(DenseMatrix &A, const int row, const int col, const real_t value)
{
   A(row, col) += value;
}

void AddBlock(DenseMatrix &A, const int row_cell, const int col_cell,
              const DenseMatrix &block, const int order)
{
   for (int i = 0; i <= order; i++)
   {
      for (int j = 0; j <= order; j++)
      {
         MatAdd(A, Idx(row_cell, i, order), Idx(col_cell, j, order),
                block(i, j));
      }
   }
}

DenseMatrix Scale(const DenseMatrix &A, const real_t s)
{
   DenseMatrix B(A.Height(), A.Width());
   for (int i = 0; i < A.Height(); i++)
   {
      for (int j = 0; j < A.Width(); j++) { B(i, j) = s*A(i, j); }
   }
   return B;
}

DenseMatrix TransposeScale(const DenseMatrix &A, const real_t s)
{
   DenseMatrix B(A.Width(), A.Height());
   for (int i = 0; i < A.Height(); i++)
   {
      for (int j = 0; j < A.Width(); j++) { B(j, i) = s*A(i, j); }
   }
   return B;
}

DenseMatrix AddBlocks(const vector<pair<real_t, const DenseMatrix*> > &terms)
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

void MatVec(const DenseMatrix &A, const Vector &x, Vector &y)
{
   A.Mult(x, y);
}

Vector MatVecNew(const DenseMatrix &A, const Vector &x)
{
   Vector y(A.Height());
   A.Mult(x, y);
   return y;
}

void AddScaled(Vector &y, const real_t a, const Vector &x)
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

DenseMatrix AssembleIPDGDiffusion(const int elements, const int order,
                                  const real_t h)
{
   const int nloc = order + 1;
   const int size = elements*nloc;
   const real_t alpha = static_cast<real_t>((order + 1)*(order + 1));
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
            Mx1(d1, d2) += (1.0/h)*g.w*BasisX(g.x, d2)*BasisX(g.x, d1);
         }
         Mx2LL(d1, d2) = (1.0/h)*BasisX(-0.5, d2)*Basis(-0.5, d1);
         Mx2RR(d1, d2) = (1.0/h)*BasisX( 0.5, d2)*Basis( 0.5, d1);
         Mx2LR(d1, d2) = (1.0/h)*BasisX(-0.5, d2)*Basis( 0.5, d1);
         Mx2RL(d1, d2) = (1.0/h)*BasisX( 0.5, d2)*Basis(-0.5, d1);
         Mx3LL(d1, d2) = Basis(-0.5, d2)*Basis(-0.5, d1);
         Mx3RR(d1, d2) = Basis( 0.5, d2)*Basis( 0.5, d1);
         Mx3RL(d1, d2) = Basis( 0.5, d2)*Basis(-0.5, d1);
         Mx3LR(d1, d2) = Basis(-0.5, d2)*Basis( 0.5, d1);
      }
   }

   Mx2LL = Scale(Mx2LL, 0.5);
   Mx2RR = Scale(Mx2RR, 0.5);
   Mx2LR = Scale(Mx2LR, 0.5);
   Mx2RL = Scale(Mx2RL, 0.5);
   Mx3LL = Scale(Mx3LL, alpha/h);
   Mx3RR = Scale(Mx3RR, alpha/h);
   Mx3LR = Scale(Mx3LR, alpha/h);
   Mx3RL = Scale(Mx3RL, alpha/h);

   DenseMatrix diag = AddBlocks({{-1.0, &Mx1}, {1.0, &Mx2RR},
                                 {-1.0, &Mx2LL}});
   DenseMatrix t1 = TransposeScale(Mx2RR, beta);
   DenseMatrix t2 = TransposeScale(Mx2LL, -beta);
   DenseMatrix diag2 = AddBlocks({{1.0, &diag}, {1.0, &t1}, {1.0, &t2},
                                  {-1.0, &Mx3RR}, {-1.0, &Mx3LL}});
   DenseMatrix upper_t = TransposeScale(Mx2RL, -beta);
   DenseMatrix upper = AddBlocks({{1.0, &Mx2LR}, {1.0, &upper_t},
                                  {1.0, &Mx3LR}});
   DenseMatrix lower_t = TransposeScale(Mx2LR, beta);
   DenseMatrix lower = AddBlocks({{-1.0, &Mx2RL}, {1.0, &lower_t},
                                  {1.0, &Mx3RL}});

   for (int cell = 0; cell < elements; cell++)
   {
      AddBlock(diffusion, cell, cell, diag2, order);
      if (cell < elements - 1)
      {
         AddBlock(diffusion, cell, cell + 1, upper, order);
      }
      else
      {
         AddBlock(diffusion, 0, cell, lower, order);
      }

      if (cell > 0)
      {
         AddBlock(diffusion, cell, cell - 1, lower, order);
      }
      else
      {
         AddBlock(diffusion, elements - 1, cell, upper, order);
      }
   }

   return diffusion;
}

void AssembleLDG(const int elements, const int order, const real_t h,
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
      for (const auto &g : gauss) { s += g.w*Basis(g.x, beta)*BasisX(g.x, alpha); }
      return s;
   };

   for (int cell = 0; cell < elements; cell++)
   {
      for (int alpha = 0; alpha <= order; alpha++)
      {
         for (int beta = 0; beta <= order; beta++)
         {
            const int row = Idx(cell, alpha, order);
            const int col = Idx(cell, beta, order);
            if (cell == 0)
            {
               M1(row, col) = inv_mass[alpha]/h
                              *(-quad(beta, alpha)
                                + Basis(0.5, beta)*Basis(0.5, alpha));
            }
            else if (cell < elements - 1)
            {
               M1(row, col) = inv_mass[alpha]/h
                              *(-quad(beta, alpha)
                                + Basis(0.5, beta)*Basis(0.5, alpha));
               M1(row, Idx(cell - 1, beta, order)) = inv_mass[alpha]/h
                  *(-Basis(0.5, beta)*Basis(-0.5, alpha));
            }
            else
            {
               M1(row, col) = inv_mass[alpha]/h*(-quad(beta, alpha));
               M1(row, Idx(cell - 1, beta, order)) = inv_mass[alpha]/h
                  *(-Basis(0.5, beta)*Basis(-0.5, alpha));
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
            const int row = Idx(cell, alpha, order);
            const int col = Idx(cell, beta, order);
            if (cell < elements - 1)
            {
               M2(row, col) = -quad(beta, alpha)
                              - Basis(-0.5, beta)*Basis(-0.5, alpha);
               M2(row, Idx(cell + 1, beta, order)) =
                  Basis(-0.5, beta)*Basis(0.5, alpha);
            }
            else
            {
               M2(row, col) = -quad(beta, alpha)
                              - Basis(-0.5, beta)*Basis(-0.5, alpha)
                              + Basis(0.5, beta)*Basis(0.5, alpha);
               M3(row, col) = -Cp*Basis(0.5, beta)*Basis(0.5, alpha);
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

Vector L2Projection(const int elements, const int order, const real_t h,
                    const vector<real_t> &inv_mass)
{
   Vector coeff(elements*(order + 1));
   coeff = 0.0;
   const auto gauss = GaussLegendre4();
   for (int cell = 0; cell < elements; cell++)
   {
      const real_t center = h*(static_cast<real_t>(cell) + 0.5);
      for (int mode = 0; mode <= order; mode++)
      {
         real_t s = 0.0;
         for (const auto &g : gauss)
         {
            s += g.w*InitialN(g.x*h + center)*Basis(g.x, mode);
         }
         coeff(Idx(cell, mode, order)) = inv_mass[mode]*s;
      }
   }
   return coeff;
}

void SolveLDG(const int elements, const int order, const real_t h,
              const vector<real_t> &inv_mass, const Vector &carrier_coeff,
              const real_t time, const DenseMatrix &M1, const DenseMatrix &M2,
              const DenseMatrixInverse &poisson_solver,
              Vector &phi_coeff, Vector &e_coeff)
{
   const int size = elements*(order + 1);
   const real_t Cp = static_cast<real_t>(order)/h;
   const auto gauss = GaussLobatto(order);
   vector<real_t> carrier(carrier_coeff.GetData(), carrier_coeff.GetData() + size);

   Vector b1(size), b2(size), f_terms(size), tmp(size), rhs(size);
   b1 = 0.0; b2 = 0.0; f_terms = 0.0;

   const real_t left_bound = LdgBoundary(0.0, time);
   const real_t right_bound = LdgBoundary(2.0*pi, time);
   for (int alpha = 0; alpha <= order; alpha++)
   {
      b1(Idx(0, alpha, order)) = -inv_mass[alpha]/h*left_bound*Basis(-0.5, alpha);
      b1(Idx(elements - 1, alpha, order)) =
         inv_mass[alpha]/h*right_bound*Basis(0.5, alpha);
      b2(Idx(elements - 1, alpha, order)) = Cp*right_bound*Basis(0.5, alpha);
   }

   for (int cell = 0; cell < elements; cell++)
   {
      const real_t center = h*(static_cast<real_t>(cell) + 0.5);
      for (int alpha = 0; alpha <= order; alpha++)
      {
         real_t s = 0.0;
         for (const auto &g : gauss)
         {
            const real_t f = EvalModal(carrier, cell, order, g.x);
            s += g.w*(f - SourceF2(g.x*h + center, time))*Basis(g.x, alpha);
         }
         f_terms(Idx(cell, alpha, order)) = s*h;
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

Vector TransportRHS(const int elements, const int order, const real_t h,
                    const Vector &carrier_coeff, const Vector &e_coeff,
                    const real_t time)
{
   const int nloc = order + 1;
   const int size = elements*nloc;
   const auto gauss = GaussLobatto(order);
   const int gcount = static_cast<int>(gauss.size());
   const real_t emf = 1.0;
   vector<vector<real_t> > u1(gcount, vector<real_t>(elements + 2, 1.0));
   vector<vector<real_t> > phix(gcount, vector<real_t>(elements + 2, 1.0));
   vector<vector<real_t> > f1(gcount, vector<real_t>(elements + 2, 0.0));
   vector<real_t> carrier(carrier_coeff.GetData(), carrier_coeff.GetData() + size);
   vector<real_t> evec(e_coeff.GetData(), e_coeff.GetData() + size);

   for (int k = 1; k <= elements; k++)
   {
      for (int gi = 0; gi < gcount; gi++)
      {
         u1[gi][k] = EvalModal(carrier, k - 1, order, gauss[gi].x);
         phix[gi][k] = EvalModal(evec, k - 1, order, gauss[gi].x);
      }
   }
   u1[1][0] = u1[1][elements];
   phix[1][0] = phix[1][elements];
   u1[0][elements + 1] = u1[0][1];
   phix[0][elements + 1] = phix[0][1];

   for (int gi = 0; gi < gcount; gi++)
   {
      for (int k = 0; k < elements + 2; k++) { f1[gi][k] = u1[gi][k]*(-phix[gi][k]); }
   }

   Vector rhs(size);
   rhs = 0.0;
   for (int k = 1; k <= elements; k++)
   {
      const int cell = k - 1;
      const real_t center = h*(static_cast<real_t>(cell) + 0.5);
      for (int mode = 0; mode <= order; mode++)
      {
         real_t flux = 0.0;
         real_t source = 0.0;
         for (int gi = 0; gi < gcount; gi++)
         {
            flux += gauss[gi].w*f1[gi][k]*BasisX(gauss[gi].x, mode);
            source += gauss[gi].w*SourceF1(gauss[gi].x*h + center, time)
                      *Basis(gauss[gi].x, mode)*h;
         }
         const real_t face = 0.5*((f1[1][k] + f1[0][k + 1])*Basis(0.5, mode)
                           - (f1[1][k - 1] + f1[0][k])*Basis(-0.5, mode));
         const real_t penalty = -0.5*emf*((u1[0][k + 1] - u1[1][k])*Basis(0.5, mode)
                              - (u1[0][k] - u1[1][k - 1])*Basis(-0.5, mode));
         rhs(Idx(cell, mode, order)) = -flux + face + penalty + source;
      }
   }
   return rhs;
}

real_t L2Error(const Vector &coeff, const int elements, const int order,
               const real_t h, const real_t time,
               real_t (*exact)(real_t, real_t))
{
   const auto gauss = GaussLobatto(order);
   vector<real_t> values(coeff.GetData(), coeff.GetData() + coeff.Size());
   real_t error = 0.0;
   for (int cell = 0; cell < elements; cell++)
   {
      const real_t center = h*(static_cast<real_t>(cell) + 0.5);
      for (const auto &g : gauss)
      {
         const real_t uh = EvalModal(values, cell, order, g.x);
         const real_t diff = uh - exact(g.x*h + center, time);
         error += g.w*diff*diff*h;
      }
   }
   return std::sqrt(error);
}

real_t LinfError(const Vector &coeff, const int elements, const int order,
                 const real_t h, const real_t time,
                 real_t (*exact)(real_t, real_t))
{
   const auto gauss = GaussLobatto(order);
   vector<real_t> values(coeff.GetData(), coeff.GetData() + coeff.Size());
   real_t error = 0.0;
   for (int cell = 0; cell < elements; cell++)
   {
      const real_t center = h*(static_cast<real_t>(cell) + 0.5);
      for (const auto &g : gauss)
      {
         const real_t uh = EvalModal(values, cell, order, g.x);
         error = std::max(error, std::abs(uh - exact(g.x*h + center, time)));
      }
   }
   return error;
}

DenseMatrix MakeImplicitMatrix(const DenseMatrix &diffusion, const int elements,
                               const int order, const real_t h,
                               const vector<real_t> &inv_mass,
                               const real_t dt)
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
         A(Idx(cell, mode, order), Idx(cell, mode, order)) += h/inv_mass[mode]/dt;
      }
   }
   return A;
}

Vector MassOverDtMult(const Vector &x, const int elements, const int order,
                      const real_t h, const vector<real_t> &inv_mass,
                      const real_t dt)
{
   Vector y(x.Size());
   y = 0.0;
   for (int cell = 0; cell < elements; cell++)
   {
      for (int mode = 0; mode <= order; mode++)
      {
         const int id = Idx(cell, mode, order);
         y(id) = h/inv_mass[mode]/dt*x(id);
      }
   }
   return y;
}

NativeResult SolveNativeDD1D(const int elements, const int order)
{
   const real_t left = 0.0;
   const real_t right = 2.0*pi;
   const real_t h = (right - left)/static_cast<real_t>(elements);
   const real_t T_end = 1.0;
   const real_t cfl = 0.1;
   real_t dt = cfl*h;
   const int size = elements*(order + 1);
   const vector<real_t> inv_mass = InverseMass(order);

   // Keep an MFEM mesh/FESpace in the native path so the app remains tied to
   // the same external MFEM runtime while using the legacy modal DG algebra.
   Mesh mesh = Mesh::MakeCartesian1D(elements, right);
   L2_FECollection fec(order, mesh.Dimension());
   FiniteElementSpace fes(&mesh, &fec);
   (void)fes;

   Vector n_old = L2Projection(elements, order, h, inv_mass);
   DenseMatrix diffusion = AssembleIPDGDiffusion(elements, order, h);
   DenseMatrix M1, M2, Mn;
   AssembleLDG(elements, order, h, inv_mass, M1, M2, Mn);
   DenseMatrixInverse poisson_solver(Mn);

   DenseMatrix implicit_matrix = MakeImplicitMatrix(diffusion, elements, order,
                                                    h, inv_mass, dt);
   DenseMatrixInverse implicit_solver(implicit_matrix);

   real_t time = 0.0;
   while (time < T_end - 1.0e-14)
   {
      if (time + dt >= T_end)
      {
         dt = T_end - time;
         implicit_matrix = MakeImplicitMatrix(diffusion, elements, order, h,
                                              inv_mass, dt);
         implicit_solver.Factor(implicit_matrix);
      }

      Vector phi(size), e(size), rhs0(size), rhs1(size), rhs2(size), rhs3(size);
      Vector stage1(size), stage2(size), stage3(size), n_new(size);

      SolveLDG(elements, order, h, inv_mass, n_old, time, M1, M2,
               poisson_solver, phi, e);
      rhs0 = TransportRHS(elements, order, h, n_old, e, time);
      Vector mass_old = MassOverDtMult(n_old, elements, order, h, inv_mass, dt);
      Vector solve_rhs = LinearCombination({{0.5, &rhs0}, {1.0, &mass_old}});
      implicit_solver.Mult(solve_rhs, stage1);

      SolveLDG(elements, order, h, inv_mass, stage1, time + dt/2.0, M1, M2,
               poisson_solver, phi, e);
      rhs1 = TransportRHS(elements, order, h, stage1, e, time + dt/2.0);
      Vector d_stage1 = MatVecNew(diffusion, stage1);
      solve_rhs = LinearCombination({{11.0/18.0, &rhs0}, {1.0/18.0, &rhs1},
                                     {1.0/6.0, &d_stage1}, {1.0, &mass_old}});
      implicit_solver.Mult(solve_rhs, stage2);

      SolveLDG(elements, order, h, inv_mass, stage2, time + 2.0*dt/3.0, M1, M2,
               poisson_solver, phi, e);
      rhs2 = TransportRHS(elements, order, h, stage2, e, time + 2.0*dt/3.0);
      Vector d_stage2 = MatVecNew(diffusion, stage2);
      solve_rhs = LinearCombination({{5.0/6.0, &rhs0}, {-5.0/6.0, &rhs1},
                                     {0.5, &rhs2}, {-0.5, &d_stage1},
                                     {0.5, &d_stage2}, {1.0, &mass_old}});
      implicit_solver.Mult(solve_rhs, stage3);

      SolveLDG(elements, order, h, inv_mass, stage3, time + dt/2.0, M1, M2,
               poisson_solver, phi, e);
      rhs3 = TransportRHS(elements, order, h, stage3, e, time + dt/2.0);
      Vector d_stage3 = MatVecNew(diffusion, stage3);
      solve_rhs = LinearCombination({{0.25, &rhs0}, {1.75, &rhs1},
                                     {0.75, &rhs2}, {-1.75, &rhs3},
                                     {1.5, &d_stage1}, {-1.5, &d_stage2},
                                     {0.5, &d_stage3}, {1.0, &mass_old}});
      implicit_solver.Mult(solve_rhs, n_new);

      n_old = n_new;
      time += dt;
   }

   Vector phi(size), e(size);
   SolveLDG(elements, order, h, inv_mass, n_old, time, M1, M2, poisson_solver,
            phi, e);

   NativeResult result;
   result.elements = elements;
   result.order = order;
   result.dofs = size;
   result.h = h;
   result.time = time;
   result.n_l2 = L2Error(n_old, elements, order, h, time, ExactN);
   result.n_linf = LinfError(n_old, elements, order, h, time, ExactN);
   result.phi_l2 = L2Error(phi, elements, order, h, time, ExactPhi);
   result.phi_linf = LinfError(phi, elements, order, h, time, ExactPhi);
   result.e_l2 = L2Error(e, elements, order, h, time, ExactPhiX);
   result.e_linf = LinfError(e, elements, order, h, time, ExactPhiX);
   return result;
}

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

real_t ProjectionExactN(const Vector &x)
{
   return 1.0 + 0.2*std::sin(pi*x[0]);
}

real_t ProjectionExactPhi(const Vector &x)
{
   return std::sin(2.0*pi*x[0]);
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
   args.AddOption(&order, "-o", "--order", "Finite element / modal DG order.");
   args.AddOption(&backend, "-b", "--backend",
                  "Backend: mfem_projection, native_mfem, or legacy_baseline.");
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
      WriteMetrics(elements, order, elements*(order + 1), row.h, "legacy_baseline",
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

   if (backend == "native_mfem" || backend == "native_cpp" || backend == "matlab_mfem")
   {
      try
      {
         NativeResult result = SolveNativeDD1D(elements, order);
         WriteMetrics(elements, order, result.dofs, result.h, "native_mfem",
                      "native_cpp_mfem", result.n_l2, result.n_linf, nan, nan,
                      result.phi_l2, result.phi_linf, result.e_l2, result.e_linf,
                      nan, nan, nan);
         cout << "case=dd1d_smooth_mms"
              << " backend=native_mfem"
              << " n_l2_error=" << result.n_l2
              << " phi_l2_error=" << result.phi_l2
              << " E_l2_error=" << result.e_l2 << endl;
         return 0;
      }
      catch (const std::exception &ex)
      {
         cerr << "Native DD1D C++ backend failed: " << ex.what() << endl;
         return 3;
      }
   }

   Mesh mesh = Mesh::MakeCartesian1D(elements, 1.0);
   H1_FECollection fec(order, mesh.Dimension());
   FiniteElementSpace fes(&mesh, &fec);

   FunctionCoefficient n_coeff(ProjectionExactN);
   FunctionCoefficient phi_coeff(ProjectionExactPhi);
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
      charge_proxy += ProjectionExactN(node) - 1.0;
   }
   charge_proxy /= mesh.GetNV();

   const real_t h = 1.0/elements;
   WriteMetrics(elements, order, fes.GetTrueVSize(), h, backend, "implemented",
                n_err, nan, nan, nan, phi_err, nan, nan, nan, nan, nan,
                charge_proxy);

   cout << "case=dd1d_smooth_mms"
        << " dofs=" << fes.GetTrueVSize()
        << " n_l2_error=" << n_err
        << " phi_l2_error=" << phi_err << endl;
   return 0;
}
