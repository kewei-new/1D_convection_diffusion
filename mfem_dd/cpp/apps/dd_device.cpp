#include "mfem.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
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

   if (backend == "matlab_mfem" || backend == "native_matlab"
       || backend == "matlab_native_bridge")
   {
      return RunMatlabDeviceBridge(elements, order);
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
