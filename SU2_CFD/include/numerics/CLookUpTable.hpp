#pragma once

#include "../../Common/include/option_structure.hpp"
#include "CFileReaderLUT.hpp"
#include "CTrapezoidalMap.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

class CLookUpTable{

protected:

  int rank;  /*!< \brief MPI Rank. */

  string file_name_lut;
  string type_lut;
  string version_lut;
  string version_reader;
  unsigned long n_points;
  unsigned long n_triangles;
  unsigned long n_variables;
  unsigned long n_hull_points;

  su2double limits_table_enth[2];
  su2double limits_table_prog[2];

  /* !brief
   * Holds the variable names stored in the table file. 
   * Order is in sync with data
   */
  vector< string > names_var;

    /* !brief
   * Holds all data stored in the table. First index addresses the variable
   * while second index addresses the point.
   */
  vector< vector< su2double > > table_data;
    
  vector< vector< unsigned long > > triangles;
  vector< vector< unsigned long > > edges;
  vector< vector< unsigned long > > edge_to_triangle;

  vector< unsigned long > hull;

  CTrapezoidalMap trap_map_prog_enth;

  vector< vector< unsigned long > > interp_points;

  vector< vector< vector< su2double > > > interp_mat_inv_prog_enth;

  inline int GetIndexOfVar(string nameVar) {
    int index;
    int endoflist;
    // FIXME dan: there should be an error handling when nameVar is not found in names_var
    index =  (int)(find(names_var.begin(), names_var.end(), nameVar) - names_var.begin());
    endoflist = names_var.size();
    if (index == endoflist){
      index = -1;
      string error_msg = "Variable '";
      error_msg.append(nameVar);
      error_msg.append("' is not in the lookup table.");
      SU2_MPI::Error(error_msg, CURRENT_FUNCTION);
    }
    return index;
  }
  
  inline const vector<su2double> &GetData(string name_var) {
    int ix_var = GetIndexOfVar(name_var);
    return table_data.at(ix_var);
  }
  
  inline const vector< vector< unsigned long > > &GetEdges() const {
    return edges;
  }
  
  inline const vector< vector< unsigned long > > &GetEdgeToTriangle() const {
    return edge_to_triangle;
  }
  
  void FindTableLimits(string name_prog, string name_enth);

  void IdentifyUniqueEdges();

  void LoadTableRaw(string file_name_lut);

  void ComputeInterpCoeffs(string name_prog, string name_enth);
  
  void GetInterpMatInv(const vector<su2double>       &vec_x,
                       const vector<su2double>       &vec_y,
                       vector<unsigned long>         &point_ids,
                       vector< vector< su2double > > &interp_mat_inv);

  void GetInterpCoeffs(su2double                     val_x,
                       su2double                     val_y,
                       vector< vector< su2double > > &interp_mat_inv,
                       vector<su2double>             &interp_coeffs);

  void GaussianInverse(vector<vector<su2double> > &mat,
                       vector<vector<su2double> > &mat_inv);

  su2double Interpolate(const vector<su2double> &val_samples,
                        vector<unsigned long>   &val_triangle,
                        vector<su2double>       &val_interp_coeffs);

  unsigned long FindNearestNeighborOnHull(su2double val_enth,
                                          su2double val_prog, string name_prog,string name_enth);

  bool IsInTriangle(su2double val_x, su2double val_y, unsigned long val_id_triangle, string name_prog,string name_enth);

  inline su2double TriArea(su2double x1, su2double y1, su2double x2, su2double y2, su2double x3, su2double y3) {
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) * 0.5);
  }

 public:
  CLookUpTable(string file_name_lut, string name_prog, string name_enth);

  void PrintTableInfo();

  unsigned long LookUp_ProgEnth(string     val_name_var,
                                su2double *val_var,
                                su2double  val_prog,
                                su2double  val_enth, string name_prog, string name_enth);

  unsigned long LookUp_ProgEnth(vector<string>     &val_names_var,
                                vector<su2double*> &val_vars,
                                su2double           val_prog,
                                su2double           val_enth, string name_prog, string name_enth);

  unsigned long LookUp_ProgEnth(vector<string>     &val_names_var,
                                vector<su2double>  &val_vars,
                                su2double           val_prog,
                                su2double           val_enth, string name_prog, string name_enth);

  inline pair<su2double, su2double> GetTableLimitsEnth() {
    return make_pair(limits_table_enth[0], limits_table_enth[1]);
}

inline pair<su2double, su2double> GetTableLimitsProg() {
  return make_pair(limits_table_prog[0], limits_table_prog[1]);
}

};
