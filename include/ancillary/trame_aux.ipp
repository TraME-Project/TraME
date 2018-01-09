/*################################################################################
  ##
  ##   Copyright (C) 2015-2018 the TraME Team:
  ##      Alfred Galichon
  ##      Keith O'Hara
  ##
  ##   This file is part of TraME.
  ##
  ##   TraME is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   TraME is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

/*
 * auxiliary functions
 *
 * Keith O'Hara
 * 08/08/2016
 *
 * This version:
 * 01/09/2018
 */

//
// LogSumExp

inline
double
lse(const arma::mat& X) {
    double max_elem = X.max();

    return max_elem + std::log( arma::accu( arma::exp(X - max_elem) ) );
}

//
// Find indices that correspond to maximum values

inline
arma::uvec
which_max(const arma::mat& X, const int which_dim)
{
    const int n = X.n_rows;
    const int k = X.n_cols;
    //
    int max_ind = 0;
    double max_val = 0;
    arma::uvec max_vec;

    if (which_dim==0) { // each column
        max_vec.set_size(k);

        for (int j=0; j < k; j++) {
            max_val = X(0,j);
            max_ind = 0;

            for (int i=1; i < n; i++) {
                if (X(i,j) > max_val) {
                    max_val = X(i,j);
                    max_ind = i;
                }
            }

            max_vec(j) = max_ind;
        }
    } else { // each row
        max_vec.set_size(n);

        for (int i=0; i<n; i++) {
            max_val = X(i,0);
            max_ind = 0;

            for (int j=0; j < k; j++) {
                if (X(i,j) > max_val) {
                    max_val = X(i,j);
                    max_ind = j;
                }
            }

            max_vec(i) = max_ind;
        }
    }
    //
    return max_vec;
}

//
// unit vector

inline
arma::vec
unit_vec(const int j, const int n)
{
    arma::vec ret = arma::zeros(n,1);
    ret(j) = 1;

    return ret;
}

//
// Generate an equi-spaced grid of integers

inline
arma::uvec
uvec_linspace(const int a, const int b)
{
    const int n_points = b - a + 1;

    arma::uvec ret(n_points);
    //
    for (int i=0; i < n_points; i++) {
        ret(i) = a + i;
    }
    //
    return ret;
}

//
// convert uword object (const long long int*) to int*
// note: a simple reinterpret_cast will NOT work here

inline
int*
uword_to_int(const arma::uword* var_inp, const int n_elem)
{
    int* var_out = new int[n_elem];

    for (int jj=0; jj < n_elem; jj++) {
        var_out[jj] = var_inp[jj];
    }

    return var_out;
}

//
// by-row reconstruction of a matrix (mimics R's 'byrow=TRUE' option when constructing a matrix)

inline
arma::mat
byrow(const arma::mat& mat_inp, const int n_rows, const int n_cols)
{
    arma::vec vec_tmp = arma::vectorise(mat_inp);
    const int n_vals = vec_tmp.n_elem;

    arma::mat ret(n_rows,n_cols);

    if (n_vals != n_rows*n_cols) { // use repmat in this case
        if (n_vals == n_rows) {
            ret = arma::repmat(vec_tmp,1,n_cols);
        } else if (n_vals == n_cols) {
            ret = arma::repmat(vec_tmp.t(),n_rows,1);
        } else {
            printf("error in byrow\n");
            return ret;
        }
    } else { // otherwise rebuild the input matrix using byrow

        int kr = 0, kc = 0;

        for (int i=0; i < n_rows*n_cols; i++) {
            ret(kr,kc) = vec_tmp(i);

            kc++;

            if (kc >= n_cols) {
                kc = 0;
                kr++;
            }
        }
    }

    return ret;
}

//
// element-by-element addition;
// written to mimic R's approach to element-by-element addition in vector-matrix cases

inline
arma::mat
elem_add(const arma::mat& mat_1, const arma::mat& mat_2)
{
    const int rows_1 = mat_1.n_rows;
    const int rows_2 = mat_2.n_rows;

    const int cols_1 = mat_1.n_cols;
    const int cols_2 = mat_2.n_cols;

    const bool same_rows = (rows_1==rows_2) ? true : false;
    const bool same_cols = (cols_1==cols_2) ? true : false;
    //
    arma::mat ret;

    if (!same_rows && !same_cols) {
        printf("elem_add: need matrices to agree on at least one dimension\n");
        return ret;
    }

    //
    // instead of loops we could also use arma::repmat, not sure which is quicker...

    if (same_rows && same_cols) {
        ret = mat_1 + mat_2;
    } else if (same_rows && !same_cols) { // same #rows, different #cols
        if (cols_1==1) {
            ret.set_size(rows_1,cols_2);
            for (int i=0; i < cols_2; i++) {
                ret.col(i) = mat_1 + mat_2.col(i);
            }
        } else if (cols_2==1) {
            ret.set_size(rows_1,cols_1);
            for (int i=0; i < cols_1; i++) {
                ret.col(i) = mat_1.col(i) + mat_2;
            }
        } else {
            printf("elem_add: need one of the matrices to be a column vector\n");
            return ret;
        }
    } else if (!same_rows && same_cols) { // same #cols, different #rows
        if (rows_1==1) {
            ret.set_size(rows_2,cols_1);
            for (int i=0; i < rows_2; i++) {
                ret.row(i) = mat_1 + mat_2.row(i);
            }
        } else if (rows_2==1) {
            ret.set_size(rows_1,cols_1);
            for (int i=0; i < rows_1; i++) {
                ret.row(i) = mat_1.row(i) + mat_2;
            }
        } else {
            printf("elem_add: need one of the matrices to be a row vector\n");
            return ret;
        }
    } else {
        printf("elem_add: unknown error\n");
        return ret;
    }

    //

    return ret;
}

//
// element-by-element subtraction;
// written to mimic R's approach to element-by-element subtraction in vector-matrix cases

inline
arma::mat
elem_sub(const arma::mat& mat_1, const arma::mat& mat_2)
{
    const int rows_1 = mat_1.n_rows;
    const int rows_2 = mat_2.n_rows;

    const int cols_1 = mat_1.n_cols;
    const int cols_2 = mat_2.n_cols;

    const bool same_rows = (rows_1==rows_2) ? true : false;
    const bool same_cols = (cols_1==cols_2) ? true : false;
    //
    arma::mat ret;

    if (!same_rows && !same_cols) {
        printf("elem_sub: need matrices to agree on at least one dimension\n");
        return ret;
    }

    //
    // instead of loops we could also use arma::repmat, not sure which is quicker...

    if (same_rows && same_cols) {
        ret = mat_1 - mat_2;
    } else if (same_rows && !same_cols) {
        if (cols_1==1) {
            ret.set_size(rows_1,cols_2);
            for (int i=0; i < cols_2; i++) {
                ret.col(i) = mat_1 - mat_2.col(i);
            }
        } else if (cols_2==1) {
            ret.set_size(rows_1,cols_1);
            for (int i=0; i < cols_1; i++) {
                ret.col(i) = mat_1.col(i) - mat_2;
            }
        } else {
            printf("elem_sub: need one of the matrices to be a column vector\n");
            return ret;
        }
    } else if (!same_rows && same_cols) {
        if (rows_1==1) {
            ret.set_size(rows_2,cols_1);
            for (int i=0; i < rows_2; i++) {
                ret.row(i) = mat_1 - mat_2.row(i);
            }
        } else if (rows_2==1) {
            ret.set_size(rows_1,cols_1);
            for (int i=0; i < rows_1; i++) {
                ret.row(i) = mat_1.row(i) - mat_2;
            }
        } else {
            printf("elem_sub: need one of the matrices to be a row vector\n");
            return ret;
        }
    } else {
        printf("elem_sub: unknown error\n");
        return ret;
    }
    //
    return ret;
}

//
// Hadamard product;
// written to mimic R's approach to element-by-element multiplication in vector-matrix cases

inline
arma::mat
elem_prod(const arma::mat& mat_1, const arma::mat& mat_2)
{
    const int rows_1 = mat_1.n_rows;
    const int rows_2 = mat_2.n_rows;

    const int cols_1 = mat_1.n_cols;
    const int cols_2 = mat_2.n_cols;

    const bool same_rows = (rows_1==rows_2) ? true : false;
    const bool same_cols = (cols_1==cols_2) ? true : false;
    //
    arma::mat ret;

    if (!same_rows && !same_cols) {
        printf("elem_prod: need matrices to agree on at least one dimension\n");
        return ret;
    }

    //
    // instead of loops we could also use arma::repmat, not sure which is quicker...

    if (same_rows && same_cols) {
        ret = mat_1 % mat_2;
    } else if (same_rows && !same_cols) {
        if (cols_1==1) {
            ret.set_size(rows_1,cols_2);
            for (int i=0; i < cols_2; i++) {
                ret.col(i) = mat_1 % mat_2.col(i);
            }
        } else if (cols_2==1) {
            ret.set_size(rows_1,cols_1);
            for (int i=0; i < cols_1; i++) {
                ret.col(i) = mat_1.col(i) % mat_2;
            }
        } else {
            printf("elem_prod: need one of the matrices to be a column vector\n");
            return ret;
        }
    } else if (!same_rows && same_cols) {
        if (rows_1==1) {
            ret.set_size(rows_2,cols_1);
            for (int i=0; i < rows_2; i++) {
                ret.row(i) = mat_1 % mat_2.row(i);
            }
        } else if (rows_2==1) {
            ret.set_size(rows_1,cols_1);
            for (int i=0; i < rows_1; i++) {
                ret.row(i) = mat_1.row(i) % mat_2;
            }
        } else {
            printf("elem_prod: need one of the matrices to be a row vector\n");
            return ret;
        }
    } else {
        printf("elem_prod: unknown error\n");
        return ret;
    }
    //
    return ret;
}

//
// element-by-element division;
// written to mimic R's approach to element-by-element division in vector-matrix cases

inline
arma::mat
elem_div(const arma::mat& mat_1, const arma::mat& mat_2)
{
    const int rows_1 = mat_1.n_rows;
    const int rows_2 = mat_2.n_rows;

    const int cols_1 = mat_1.n_cols;
    const int cols_2 = mat_2.n_cols;

    const bool same_rows = (rows_1==rows_2) ? true : false;
    const bool same_cols = (cols_1==cols_2) ? true : false;
    //
    arma::mat ret;

    if (!same_rows && !same_cols) {
        printf("elem_div: need matrices to agree on at least one dimension\n");
        return ret;
    }

    //
    // instead of loops we could also use arma::repmat, not sure which is quicker...

    if (same_rows && same_cols) {
        ret = mat_1 / mat_2;
    } else if (same_rows && !same_cols) {
        if (cols_1==1) {
            ret.set_size(rows_1,cols_2);
            for (int i=0; i < cols_2; i++) {
                ret.col(i) = mat_1 / mat_2.col(i);
            }
        } else if (cols_2==1) {
            ret.set_size(rows_1,cols_1);
            for (int i=0; i < cols_1; i++) {
                ret.col(i) = mat_1.col(i) / mat_2;
            }
        } else {
            printf("elem_div: need one of the matrices to be a column vector\n");
            return ret;
        }
    } else if (!same_rows && same_cols) {
        if (rows_1==1) {
            ret.set_size(rows_2,cols_1);
            for (int i=0; i < rows_2; i++) {
                ret.row(i) = mat_1 / mat_2.row(i);
            }
        } else if (rows_2==1) {
            ret.set_size(rows_1,cols_1);
            for (int i=0; i < rows_1; i++) {
                ret.row(i) = mat_1.row(i) / mat_2;
            }
        } else {
            printf("elem_div: need one of the matrices to be a row vector\n");
            return ret;
        }
    } else {
        printf("elem_div: unknown error\n");
        return ret;
    }
    //
    return ret;
}

//
// return the minimal-value in the matrix

inline
double
elem_min(const arma::mat& mat_1)
{
    return mat_1.min();
}

//
// pointwise minimum between two matrices

inline
arma::mat
elem_min(const arma::mat& mat_1, const arma::mat& mat_2)
{
    return arma::min(mat_1,mat_2);
}

//
// pointwise minimum between a matrix and a single value

inline
arma::mat
elem_min(const arma::mat& mat_1, const double comp_val)
{
    const int rows_1 = mat_1.n_rows;
    const int cols_1 = mat_1.n_cols;
    //
    arma::mat comp_mat(rows_1,cols_1);
    comp_mat.fill(comp_val);

    return arma::min(mat_1,comp_mat);
}

inline
arma::mat
elem_min(const double comp_val, const arma::mat& mat_1)
{
    return elem_min(mat_1,comp_val);
}

//
// minimum between two values

inline
double
elem_min(const double comp_val_1, const double comp_val_2)
{
    return std::min(comp_val_1,comp_val_2);
}

//
// return the maximal-value in the matrix

inline
double
elem_max(const arma::mat& mat_1)
{
    return mat_1.max();
}

//
// pointwise maximum between two matrices

inline
arma::mat
elem_max(const arma::mat& mat_1, const arma::mat& mat_2)
{
    return arma::max(mat_1,mat_2);
}

//
// pointwise maximum between a matrix and a single value

inline
arma::mat
elem_max(const arma::mat& mat_1, const double comp_val)
{
    const int rows_1 = mat_1.n_rows;
    const int cols_1 = mat_1.n_cols;
    //
    arma::mat comp_mat(rows_1,cols_1);
    comp_mat.fill(comp_val);

    return arma::max(mat_1,comp_mat);
}

inline
arma::mat
elem_max(const double comp_val, const arma::mat& mat_1)
{
    return elem_max(mat_1,comp_val);
}

//
// maximum between two values

inline
double
elem_max(const double comp_val_1, const double comp_val_2)
{
    return std::max(comp_val_1,comp_val_2);
}

//
// sum over a dimension of a cube

inline
arma::mat
cube_sum(const arma::cube& cube_inp, const int which_dim)
{
    if (which_dim > 1 || which_dim < 0) {
        printf("error: cube_sum: unrecognized dim value; should be in (0,1)\n");
    }
    //
    const int dim_0 = cube_inp.n_rows;
    const int dim_1 = cube_inp.n_cols;
    const int dim_2 = cube_inp.n_slices;

    arma::mat ret;

    if (which_dim == 0) { // over rows
        ret.set_size(dim_0,dim_2);
        for (int i=0; i < dim_2; i++) {
            arma::mat mat_s = cube_inp.slice(i);
            ret.col(i) = arma::sum(mat_s,1);
        }
    } else { // dim == 1
        ret.set_size(dim_1,dim_2);
        for (int i=0; i < dim_2; i++) {
            arma::mat mat_s = cube_inp.slice(i);
            ret.col(i) = arma::trans(arma::sum(mat_s,0));
        }
    }
    //
    return ret;
}

//
// mirror R's approach to creating a matrix from an array:
// take each slice and vectorise that matrix

inline
arma::mat
cube_to_mat(const arma::cube& cube_inp)
{
    const int dim_0 = cube_inp.n_rows;
    const int dim_1 = cube_inp.n_cols;
    const int dim_2 = cube_inp.n_slices;

    arma::mat mat_out(dim_0*dim_1,dim_2);

    for(int k = 0; k < dim_2; k++) {
        mat_out.col(k) = arma::vectorise(cube_inp.slice(k));
    }

    return mat_out;
}
