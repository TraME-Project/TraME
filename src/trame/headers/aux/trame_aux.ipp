/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2016 the TraME Team:
  ##      Alfred Galichon
  ##      Keith O'Hara
  ##      Simon Weber
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
 * Aux functions
 *
 * Keith O'Hara
 * 08/08/2016
 */

// Find indices that correspond to maximum values 
inline arma::uvec which_max(const arma::mat* X, int which_dim)
{
    int i,j;
    
    int n = X->n_rows;
    int k = X->n_cols;

    arma::uvec max_vec;
     
    int max_ind = 0;
    double max_val = 0;
     
    if (which_dim==0) { // each column
        max_vec.set_size(k);
        for (j=0; j<k; j++) {
            max_val = (*X)(0,j);
            max_ind = 0;
            for (i=1; i<n; i++) {
                if ((*X)(i,j)>max_val) {
                    max_val = (*X)(i,j);
                    max_ind = i;
                }
            }
            max_vec(j) = max_ind;
        }
    } else { // each row
        max_vec.set_size(n);
        for (i=0; i<n; i++) {
            max_val = (*X)(i,0);
            max_ind = 0;
            for (j=0; j<k; j++) {
                if ((*X)(i,j) > max_val) {
                    max_val = (*X)(i,j);
                    max_ind = j;
                }
            }
            max_vec(i) = max_ind;
        }
    }
    //
    return max_vec;
}

// Generate an equally-spaced grid of integers
inline arma::uvec uvec_linspace (int a, int b)
{
    int i;
    int n_points = b - a + 1;
    
    arma::uvec ret(n_points);
    //
    for (i = 0; i<n_points; i++) {
        ret(i) = a + i;
    }
    //
    return ret;
}

/*
 * element-by-element addition;
 * written to mimic R's approach to element-by-element addition in vector-matrix cases
 *
 * Keith O'Hara
 * 08/18/2016
 */

inline arma::mat elem_add(const arma::mat& mat_1, const arma::mat& mat_2)
{
    int rows_1 = mat_1.n_rows;
    int rows_2 = mat_2.n_rows;

    int cols_1 = mat_1.n_cols;
    int cols_2 = mat_2.n_cols; 
    //
    bool same_rows, same_cols;

    if (rows_1 == rows_2) {
        same_rows = true;
    } else {
        same_rows = false;
    }

    if (cols_1 == cols_2) {
        same_cols = true;
    } else {
        same_cols = false;
    }
    //
    arma::mat ret;

    if (!same_rows && !same_cols) {
        printf("elem_sub: need matrices to agree on at least one dimension\n");
        return ret;
    }
    //
    // instead of loops we could also use arma::repmat, not sure which is quicker...
    int i;

    if (same_rows && same_cols) {
        ret = mat_1 + mat_2;
    } else if (same_rows && !same_cols) {
        if (cols_1==1) {
            ret.set_size(rows_1,cols_2);
            for (i=0; i<cols_2; i++) {
                ret.col(i) = mat_1 + mat_2.col(i);
            }
        } else if (cols_2==1) {
            ret.set_size(rows_1,cols_1);
            for (i=0; i<cols_1; i++) {
                ret.col(i) = mat_1.col(i) + mat_2;
            }
        } else {
            printf("elem_sub: need one of the matrices to be a column vector\n");
            return ret;
        }
    } else if (!same_rows && same_cols) {
        if (rows_1==1) {
            ret.set_size(rows_2,cols_1);
            for (i=0; i<rows_2; i++) {
                ret.row(i) = mat_1 + mat_2.row(i);
            }
        } else if (rows_2==1) {
            ret.set_size(rows_1,cols_1);
            for (i=0; i<rows_1; i++) {
                ret.row(i) = mat_1.row(i) + mat_2;
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

/*
 * element-by-element subtraction;
 * written to mimic R's approach to element-by-element subtraction in vector-matrix cases
 *
 * Keith O'Hara
 * 08/18/2016
 */

inline arma::mat elem_sub(const arma::mat& mat_1, const arma::mat& mat_2)
{
    int rows_1 = mat_1.n_rows;
    int rows_2 = mat_2.n_rows;

    int cols_1 = mat_1.n_cols;
    int cols_2 = mat_2.n_cols; 
    //
    bool same_rows, same_cols;

    if (rows_1 == rows_2) {
        same_rows = true;
    } else {
        same_rows = false;
    }

    if (cols_1 == cols_2) {
        same_cols = true;
    } else {
        same_cols = false;
    }
    //
    arma::mat ret;

    if (!same_rows && !same_cols) {
        printf("elem_sub: need matrices to agree on at least one dimension\n");
        return ret;
    }
    //
    // instead of loops we could also use arma::repmat, not sure which is quicker...
    int i;

    if (same_rows && same_cols) {
        ret = mat_1 - mat_2;
    } else if (same_rows && !same_cols) {
        if (cols_1==1) {
            ret.set_size(rows_1,cols_2);
            for (i=0; i<cols_2; i++) {
                ret.col(i) = mat_1 - mat_2.col(i);
            }
        } else if (cols_2==1) {
            ret.set_size(rows_1,cols_1);
            for (i=0; i<cols_1; i++) {
                ret.col(i) = mat_1.col(i) - mat_2;
            }
        } else {
            printf("elem_sub: need one of the matrices to be a column vector\n");
            return ret;
        }
    } else if (!same_rows && same_cols) {
        if (rows_1==1) {
            ret.set_size(rows_2,cols_1);
            for (i=0; i<rows_2; i++) {
                ret.row(i) = mat_1 - mat_2.row(i);
            }
        } else if (rows_2==1) {
            ret.set_size(rows_1,cols_1);
            for (i=0; i<rows_1; i++) {
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

/*
 * Hadamard product;
 * written to mimic R's approach to element-by-element multiplication in vector-matrix cases
 *
 * Keith O'Hara
 * 08/18/2016
 */

inline arma::mat elem_prod(const arma::mat& mat_1, const arma::mat& mat_2)
{
    int rows_1 = mat_1.n_rows;
    int rows_2 = mat_2.n_rows;

    int cols_1 = mat_1.n_cols;
    int cols_2 = mat_2.n_cols; 
    //
    bool same_rows, same_cols;

    if (rows_1 == rows_2) {
        same_rows = true;
    } else {
        same_rows = false;
    }

    if (cols_1 == cols_2) {
        same_cols = true;
    } else {
        same_cols = false;
    }
    //
    arma::mat ret;

    if (!same_rows && !same_cols) {
        printf("elem_prod: need matrices to agree on at least one dimension\n");
        return ret;
    }
    //
    // instead of loops we could also use arma::repmat, not sure which is quicker...
    int i;

    if (same_rows && same_cols) {
        ret = mat_1 % mat_2;
    } else if (same_rows && !same_cols) {
        if (cols_1==1) {
            ret.set_size(rows_1,cols_2);
            for (i=0; i<cols_2; i++) {
                ret.col(i) = mat_1 % mat_2.col(i);
            }
        } else if (cols_2==1) {
            ret.set_size(rows_1,cols_1);
            for (i=0; i<cols_1; i++) {
                ret.col(i) = mat_1.col(i) % mat_2;
            }
        } else {
            printf("elem_prod: need one of the matrices to be a column vector\n");
            return ret;
        }
    } else if (!same_rows && same_cols) {
        if (rows_1==1) {
            ret.set_size(rows_2,cols_1);
            for (i=0; i<rows_2; i++) {
                ret.row(i) = mat_1 % mat_2.row(i);
            }
        } else if (rows_2==1) {
            ret.set_size(rows_1,cols_1);
            for (i=0; i<rows_1; i++) {
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

/*
 * element-by-element division;
 * written to mimic R's approach to element-by-element division in vector-matrix cases
 *
 * Keith O'Hara
 * 08/18/2016
 */

inline arma::mat elem_div(const arma::mat& mat_1, const arma::mat& mat_2)
{
    int rows_1 = mat_1.n_rows;
    int rows_2 = mat_2.n_rows;

    int cols_1 = mat_1.n_cols;
    int cols_2 = mat_2.n_cols; 
    //
    bool same_rows, same_cols;

    if (rows_1 == rows_2) {
        same_rows = true;
    } else {
        same_rows = false;
    }

    if (cols_1 == cols_2) {
        same_cols = true;
    } else {
        same_cols = false;
    }
    //
    arma::mat ret;

    if (!same_rows && !same_cols) {
        printf("elem_div: need matrices to agree on at least one dimension\n");
        return ret;
    }
    //
    // instead of loops we could also use arma::repmat, not sure which is quicker...
    int i;

    if (same_rows && same_cols) {
        ret = mat_1 / mat_2;
    } else if (same_rows && !same_cols) {
        if (cols_1==1) {
            ret.set_size(rows_1,cols_2);
            for (i=0; i<cols_2; i++) {
                ret.col(i) = mat_1 / mat_2.col(i);
            }
        } else if (cols_2==1) {
            ret.set_size(rows_1,cols_1);
            for (i=0; i<cols_1; i++) {
                ret.col(i) = mat_1.col(i) / mat_2;
            }
        } else {
            printf("elem_div: need one of the matrices to be a column vector\n");
            return ret;
        }
    } else if (!same_rows && same_cols) {
        if (rows_1==1) {
            ret.set_size(rows_2,cols_1);
            for (i=0; i<rows_2; i++) {
                ret.row(i) = mat_1 / mat_2.row(i);
            }
        } else if (rows_2==1) {
            ret.set_size(rows_1,cols_1);
            for (i=0; i<rows_1; i++) {
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

inline double elem_min(const arma::mat& mat_1)
{
    double ret = mat_1.min();
    //
    return ret;
}

inline double elem_max(const arma::mat& mat_1)
{
    double ret = mat_1.max();
    //
    return ret;
}
