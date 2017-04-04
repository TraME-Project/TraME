/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
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
 * Kronecker product of sparse matrices
 *
 * Keith O'Hara
 * 05/30/2016
 */

arma::sp_mat kron_sp(const arma::sp_mat& A, const arma::sp_mat& B)
{
    int A_rows = A.n_rows;
    int A_cols = A.n_cols;
    int B_rows = B.n_rows;
    int B_cols = B.n_cols;
    
    arma::sp_mat ret_mat(A_rows*B_rows, A_cols*B_cols);
    
    int i,j;
    double A_val;
    
    arma::sp_mat::iterator it     = A.begin();
    arma::sp_mat::iterator it_end = A.end();

    for (; it != it_end; ++it) {
        i = it.row();
        j = it.col();
            
        A_val = (*it);
            
        if (A_val!=0){
            ret_mat.submat(i*B_rows, j*B_cols, (i+1)*B_rows-1, (j+1)*B_cols-1) = A_val * B;
        }
    }
    
    return ret_mat;
}

arma::sp_mat kron_sp(const arma::mat& A, const arma::sp_mat& B)
{
    int A_rows = A.n_rows;
    int A_cols = A.n_cols;
    int B_rows = B.n_rows;
    int B_cols = B.n_cols;
    
    arma::sp_mat ret_mat(A_rows*B_rows, A_cols*B_cols);
    
    int i,j;
    
    //#pragma omp parallel for private(i,j)
        for (j = 0; j<A_cols; j++) {
            for (i = 0; i < A_rows; i++) {
                if (A.at(i,j)!=0.0) {
                    ret_mat.submat(i*B_rows, j*B_cols, (i+1)*B_rows-1, (j+1)*B_cols-1) = A.at(i,j) * B;
                }
            }
        }
    
    return ret_mat;
}

arma::sp_mat kron_sp(const arma::sp_mat& A, const arma::mat& B)
{
    int A_rows = A.n_rows;
    int A_cols = A.n_cols;
    int B_rows = B.n_rows;
    int B_cols = B.n_cols;
    
    arma::sp_mat ret_mat(A_rows*B_rows, A_cols*B_cols);
    
    int i,j;
    double A_val;
    
    arma::sp_mat::iterator it     = A.begin();
    arma::sp_mat::iterator it_end = A.end();

        for (; it != it_end; ++it) {
            i = it.row();
            j = it.col();
            
            A_val = (*it);
            
            if (A_val!=0) {
                ret_mat.submat(i*B_rows, j*B_cols, (i+1)*B_rows-1, (j+1)*B_cols-1) = A_val * B;
            }
        }
    
    return ret_mat;
}

arma::sp_mat kron_sp(const arma::mat& A, const arma::mat& B)
{
    int A_rows = A.n_rows;
    int A_cols = A.n_cols;
    int B_rows = B.n_rows;
    int B_cols = B.n_cols;
    
    arma::sp_mat ret_mat(A_rows*B_rows, A_cols*B_cols);
    
    int i,j;

    //#pragma omp parallel for private(i,j)
        for (j = 0; j<A_cols; j++) {
            for (i = 0; i < A_rows; i++) {
                if (A.at(i,j)!=0.0) {
                    ret_mat.submat(i*B_rows, j*B_cols, (i+1)*B_rows-1, (j+1)*B_cols-1) = A.at(i,j) * B;
                }
            }
        }
    
    return ret_mat;
}
