
inline
arma::mat
outer_prod_minus(const arma::vec& mat_1, const arma::vec& mat_2)
{
    arma::mat ret(mat_1.n_elem,mat_2.n_elem);

    for (int j=0; j < (int) mat_1.n_elem; j++) {
        ret.row(j) =  mat_2.t() - mat_1(j);
    }
    //
    return ret;
}

inline
arma::mat
outer_prod_plus(const arma::vec& mat_1, const arma::vec& mat_2)
{
    arma::mat ret(mat_1.n_elem,mat_2.n_elem);

    for (int j=0; j < (int) mat_1.n_elem; j++) {
        ret.row(j) = mat_1(j) + mat_2.t();
    }
    //
    return ret;
}

//
//

inline
arma::vec
vec_bind(const arma::mat& mat_1, const arma::mat& mat_2)
{
    return arma::join_cols( arma::vectorise(mat_1), arma::vectorise(mat_2) );
}

inline
arma::vec
vec_bind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3)
{
    return arma::join_cols( arma::join_cols(arma::vectorise(mat_1), arma::vectorise(mat_2)), arma::vectorise(mat_3) );
}

inline
arma::vec
vec_bind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3, const arma::mat& mat_4)
{
    return vec_bind( mat_1,mat_2, arma::join_cols(arma::vectorise(mat_3), arma::vectorise(mat_4)) );
}

inline
arma::vec
vec_bind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3, const arma::mat& mat_4, const arma::mat& mat_5)
{
    return vec_bind( mat_1,mat_2,mat_3, arma::join_cols(arma::vectorise(mat_4), arma::vectorise(mat_5)) );
}

inline
arma::vec
vec_bind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3, const arma::mat& mat_4, const arma::mat& mat_5, const arma::mat& mat_6)
{
    return vec_bind( mat_1,mat_2,mat_3,mat_4, arma::join_cols(arma::vectorise(mat_5), arma::vectorise(mat_6)) );
}

inline
arma::vec
vec_bind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3, const arma::mat& mat_4, const arma::mat& mat_5, const arma::mat& mat_6, const arma::mat& mat_7)
{
    return vec_bind( mat_1,mat_2,mat_3,mat_4,mat_5, arma::join_cols(arma::vectorise(mat_6), arma::vectorise(mat_7)) );
}

inline
arma::vec
vec_bind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3, const arma::mat& mat_4, const arma::mat& mat_5, const arma::mat& mat_6, const arma::mat& mat_7, const arma::mat& mat_8)
{
    return vec_bind( mat_1,mat_2,mat_3,mat_4,mat_5,mat_6, arma::join_cols(arma::vectorise(mat_7), arma::vectorise(mat_8)) );
}

//
//

inline
arma::cube
cube_bind(const arma::mat& mat_1, const arma::mat& mat_2)
{
    arma::cube ret(mat_1.n_rows,mat_1.n_cols,2);

    ret.slice(0) = mat_1;
    ret.slice(1) = mat_2;

    return ret;
}

inline
arma::cube
cube_bind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3)
{
    arma::cube ret(mat_1.n_rows,mat_1.n_cols,3);

    ret.slice(0) = mat_1;
    ret.slice(1) = mat_2;
    ret.slice(2) = mat_3;

    return ret;
}

inline
arma::cube
cube_bind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3, const arma::mat& mat_4)
{
    arma::cube ret(mat_1.n_rows,mat_1.n_cols,4);

    ret.slice(0) = mat_1;
    ret.slice(1) = mat_2;
    ret.slice(2) = mat_3;
    ret.slice(3) = mat_4;

    return ret;
}

inline
arma::cube
cube_bind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3, const arma::mat& mat_4, const arma::mat& mat_5)
{
    arma::cube ret(mat_1.n_rows,mat_1.n_cols,5);

    ret.slice(0) = mat_1;
    ret.slice(1) = mat_2;
    ret.slice(2) = mat_3;
    ret.slice(3) = mat_4;
    ret.slice(4) = mat_5;

    return ret;
}

inline
arma::cube
cube_bind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3, const arma::mat& mat_4, const arma::mat& mat_5, const arma::mat& mat_6)
{
    arma::cube ret(mat_1.n_rows,mat_1.n_cols,6);

    ret.slice(0) = mat_1;
    ret.slice(1) = mat_2;
    ret.slice(2) = mat_3;
    ret.slice(3) = mat_4;
    ret.slice(4) = mat_5;
    ret.slice(5) = mat_6;

    return ret;
}

inline
arma::cube
cube_bind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3, const arma::mat& mat_4, const arma::mat& mat_5, const arma::mat& mat_6, const arma::mat& mat_7)
{
    arma::cube ret(mat_1.n_rows,mat_1.n_cols,7);

    ret.slice(0) = mat_1;
    ret.slice(1) = mat_2;
    ret.slice(2) = mat_3;
    ret.slice(3) = mat_4;
    ret.slice(4) = mat_5;
    ret.slice(5) = mat_6;
    ret.slice(6) = mat_7;

    return ret;
}

inline
arma::cube
cube_bind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3, const arma::mat& mat_4, const arma::mat& mat_5, const arma::mat& mat_6, const arma::mat& mat_7, const arma::mat& mat_8)
{
    arma::cube ret(mat_1.n_rows,mat_1.n_cols,8);

    ret.slice(0) = mat_1;
    ret.slice(1) = mat_2;
    ret.slice(2) = mat_3;
    ret.slice(3) = mat_4;
    ret.slice(4) = mat_5;
    ret.slice(5) = mat_6;
    ret.slice(6) = mat_7;
    ret.slice(7) = mat_8;

    return ret;
}

