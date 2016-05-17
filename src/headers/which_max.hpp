/*
 * Find indices the correspond to maximum values
 *
 * Keith O'Hara
 * 05/08/2016
 */
 
arma::uvec which_max(const arma::mat* X, int which_dim)
{
    int i,j;
     
    int n = X->n_rows;
    int k = X->n_cols;

    arma::uvec max_vec;
     
    int max_ind = 0;
    double max_val = 0;
     
    if(which_dim==0){ // each column
        max_vec.set_size(k);
        for(j=0; j<k; j++){
            max_val = (*X)(0,j);
            max_ind = 0;
            for(i=1; i<n; i++){
                if((*X)(i,j)>max_val){
                    max_val = (*X)(i,j);
                    max_ind = i;
                }
            }
            max_vec(j) = max_ind;
        }
    }else{ // each row
        max_vec.set_size(n);
        for(i=0; i<n; i++){
            max_val = (*X)(i,0);
            max_ind = 0;
            for(j=0; j<k; j++){
                if((*X)(i,j)>max_val){
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