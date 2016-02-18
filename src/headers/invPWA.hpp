arma::vec invPWA (arma::vec a, arma::mat B, arma::mat C)
{
    int nb_X = a.n_elem;
    int nb_Y = B.n_cols;
    //
    arma::vec Bx, vals(nb_X), small_C, b, b_mid(nb_Y), b_low(nb_Y), temp_vec;
    vals.zeros();
    arma::uvec Bx_sort_ind;
    
    int x, y_low, y_up, y_mid, y_incl;
    double small_C_sum, lhs, term_1, term_2;
    //
    for(x = 0; x < nb_X; x++){
        Bx = B.row(x);
        Bx_sort_ind = arma::sort_index(Bx);
        b = Bx.elem(Bx_sort_ind);
        //
        //small_C = C(x,arma::span(Bx_sort_ind));
        small_C = C.row(x);
        small_C = small_C.elem(Bx_sort_ind);
        small_C_sum = arma::sum(small_C);
        //
        y_low = 1;
        y_up  = nb_Y;
        
        while(y_up > y_low){
            y_mid = y_low + (y_up - y_low)/2.0;
            
            b_low.fill(b(y_low-1));
            b_mid.fill(b(y_mid-1));
            
            lhs = b(y_mid-1) + arma::sum(small_C % arma::min(b_mid,b));
            //
            if(lhs == a(x)){
                y_low = y_mid;
                y_up = y_mid;
            }
            else if(lhs > a(x)){
                y_up = y_mid;
            }
            else{
                y_low = y_mid + 1;
            }
            //
            if((y_low==1) && ( b(y_low-1) + arma::sum(small_C % arma::min(b_low,b)) >= a(x) )){
                vals(x) = a(x) / (1 + small_C_sum);
            }
            else{
                y_incl = y_low - 1;
                
                temp_vec = (small_C % b);
                term_1 = a(x) - arma::sum(temp_vec(arma::span(0,y_incl)));
                term_2 = 1 + sum(small_C) - arma::sum(small_C(arma::span(0,y_incl)));
                
                vals(x) = term_1 / term_2;
            }
        }
    }
    //
    return vals;
}