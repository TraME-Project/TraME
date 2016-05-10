/*
 * Structures
 */
 
// for use in zeroin()
struct trame_opt_data {
    arma::mat* expUbarX;
    arma::mat* mubarX;
};

class trame_zeroin_data
{
    public:
        //
        arma::mat expUbarX;
        arma::mat mubarX;
        
    //private:
};
