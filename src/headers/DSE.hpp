/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2016 Alfred Galichon and the TraME Team
  ##
  ##   This file is part of the R package TraME.
  ##
  ##   The R package TraME is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   The R package TraME is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

// DSE market

class DSE
{
    public:
        // build objects
        bool ETU = false;
        bool LTU = false;
        bool NTU = false;
        bool TU  = false;

        bool need_norm;
        bool outsideOption;

        int nbX;
        int nbY;

        double sigma;

        arma::vec n;
        arma::vec m;

        mmf mmf_obj;
        transfers trans_obj;

        logit arums_G_none;
        logit arums_H_none;

        logit arums_G_logit;
        logit arums_H_logit;

        // member functions
        void build_ETU(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, arma::mat tau_inp, double* sigma_inp, bool need_norm_inp);
        void build_LTU(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, double* sigma_inp, bool need_norm_inp);
        void build_NTU(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, double* sigma_inp, bool need_norm_inp);
        void build_TU(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, double* sigma_inp, bool need_norm_inp);

        void trans();

    private:
        bool arum_general;
        bool arum_empirical;
        bool arum_none;
};

void MFE::build_TU(arma::vec n_inp, arma::vec m_inp, arma::mat phi, bool need_norm_inp, int arum_type)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    if (need_norm) {
        outsideOption = false;
    } else {
        outsideOption = false;
    }

    trans_obj.build_TU(phi);

    if (arum_type==1) {
        arums_G.build(nbX,nbY,sigma,outsideOption);
        arums_H.build(nbY,nbX,sigma,outsideOption);
    } else {
        
    }
}