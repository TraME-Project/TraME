################################################################################
##
##   Copyright (C) 2015 - 2016 Alfred Galichon
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
################################################################################

inversePWA_old <- function(a, B, C, k = 1)
{
    nbX = length(a)
    nbY = dim(B)[2]
    #
    vals = rep(0,nbX)
    for(x in 1:nbX){    
        sortB = sort(B[x,],index.return=T)
        sigma = sortB$ix
        
        b = sortB$x
        small_C = C[x,sigma]
        
        ylow = 1
        yup = nbY
        while(yup > ylow){
            ymid = ylow + (yup - ylow) %/% 2
            lhs = k * b[ymid] + sum(small_C * pmin(b[ymid],b))
            if(lhs == a[x]){
                yup = ylow = ymid
            }else if(lhs > a[x]){
                yup = ymid
            }else{
                ylow= ymid + 1
            }
        }
        if((ylow==1) & ( k*b[ylow]+sum(small_C * pmin(b[ylow],b)) >= a[x])){
            vals[x] = a[x] / (k+sum(small_C))
        }else{
            ysincluded = which((1:nbY) <= ylow)
            vals[x]= (a[x] - sum((small_C*b)[ysincluded])) / (k + sum(small_C) - sum(small_C[ysincluded]))
        }
    }
    #
    return(vals)
}

inversePWA <- function(a, B, C, k=1.0)
{
    #
    vals <- .Call("invPWA_R", a,B,C,k, PACKAGE = "TraME")$vals
    #
    return(c(vals))
}

tests_TraME <- function(nbDraws = 1e3)
{
    ptm = proc.time()
    #
    hash_arum <- tests_arum(notifications=FALSE,nbDraws=10*nbDraws)
    hash_equilibrium <- tests_equilibrium(notifications=FALSE,nbDraws=nbDraws)
    hash_estimation <- tests_estimation(notifications=FALSE)
    #
    time = proc.time() - ptm
    message(paste0('All tests completed. Overall time elapsed = ', round(time["elapsed"],5), 's.'))
    #
    hash_vals <- list(hash_arum=hash_arum,hash_equilibrium=hash_equilibrium,hash_estimation=hash_estimation)
    message(compare_hashvals(hash_vals))
    #
    return(hash_vals)
}

compare_hashvals <- function(hash_vals)
{
    all_hash <- .combine_hashvals(hash_vals)
    
    main_hash  <- all_hash$main_hash
    arum_hash  <- all_hash$arum_hash
    equil_hash <- all_hash$equil_hash
    estim_hash <- all_hash$estim_hash
    #
    temp_paste <- character(0)
    
    if(identical(main_hash$true,main_hash$actual)){
        conclusion = 'Test results are correct!\n'
    }else{
        conclusion = '*** CAUTION *** There is a problem with the results of:\n'
        #
        if(!identical(main_hash$true[1],main_hash$actual[1])){
            conclusion = paste0(conclusion,'arums tests, specifically: ')
            
            temp_paste <- character(0)
            for(jj in 1:length(arum_hash$test_names)){
                if(!identical(arum_hash$true[jj],arum_hash$actual[jj])){
                    if(length(temp_paste)==0){
                        temp_paste <- c(temp_paste, arum_hash$test_names[jj])
                    }else{
                        temp_paste <- c(temp_paste,", ",arum_hash$test_names[jj])
                    }
                }
                temp_paste <- c(temp_paste,".\n")
            }
            conclusion = paste0(conclusion,temp_paste)
        }
        
        if(!identical(main_hash$true[2],main_hash$actual[2])){
            conclusion = paste0(conclusion,'equilibrium tests, specifically: ')
            
            temp_paste <- character(0)
            for(jj in 1:length(equil_hash$test_names)){
                if(!identical(equil_hash$true[jj],equil_hash$actual[jj])){
                    if(length(temp_paste)==0){
                        temp_paste <- c(temp_paste, equil_hash$test_names[jj])
                    }else{
                        temp_paste <- c(temp_paste,", ",equil_hash$test_names[jj])
                    }
                }
                temp_paste <- c(temp_paste,".\n")
            }
            conclusion = paste0(conclusion,temp_paste)
        }
        
        if(!identical(main_hash$true[3],main_hash$actual[3])){
            conclusion = paste0(conclusion,'estimation tests, specifically: ')
            
            temp_paste <- character(0)
            for(jj in 1:length(estim_hash$test_names)){
                if(!identical(estim_hash$true[jj],estim_hash$actual[jj])){
                    if(length(temp_paste)==0){
                        temp_paste <- c(temp_paste, estim_hash$test_names[jj])
                    }else{
                        temp_paste <- c(temp_paste,", ",estim_hash$test_names[jj])
                    }
                }
                temp_paste <- c(temp_paste,".\n")
            }
            conclusion = paste0(conclusion,temp_paste)
        }
        conclusion = paste0(conclusion,'Please check.\n')
    }
    return(conclusion)
}

verify_signature <- function(save_output=FALSE,output_file=NA)
{
    #
    if(save_output==TRUE){
        if(is.character(output_file)==TRUE){
            output_hide <- capture.output(hash_vals <- suppressMessages(tests_TraME()),file=output_file)
        }else{
            output_hide <- capture.output(hash_vals <- suppressMessages(tests_TraME()),file="TraME_test_results.txt")
        }
    }else{
        output_hide <- capture.output(hash_vals <- suppressMessages(tests_TraME()))
    }
    #
    message(compare_hashvals(hash_vals))
}

.combine_hashvals <- function(hash_vals)
{
    # True values
    main_true_hash <- c("456eeafce1147f6f5de6b09158004f5f", "021138f554d6b48d72965d039ffb6f97", "1291c1bcbcfda7db348a15b247939bee")
    
    arum_logit_true_hash  <- "c000a71c4d4b72f0647767e11c744e8e"
    arum_probit_true_hash <- "dae91616d83fc496b5d427f5f545b5fe"
    arum_RUSC_true_hash   <- "28ffac81b67402b971fbad45d543b7a0"
    arum_RSC_true_hash    <- "7478b3d249b9d71048b4c0c1f3f01d05"
    arum_c_true_hash <- c(arum_logit_true_hash,arum_probit_true_hash,arum_RUSC_true_hash,arum_RSC_true_hash)
    
    equil_darum_true_hash       <- "72ebe32082ba64375136c822e6b16f00"
    equil_ipfp_true_hash        <- "b65e3d957e869c2b97bda61da11931f8"
    equil_nodalNewton_true_hash <- "bd5a46e73eca8c8f22df23dd71b6cd59"
    equil_arcNewton_true_hash   <- "bd5a46e73eca8c8f22df23dd71b6cd59"
    equil_maxW_true_hash        <- "020e418de3412180741c79b416ca709d"
    equil_jacobi_true_hash      <- "56d90c68734b45216e7189d7e4143336"
    equil_CLP_true_hash         <- "563827536851a806fa7e89aba4d36d21"
    equil_oapLP_true_hash       <- "71baf868ceac93dc3dde74e71c62c25a"
    equil_nash_true_hash        <- "46c423efaa759127592f9b7b18013609"
    equil_c_true_hash <- c(equil_darum_true_hash,equil_ipfp_true_hash,equil_nodalNewton_true_hash,equil_arcNewton_true_hash,
                           equil_maxW_true_hash,equil_jacobi_true_hash,equil_CLP_true_hash,equil_oapLP_true_hash,equil_nash_true_hash,
                           equil_nash_true_hash)
    # names of test functions
    arum_c_test_names <- c("test_Logit","test_Probit","test_RUSC","test_RSC")
    
    equil_c_test_names <- c("test_ipfp","test_nodalNewton","test_arcNewton","test_maxWelfare","test_jacobi","test_darum",
                            "test_cupidsLP","test_oapLP","test_eapNash")
    
    estim_c_test_names <- c("test_loglikelihood","test_mle","test_mme")
    # computed hash numbers
    main_hash_vals <- c(hash_vals$hash_arum$res_all_md5,hash_vals$hash_equilibrium$res_all_md5,hash_vals$hash_estimation$res_all_md5)
    
    arum_logit_hash_val <- hash_vals$hash_arum$res_logit_md5
    arum_probit_hash_val <- hash_vals$hash_arum$res_probit_md5
    arum_RUSC_hash_val <- hash_vals$hash_arum$res_RUSC_md5
    arum_RSC_hash_val <- hash_vals$hash_arum$res_RSC_md5
    arum_c_hash_val <- c(arum_logit_hash_val,arum_probit_hash_val,arum_RUSC_hash_val,arum_RSC_hash_val)
    
    equil_darum_hash_val <- hash_vals$hash_equilibrium$res_darum_md5
    equil_ipfp_hash_val <- hash_vals$hash_equilibrium$res_ipfp_md5
    equil_nodalNewton_hash_val <- hash_vals$hash_equilibrium$res_nodalNewton_md5
    equil_arcNewton_hash_val <- hash_vals$hash_equilibrium$res_arcNewton_md5
    equil_maxW_hash_val <- hash_vals$hash_equilibrium$res_maxW_md5
    equil_jacobi_hash_val <- hash_vals$hash_equilibrium$res_jacobi_md5
    equil_CLP_hash_val <- hash_vals$hash_equilibrium$res_CLP_md5
    equil_oapLP_hash_val <- hash_vals$hash_equilibrium$res_oapLP_md5
    equil_nash_hash_val <- hash_vals$hash_equilibrium$res_nash_md5
    equil_c_hash_val <- c(equil_darum_hash_val,equil_ipfp_hash_val,equil_nodalNewton_hash_val,equil_arcNewton_hash_val,
                           equil_maxW_hash_val,equil_jacobi_hash_val,equil_CLP_hash_val,equil_oapLP_hash_val,equil_nash_hash_val,
                           equil_nash_hash_val)
    #
    ret_main <- list(true = main_true_hash, actual = main_hash_vals)
    ret_arum <- list(true = arum_c_true_hash, actual = arum_c_hash_vals, test_names = arum_c_test_names)
    ret_equil <- list(true = equil_c_true_hash, actual = equil_c_hash_vals, test_names = equil_c_test_names)
    ret_estim <- list(test_names = estim_c_test_names)
    #
    ret <- list(main_hash = ret_main, arum_hash = ret_arum, equil_hash = ret_equil, estim_hash = ret_estim)
    return(ret)
} 