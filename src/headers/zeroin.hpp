/*
 * Slight modification of the zeroin.c file from Netlib
 *
 * Keith O'Hara
 * 05/10/2016
 */

#include <cmath> // for abs value

double zeroin(double ax, double bx, double (*f)(double x, const trame_zeroin_data& opt_data), const trame_zeroin_data& zeroin_data, double tol, int max_iter)
{
	double a,b,c;
	double fa;
	double fb;
	double fc;
		
	a = ax;  b = bx;  fa = (*f)(a,zeroin_data);  fb = (*f)(b,zeroin_data);
	c = a;   fc = fa;
	
	// check endpoints
	if(fa == 0.0){
		return b;
	}
	if(fb == 0.0){
		return b;
	}
	
	// otherwise begin iterations
	int iter = 0;
	
	double eps_temp = std::numeric_limits<double>::epsilon();
	double tol_act = 2*eps_temp*fabs(b) + tol/2;

	double p, q, prev_step, new_step;
	//register double t1,cb,t2;
	double t1,cb,t2; // Keith: register is deprecated as of C++-11
	
	new_step = (c - b)/2;
	
	while(fabs(new_step) > tol_act && iter < max_iter){
		iter++;
		prev_step = b-a;
			
		if( fabs(fc) < fabs(fb) ){
			a = b;  b = c;  c = a;
			fa=fb;  fb=fc;  fc=fa;
		}
		
		new_step = (c-b)/2;

		if( fabs(prev_step) >= tol_act	&& fabs(fa) > fabs(fb) ){
			
			cb = c - b;
						
			if( a==c ){
				t1 = fb/fa;
				p = cb*t1;
				q = 1.0 - t1;
			}else{
				q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
				p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
				q = (q-1.0) * (t1-1.0) * (t2-1.0);
			}
						
			if( p > 0.0 ){
				q = -q;
			}else{
				p = -p;
			}

			if( p < (0.75*cb*q-fabs(tol_act*q)/2) && p < fabs(prev_step*q/2) ){
				new_step = p/q;
			}
		}
				
		if( fabs(new_step) < tol_act ){
			if( new_step > 0.0 ){
				new_step = tol_act;
			}else{
				new_step = -tol_act;
			}
		}
						
		a = b;  fa = fb;
		b += new_step;  fb = (*f)(b,zeroin_data);
		if( (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) ){
			c = a;  fc = fa;
		}
	}
	//
	return b;
}

/*double zeroin_old(double ax, double bx, double (*f)(double x, void* opt_data), void* zeroin_data, double tol, int max_iter)
{
	double a,b,c;				// Abscissae
	double fa;					// f(a)
	double fb;					// f(b)
	double fc;					// f(c)
		
	a = ax;  b = bx;  fa = (*f)(a,&zeroin_data);  fb = (*f)(b,&zeroin_data);
	c = a;   fc = fa;
	
	// check endpoints
	if(fa == 0.0){
		return b;
	}
	if(fb == 0.0){
		return b;
	}
	
	// otherwise begin iterations
	int iter = 0;
	
	double eps_temp = std::numeric_limits<double>::epsilon();
	double tol_act = 2*eps_temp*fabs(b) + tol/2;

	double p, q, prev_step, new_step;
	register double t1,cb,t2;
	
	new_step = (c - b)/2;
	
	while(fabs(new_step) > tol_act && iter < max_iter){
		iter++;
		prev_step = b-a;
			
		if( fabs(fc) < fabs(fb) ){
			a = b;  b = c;  c = a;
			fa=fb;  fb=fc;  fc=fa;
		}
		
		new_step = (c-b)/2;

		if( fabs(prev_step) >= tol_act	&& fabs(fa) > fabs(fb) ){
			
			cb = c - b;
						
			if( a==c ){
				t1 = fb/fa;
				p = cb*t1;
				q = 1.0 - t1;
			}else{
				q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
				p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
				q = (q-1.0) * (t1-1.0) * (t2-1.0);
			}
						
			if( p > 0.0 ){
				q = -q;
			}else{
				p = -p;
			}

			if( p < (0.75*cb*q-fabs(tol_act*q)/2) && p < fabs(prev_step*q/2) ){
				new_step = p/q;
			}
		}
				
		if( fabs(new_step) < tol_act ){
			if( new_step > 0.0 ){
				new_step = tol_act;
			}else{
				new_step = -tol_act;
			}
		}
						
		a = b;  fa = fb;
		b += new_step;  fb = (*f)(b,&zeroin_data);
		if( (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) ){
			c = a;  fc = fa;
		}
	}
	
	return b;
}*/