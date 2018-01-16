/*!
 *  \file gsta_opt.cpp gsta_opt.h
 *	\brief Generalized Statistical Thermodynamic Adsorption (GSTA) Optimization Routine
 *  \author Austin Ladshaw
 *	\date 12/17/2013
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "gsta_opt.h"


//Function to round numbers down at 0.5 or less and up for all else
int roundIt(double d)
{
	int a=0; double intpart;
	modf(d,&intpart) > 0.5 ? a=ceil(d) : a=floor(d);
	return a;
}

// The Two-Fifths Compromise Routine
int twoFifths(int m)
{
	int n=0; double D;
	D = (2 * m) / 5 ;
	n = roundIt(D);
	return n;
}

//Function to return the order of magnitude of a number
int orderMag(double x)
{
	int m=0; double y=fabs(x);
	if (x == 0.0)
		return 1;
	if (x>=10)
	{
		while(y>=10)
		{
			y = y / 10;
			m++;
		}
	}
	else
	{
		while(y<1)
		{
			y = y * 10;
			m--;
		}
	}
	return m;
}

//Function returns the minimum int in a vector
int minValue(std::vector<int> &array)
{
	int min = array[0];
	unsigned long int length = array.size();
	for(int i = 1; i<length; i++)
	{
		if(array[i] < min) min = array[i];
	}
	return min;
}

//Function to find lowest (minimum) value in array
int minIndex(std::vector<double> &array)
{
     double min = array[0];
     unsigned long int length = array.size();
     int a = 0;

     for(int i = 1; i<length; i++)
     {
          if(array[i] < min)
          {
                min = array[i];
                a = i;
          }
     }
     return a;
}

//Function to determine the average of an array of ints and return a rounded int
int avgPar(std::vector<int> &array)
{
	int avg = 0; double avgd=0;
	unsigned long int length = array.size();
	for(int i=0; i<length; i++)
	{
		avgd = avgd + array[i];
	}
	avgd = avgd / array.size();
	avg = roundIt(avgd);
	return avg;
}

//Function to determine the average of an array of doubles
double avgValue(std::vector<double> &array)
{
	double avg=0;
	unsigned long int length = array.size();
	for (int i=0; i<length; i++)
	{
		avg = avg + array[i];
	}
	avg = avg / array.size();
	return avg;
}

//Function to determine the weighted average of an array based on the E.Norms
double weightedAvg(double *enorm, double *x, int n)
{
	double wAvg = 0, tnorm = 0;
	double frac[n], w[n];

	for (int i=0; i<n; i++)
	{
		tnorm = tnorm + enorm[i];
	}
	for (int i=0; i<n; i++)
	{
		frac[i] = enorm[i] / tnorm;
		w[i] = (1 - frac[i]) / (n - 1);
		wAvg = wAvg + (x[i] * w[i]);
	}

	return wAvg;
}

//Calculate the Coefficient of Determination (i.e. R_Squared) of a Linear Regression
double rSq(double *x, double *y, double slope, double vint, int m_dat)
{
	double rSq = 0, avgY = 0;
	double f[m_dat];
	double tot = 0, res = 0;

	for (int i=0; i<m_dat; i++)
	{
		f[i] = (x[i] * slope) + vint;
		avgY = avgY + y[i];
	}
	avgY = avgY / m_dat;

	for (int i=0; i<m_dat; i++)
	{
		tot = tot + pow( (y[i] - avgY) ,2);
		res = res + pow( (y[i] - f[i]) , 2);
	}
	rSq = 1 - (res / tot);

	return rSq;
}

//Check the smoothness of the parameter array
bool isSmooth(double *par, void *data)
{
	GSTA_OPT_DATA *dat = (GSTA_OPT_DATA *) data;
	bool ans; int count=0;
	std::vector<double> first, second;

	ans = true;
	if (dat->qmax != 0)
	{
		//checking interior points
		for(int i=1; i<(dat->n_par-1); i++)
		{
			first.push_back( (orderMag(fabs(par[(i+1)])) - orderMag(fabs(par[(i)]))) / 1.0 );
			second.push_back( (orderMag(fabs(par[(i+1)])) - (2*orderMag(fabs(par[i]))) + orderMag(fabs(par[(i-1)]))) / 1.0 );

			if (count !=0 && fabs( second[(count)] - second[(count-1)] ) > 6 )
			{
				ans = false;
				break;
			}
			count++;
		}
	}
	else
	{
		//checking interior points
		for(int i=2; i<(dat->n_par-1); i++)
		{
			first.push_back( ( (1.0*orderMag(fabs(par[(i+1)])) - 1.0*orderMag(fabs(par[(i)]))) / 1.0 ) );
			second.push_back( ( (orderMag(fabs(par[(i+1)])) -
					(2.0*orderMag(fabs(par[i]))) + orderMag(fabs(par[(i-1)]))) / 1.0 ) );

			if (count !=0 && fabs( second[(count)] - second[(count-1)] ) > 6)
			{
				ans = false;
				break;
			}
			count++;
		}
	}

	first.clear();
	second.clear();
	return ans;
}

//Function to perform orthogonal linear regression and is used in conjunction with eduGuess_GSTA
void orthoLinReg(double *x, double *y, double *par, int m_dat, int n_par)
{
	double sx=0, x2=0, y2=0, sy=0, xy=0, r=0, base=0, slope=0, vint=0, sp=0, sn=0, check=0;
	int nm=0;
	for(int i = 0; i < m_dat; i++)
	{
		nm++;
		sx = sx + x[i];
		sy = sy + y[i];
		x2 = x2 + pow(x[i],2);
		y2 = y2 + pow(y[i],2);
		xy = xy + (x[i] * y[i]);
	}
	base = (nm * xy) - (sx * sy);
	r = ( ( (nm * x2) - (pow(sx,2)) ) - ( (nm * y2) - (pow(sy,2)) ) ) / base ;
	if(base > 0)
	{
		slope = ( (-1 * r) + sqrt((pow(r,2) + 4)) ) / 2 ;
	}
	else
	{
		slope = ( (-1 * r) - sqrt((pow(r,2) + 4)) ) / 2 ;
	}
	sp = ( (-1 * r) + sqrt((pow(r,2) + 4)) ) / 2 ;
	sn = ( (-1 * r) - sqrt((pow(r,2) + 4)) ) / 2 ;
	check = sp * sn;	//Should equal -1 if done correctly
	if (check < -1.1 || check > -0.9)
	{
		mError(ortho_check_fail); std::cout << "orthoLinReg Error: Ignore if Check = -1 -> Check = " << check << std::endl;
	}
	vint = (sy - (slope * sx)) / nm;

	//Determines which parameter needs estimating based on whether or not qmax is given
	if (n_par == 1) //if qmax is given, n_par == 1
	{
		par[0] = fabs(vint / slope);
	}
	else if (n_par == 2) //if qmax is not given, n_par == 2
	{
		par[0] = fabs(1 / vint);
		par[1] = fabs(vint / slope);
	}
	else //performed if n_par is other
	{
		par[0] = vint;
		par[1] = slope;
	}
}

//Educated Guessing Algorithm
void eduGuess(double *P, double *q, double *par, int k, int m_dat, void *data)
{
	GSTA_OPT_DATA *dat = (GSTA_OPT_DATA *) data;
	double x[m_dat], y[m_dat], slope = 0, par_sum = 0;

	if (dat->n_par==1 || (dat->qmax==0 && dat->n_par==2))
	{
		for(int i = 0; i < m_dat; i++)
		{
			x[i] = 1/P[i];
			y[i] = 1/q[i];
		}
		orthoLinReg(x, y, par, m_dat, dat->n_par);
	}
	else
	{
		for(int i = 0; i < dat->n_par; i++)
		{
			//Set the new parameters to a projected order of magnitude based on other optimized parameters
			if(i < (dat->n_par - 2))
			{
				if(k < 2)
				{
					par[i] = fabs(dat->all_pars[dat->iso][(k - 1)][i]);
					par_sum = par_sum + orderMag(fabs(par[i]));
				}
				else
				{
					slope = orderMag(fabs(dat->all_pars[dat->iso][(k-1)][i])) -
							orderMag(fabs(dat->all_pars[dat->iso][(k-2)][i]));
					par[i] = pow(10, (slope + orderMag(fabs(dat->all_pars[dat->iso][(k-1)][i])) ) );
					par_sum = par_sum + orderMag(fabs(par[i]));
				}
			}
			//Set the additional parameter argument to a projection based off of other optimized parameters
			else if (i < (dat->n_par - 1))
			{
				par[i] = pow(10, (slope + orderMag(fabs(dat->all_pars[dat->iso][(k-1)][i])) ) );
				par_sum = par_sum + orderMag(fabs(par[i]));
			}
			else
			{
				if ( (dat->n_par > 2 && dat->qmax != 0 ) || (dat->n_par > 3 && dat->qmax == 0) )
				{
					slope = orderMag(fabs(par[(i-1)])) - orderMag(fabs(par[(i-2)]));
					par[i] = pow(10,(slope + orderMag(fabs(par[(i-1)])) ) );
				}
				else
				{
					//Remove qmax from par_sum if it was adjustable
					if (dat->qmax == 0) par_sum = par_sum - orderMag(fabs(par[0]));
					par[i] = pow(10, roundIt((par_sum / (dat->n_par - 1))));
				}
			}
		}
	}
}

// model function: The GSTA Isotherm q = (qmax/m) * (sum(n*Kn*P^n)/(1+sum(Kn*P^n)))
double gstaFunc( double p, const double *K, double qmax, int n_par)
{
	double Y=0, Numerator=0, Denominator=1;
	if (qmax != 0)
	{
		for(int n=1; n <= n_par; n++)
		{
			Numerator = Numerator + (n * fabs(K[(n-1)]) * pow(p,n));
			Denominator = Denominator + (fabs(K[(n-1)]) * pow(p,n));
		}
		Y = (qmax / (n_par)) * (Numerator / Denominator);
	}
	else
	{
		for(int n=1; n < n_par; n++)
		{
			Numerator = Numerator + (n * fabs(K[n]) * pow(p,n));
			Denominator = Denominator + (fabs(K[n]) * pow(p,n));
		}
		Y = (fabs(K[0]) / (n_par - 1)) * (Numerator / Denominator);
	}
	return Y;
}

//Function used to compute the gsta objective function for fitness determination
double gstaObjFunc(double *t, double *y, double *par, int m_dat, void *data)
{
	GSTA_OPT_DATA *dat = (GSTA_OPT_DATA *) data;
	double objF, sigma = 0;
	for (int i=0;i<m_dat;i++)
	{
		sigma = sigma + pow(((y[i] - gstaFunc(t[i],par,dat->qmax, dat->n_par))/y[i]),2);
	}
	objF = sqrt(sigma / (m_dat - (dat->n_par) - 1));
	return objF;
}

//Evaluation of the vectoral residuals of the GSTA isotherm model
void eval_GSTA(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
	GSTA_OPT_DATA *dat = (GSTA_OPT_DATA *) data;
	double q_eval;
	double Numerator, Denominator;
	for(int i=0; i<m_dat; i++)
	{
		Numerator = 0; Denominator = 1;
		if(dat->qmax != 0)
		{
			for(int n=1; n <= dat->n_par; n++)
			{
				Numerator = Numerator + (n * fabs(par[(n-1)]) * pow(dat->P[dat->iso][i],n));
				Denominator = Denominator + (fabs(par[(n-1)]) * pow(dat->P[dat->iso][i],n));
			}
			q_eval = (dat->qmax / (dat->n_par)) * (Numerator / Denominator);
			fvec[i] = dat->q[dat->iso][i] - q_eval;
		}
		else
		{
			for(int n=1; n < dat->n_par; n++)
			{
				Numerator = Numerator + (n * fabs(par[n]) * pow(dat->P[dat->iso][i],n));
				Denominator = Denominator + (fabs(par[n]) * pow(dat->P[dat->iso][i],n));
			}
			q_eval = (fabs(par[0]) / (dat->n_par - 1)) * (Numerator / Denominator);
			fvec[i] = dat->q[dat->iso][i] - q_eval;
		}
	}
}

int gsta_optimize(const char* fileName)
{
	//Check for NULL argument
    if (fileName==NULL)
    {
    	std::cout << "Enter file name and extension for input file in this directory\n";
    	std::cout << "(Example: 'Input.txt')\n";
    	std::string NAME;
    	std::cin >> NAME;
    	fileName = NAME.c_str();
    }
    std::ifstream myfile( fileName );
	//Check to see if file exists
    if (myfile.good()==false)
	{
		mError(file_dne);
		return -1;
	}

	//Declarations
    std::vector<int> m_dat; 						//# of data pairs to be determined by list size
    double n; int j=0, k=0;							//Declaration of n for reading and j and k for indexing
    double time;									//Used to calculate the time it took to run the program
    int n_iso;										//First to be read in. Determines number of isotherms.
    std::vector<double> iso_temp;					//Stores the system temperature for each isotherm
    std:: vector<int> num_par, best_num_par;		//Stores the best fit number of params for each isotherm
    int max_num_par;								//Stores the max number of parameters valid for all isotherms
    int best_num_all;								//Stores the best number of parameters valid for all isotherms
    GSTA_OPT_DATA dat;								//Used for large sets of information to be passed to functions
    int QMAX = 0;									//Index for the qmax loop
    std:: vector<double> dHo, dSo;					//Vectors to store the energy terms to be determined
    std::vector<std::vector<double> > dGo;			//Vector to store Gibb's energy
    int max_par;									//Flag variable for inside main loop
    double temp_sum;								//Temp variable for determination of qmax
    double temp_norm;								//Temp variable for determination of avg norm for qmax
    int temp_int;									//Temp variable for determination of qmax
    std::vector<double> avg_norm;					//Vector used for weighted average of qmax
    int abs_max;									//Stores the absolute max par by Gibbs' Phase Rule
    double abs_max_q;								//Stores the highest q value recorded from input data

    //Initializes the first column for all objects
    dat.iso = 0;
    dat.P.push_back(std::vector<double> ());
    dat.q.push_back(std::vector<double> ());
    dGo.push_back(std::vector<double> ());
    dat.total_eval = 0;
    max_par = 0;
    temp_sum = 0;
    temp_int = 0;
    temp_norm = 0;
    time = clock();
    abs_max = 6;
    abs_max_q = 0;

    //Read input from text file
    myfile >> n_iso;
    myfile >> dat.qmax;
    int a=0;
    for(int i=0; i<n_iso; i++)
    {
    	a=0;
    	myfile >> n; iso_temp.push_back(n);
    	myfile >> n; m_dat.push_back(n);

    	for(int k=0; k<(2*m_dat.at(i)); k++)
    	{
    		if(j == 0)
    		{
    			myfile >> n; dat.P[i].push_back(n);
    			j++;
    		}
    		else if(j == 1)
    		{
    			myfile >> n; dat.q[i].push_back(n);
    			j--;a++;
    		}
    		else { mError(indexing_error); return-1;}
    	}
    	dat.P.push_back(std::vector<double> ());
    	dat.q.push_back(std::vector<double> ());
    }
    myfile.close();

    //Start outer most loop for determination of a single qmax for all isotherms
    do
    {
    	//Specific initializations
    	dat.all_pars.push_back(std::vector<std::vector<double> > ());
    	dat.Fobj.push_back(std::vector<double> ());
    	dat.best_par.push_back(std::vector<double> ());
    	dat.Kno.push_back(std::vector<double> ());
    	dat.norms.push_back(std::vector<double> ());

    	//Open the output text files
    	FILE *output, *best, *norm, *params, *heats;
    	output = fopen("output/All_GSTA_Results.txt","w+");
    	best = fopen("output/Best_GSTA_Results.txt", "w+");
    	norm = fopen("output/GSTA_Norms_v_Fobj.txt","w+");
    	params = fopen("output/GSTA_Parameter_Results.txt","w+");
    	heats = fopen("output/GSTA_Energy_Results.txt","w+");
		
		if (output == nullptr)
		{
			system("mkdir output");
			output = fopen("output/All_GSTA_Results.txt","w+");
			best = fopen("output/Best_GSTA_Results.txt", "w+");
			norm = fopen("output/GSTA_Norms_v_Fobj.txt","w+");
			params = fopen("output/GSTA_Parameter_Results.txt","w+");
			heats = fopen("output/GSTA_Energy_Results.txt","w+");
		}

    	//Start outer loop for all isotherms
    	dat.iso=0;
    	do
    	{
    		temp_sum = 0; temp_int = 0; temp_norm = 0;
    		//Declaration for arrays of data pairs
    		double t[m_dat[dat.iso]],y[m_dat[dat.iso]];

    		for(int i=0;i<m_dat[dat.iso];i++)
    		{
    			t[i] = dat.P[dat.iso].at(i);
    			y[i] = dat.q[dat.iso].at(i);
    			if (dat.q[dat.iso][i] > abs_max_q) abs_max_q = dat.q[dat.iso][i];
    		}

    		// auxiliary parameters
    		lm_status_struct status;
    		lm_control_struct control = lm_control_double;
    		// monitor status (+1) and parameters (+2)
    		control.printflags = 0;

    		if (dat.qmax != 0) 	//Initializes the number of parameters
    		{
    			dat.n_par = 1;
    			QMAX = 1;
    		}
    		else 				//If qmax is unknown, it becomes a parameter
    		{
    			dat.n_par = 2;
    			QMAX = 0;
    			std::cout << "\nqmax is unknown - Code will be run twice to estimate this parameter..." << std::endl;
    		}

    		//Reinitialized indices
    		j=0,k=0;
    		//Inner loop for the increasing number of parameters for each isotherm
    		do
    		{
    			if (dat.iso > 0 && max_par < dat.n_par)
    			{
    				std::cout << "Maximum plausible number of parameters reached..." << std::endl;
    				break;
    			}
    			// Parameter Vector of changing size must be declared here to reinitialize with each loop
    			double par[dat.n_par];

    			//Make the educated guess for the parameter vector starting points
    			eduGuess(t, y, par, k, m_dat[dat.iso], (void *)&dat);

    			// perform the optimization
    			printf("-------------------------------------------------------------------\n");
    			printf( "Optimizing isotherm %i for %i parameters, please wait...\n",(dat.iso+1), dat.n_par );
    			lmmin(dat.n_par, par, m_dat[dat.iso], (void *)&dat, eval_GSTA, &control, &status, lm_printout_std);

    			//Store the number of evaluations at each call
    			dat.total_eval = dat.total_eval + status.nfev;
    			printf("\nStatus after %d function evaluations:\n  %s\n\n", status.nfev, lm_infmsg[status.info] );

    			//Checks smoothness of the parameters for all n...
    			if ((dat.qmax != 0 && dat.n_par > 2) || (dat.qmax == 0 && dat.n_par > 3))
    			{
    				if(isSmooth(par, (void *)&dat) == false)
    				{
    					std::cout << "Non-smoothness encountered in the parameter solution vector...\n" << std::endl;
    					break;
    				}
    			}

    			//Stops the loop after fit stops improving or if the increase in norms is greater than 2.5% if qmax is known
    					//and 10% if qmax is unknown
    			if((k!=0 && dat.qmax!=0 && (status.fnorm - dat.norms[dat.iso].at((k-1))) >
    					(dat.norms[dat.iso].at((k-1)) * 0.025))
    						|| (k!=0 && dat.qmax==0 && (status.fnorm - dat.norms[dat.iso].at((k-1))) >
    								(dat.norms[dat.iso].at((k-1))) * 0.100 ))
    			{
    				std::cout << "E.Norm for isotherm " << (dat.iso+1)
    						<< " is no longer improving after " << (dat.n_par-1) << " parameters.\n" <<  std::endl;
    				break;
    			}

    			//Creating an average of the adjustable qmax values
    			if(dat.qmax == 0 && fabs(par[0]) > y[(m_dat[dat.iso]-1)] )
    			{
    				temp_sum = temp_sum + fabs(par[0]);
    				temp_norm = temp_norm + status.fnorm;
    				temp_int++;
    			}

    			//Calculate a store the norms and Func objectives
    			dat.Fobj[dat.iso].push_back( gstaObjFunc(t, y, par, m_dat[dat.iso], (void *)&dat) );
    			dat.norms[dat.iso].push_back(status.fnorm);

    			//The following loop stores all the parameters into a single array and the current parameters
    			for (int i=0;i<dat.n_par;i++)
    			{
    				dat.all_pars.at(dat.iso).push_back(std::vector<double> ());
    				dat.all_pars[dat.iso][k].push_back(par[i]);
    			}

    			// print all results to a file
    			fprintf(output, "Results for isotherm %i at %12g K:\n",(dat.iso+1),iso_temp[dat.iso] );
    			fprintf(output,"\nOptimized Parameters:\n");
    			if (dat.qmax != 0)
    			{
    				fprintf(output,"qmax =\t%12g\n", (dat.qmax));

    				for (int i = 0; i < dat.n_par; ++i)
    				{
    					fprintf(output,"K[%i] =\t%12g\n", (i+1), fabs(par[i]));
    				}
    			}
    			else
    			{
    				fprintf(output,"qmax =\t%12g\n", fabs(par[0]));
    				for (int i = 1; i < dat.n_par; ++i)
    				{
    					fprintf(output, "K[%i] =\t%12g\n", i, fabs(par[i]));
    				}
    			}
    			fprintf(output,"\nEuclidean Norm:\t%12g\n", status.fnorm );
    			if (dat.qmax !=0) fprintf(output,"Objective Function:\t%12g\n", dat.Fobj[dat.iso].at((dat.n_par-1)) );
    			else fprintf(output,"Objective Function:\t%12g\n", dat.Fobj[dat.iso].at((dat.n_par-2)) );
    			fprintf(output,"\nOptimized Results are as follows:\n");
    			for (int i = 0; i < m_dat[dat.iso]; ++i)
    			{
    				fprintf(output, "P[%i]=\t%12g\t q_obs=\t%12g\t q_eval=\t%12g\t Residual=\t%12g\n",
    						i, t[i], y[i], gstaFunc(t[i],par, dat.qmax, dat.n_par),
    							y[i] - gstaFunc(t[i],par, dat.qmax, dat.n_par) );
    			}
    			fprintf(output, "\n\n");

    			dat.n_par++;
    			k++;
    			dat.all_pars.push_back(std::vector<std::vector<double> > ());
    		}while(dat.n_par<(m_dat[dat.iso]-2) && dat.n_par <= ((twoFifths(m_dat[dat.iso])) + 1) && dat.n_par<=abs_max );

    		num_par.push_back((dat.n_par-1));
    		std::cout << "Maximum number of parameters for isotherm " <<
    				(dat.iso+1) <<" is " << num_par.at(dat.iso) << std::endl << std::endl;
    		max_par = dat.n_par - 1;
    		dat.opt_qmax.push_back( (temp_sum / temp_int) );
    		avg_norm.push_back( (temp_norm /  temp_int) );

    		//Push back vector<vectors< > > for next iterations
    		dat.Fobj.push_back(std::vector<double> ());
    		dat.norms.push_back(std::vector<double> ());
    		dat.all_pars.push_back(std::vector<std::vector<double> > ());

    		dat.iso++;
    	}while(dat.iso<n_iso);

    	//Determines the maximum number of parameters that is valid for all isotherms (based off of norms)
    	max_num_par = minValue(num_par);
    	printf("-------------------------------------------------------------------\n");
    	printf("-------------------------------------------------------------------\n");
    	std::cout << "The maximum number of parameters valid for all isotherms is " << max_num_par << std::endl;

    	//Determines the best number of parameters that can describe each isotherm (based off of Fobj)
    	if (dat.qmax!=0)
    	{
    		for(int i=0; i<n_iso; i++)
    		{
    			dat.Fobj[i].resize((max_num_par));
    			best_num_par.push_back((minIndex(dat.Fobj[i])+1));
    		}
    	}
    	else
    	{
    		for(int i=0; i<n_iso; i++)
    		{
    			dat.Fobj[i].resize((max_num_par-1));
    			best_num_par.push_back((minIndex(dat.Fobj[i])+2));
    		}
    	}

    	//Determine the best number of parameters that works for all isotherms
    	best_num_all = avgPar(best_num_par);
    	std::cout << "and the best number of parameters for all isotherms is " << best_num_all << std::endl;
    	std::cout << "\nTotal Function Evaluations: " << dat.total_eval << std::endl;

    	//Reset n_par to be used correctly in characteristic equation
    	dat.n_par = best_num_all;

    	//This stores the best fit parameters into the vector object best_par[i][j]
    	if (dat.qmax!=0)
    	{
    		for(int i=0; i<n_iso; i++)
    		{
    			for(int j=0; j<best_num_all; j++)
    			{
    				dat.best_par[i].push_back(dat.all_pars[i][(best_num_all-1)].at(j));
    				dat.Kno[i].push_back(fabs(dat.best_par[i][j]) * pow(Po, (j+1)));
    				//cout << "Ko[" << (i+1) << "][" << (j+1) << "] = " << dat.Kno[i][j] << endl;
    			}
    			//cout << endl;
    			dat.best_par.push_back(std::vector<double> ());
    			dat.Kno.push_back(std::vector<double> ());
    		}
    	}
    	else
    	{
    		for(int i=0; i<n_iso; i++)
    		{
    			for(int j=0; j<best_num_all; j++)
    			{
    				dat.best_par[i].push_back(dat.all_pars[i][(best_num_all-2)].at(j));
    				dat.Kno[i].push_back(fabs(dat.best_par[i][j]) * pow(Po, (j+1)));
    			}
    			dat.best_par.push_back(std::vector<double> ());
    			dat.Kno.push_back(std::vector<double> ());
    		}
    	}

    	//Only calculate the energies on the last run through (when qmax is not adjustable)
    	if (dat.qmax!=0 && n_iso > 1)
    	{
    		fprintf(heats, "Standard Molar Enthalpy and Entropy of Adsorption with Regression Analysis\n\n");
    		fprintf(heats, "n\tdHo[n]\tdSo[n]\tdHo[n]/n\tdSo[n]/n\t\n");
    		fprintf(heats, "(molecules)\t(J/mol)\t(J/(K*mol))\t(J/mol/molecule)\t(J/(K*mol)/molecule)\tR_Squared\n");
    		fprintf(heats, "-----------------------------------------------------------"
    				"----------------------------------------------\n");
    		//Declare temp X,Y, & par objects to be used in determining energies with orthoLinReg
    		double X[n_iso], Y[n_iso], heat_par[2], R_SQ = 0;
    		//Loops for determination of the energy terms and writing to output file
    		for (int i=0; i<best_num_all; i++)
    		{
    			for (int j=0; j<n_iso; j++)
    			{
    				X[j] = 1 / iso_temp[j];
    				Y[j] = log(fabs(dat.Kno[j][i]));
    			}
    			//At this line: X & Y have been filled for one LinReg
    			//NOTE: heat_par[0] = (dSo/R) & heat_par[1] = (-dHo/R)
    			orthoLinReg(X,Y,heat_par,n_iso,0);
    			dSo.push_back(heat_par[0] * R);
    			dHo.push_back(heat_par[1] * -R);
    			R_SQ = rSq(X,Y,heat_par[1], heat_par[0], n_iso);
    			fprintf(heats, "%6i\t%12g\t%12g\t%12g\t%12g\t%12g\t\n",
    					(i+1), dHo[i], dSo[i], (dHo[i]/(i+1)), (dSo[i]/(i+1)), R_SQ);
    		}
    	}

    	j=0;
    	fprintf(params, "Parameter Results and Standard Molar Gibb's Free Energy\n\n");
    	//Loop to write out all best results to a file
    	do
    	{
    		//Reinitializes arrays to be used in gstaFunc call
    		double t[m_dat[j]],y[m_dat[j]];
    		double par[best_num_all];

    		for(int i=0;i<best_num_all;i++)
    		{
    			par[i] = dat.best_par[j][i];
    		}

    		//Writes out the best results to a file
    		fprintf(best, "Results for isotherm %i at %12g K:\n", (j+1), iso_temp[j]);
    		fprintf(best,"\nOptimized Parameters:\n");


    		if (dat.qmax != 0 && n_iso > 1)
    		{
    			fprintf(params, "Temp= %6gK\tqmax= %6g[q]\tHe= %10g[q]/kPa\tE.Norm= %10g\n",
    					iso_temp[j], dat.qmax, ( (dat.qmax * dat.Kno[j][0]) / (best_num_all * Po) ),
    						dat.norms[j][(best_num_all-1)]);
    			fprintf(params, "---------------------------------------------"
    					"---------------------------------------------\n");
    			fprintf(params, "n\tK[n](kPa^-n)\tKo[n]\tdGo[n](J/mol)\tdGo[n]/n(J/mol/molecule)\n");
    			fprintf(params, "----------------------------------------------"
    					"--------------------------------------------\n");

    			fprintf(best,"qmax =\t%12g\n", (dat.qmax));

    			for (int i = 0; i < best_num_all; ++i)
    			{
    				dGo[j].push_back( (dHo[i] - (iso_temp[j] * dSo[i])) );
    				fprintf(best,"K[%i] =\t%12g\n", (i+1), fabs(par[i]));
    				fprintf(params, "%12i\t%12g\t%12g\t%12g\t%12g\n",
    						(i+1), fabs(dat.best_par[j][i]), dat.Kno[j][i], dGo[j][i], (dGo[j][i]/(i+1)) );
    			}
    			fprintf(best,"\nEuclidean Norm:\t%12g\n", dat.norms[j][(best_num_all-1)] );
    			fprintf(best,"Objective Function:\t%12g\n", dat.Fobj[j][(best_num_all-1)]);
    		}
    		else if (dat.qmax != 0 && n_iso == 1)
    		{
    			fprintf(best,"qmax =\t%12g\n", dat.qmax);
    			for (int i = 0; i < best_num_all; ++i)
    			{
    				fprintf(best, "K[%i] =\t%12g\n", i, fabs(par[i]));

    			}
    			fprintf(best,"\nEuclidean Norm:\t%12g\n",dat.norms[j][(best_num_all-1)] );
    			fprintf(best,"Objective Function:\t%12g\n", dat.Fobj[j][(best_num_all-1)]);
    		}
    		else
    		{
    			fprintf(best,"qmax =\t%12g\n", fabs(par[0]));
    			for (int i = 1; i < best_num_all; ++i)
    			{
    				fprintf(best, "K[%i] =\t%12g\n", i, fabs(par[i]));

    			}
    			fprintf(best,"\nEuclidean Norm:\t%12g\n",dat.norms[j][(best_num_all-2)] );
    			fprintf(best,"Objective Function:\t%12g\n", dat.Fobj[j][(best_num_all-2)]);
    		}
    		fprintf(best,"\nOptimized Results are as follows:\n");
    		for (int i = 0; i < m_dat[j]; ++i)
    		{
    			t[i] = dat.P[j].at(i);
    			y[i] = dat.q[j].at(i);
    			fprintf(best, "P[%i]=\t%12g\t q_obs=\t%12g\t q_eval=\t%12g\t Residual=\t%12g\n",
    					i, t[i], y[i], gstaFunc(t[i],par,dat.qmax,dat.n_par),
    						y[i] - gstaFunc(t[i],par,dat.qmax,dat.n_par) );
    		}
    		fprintf(best,"\n\n");

    		//Writes out the norms and Fobj to a file (to be removed later)
    		fprintf(norm,"\nIsotherm %i\n",(j+1));
    		fprintf(norm,"n_par \t\t norm \t\t Fobj\n");
    		if (dat.qmax!=0)
    		{
    			for (int i=0;i<best_num_all;i++)
				{
					fprintf(norm,"%i \t %12g \t %12g\n",(i+1),dat.norms[j].at(i),dat.Fobj[j].at(i));
				}
    		}
    		else
    		{
    			for (int i=1;i<best_num_all;i++)
				{
					fprintf(norm,"%i \t %12g \t %12g\n",(i+1),dat.norms[j].at((i-1)),dat.Fobj[j].at((i-1)));
				}
    		}
    		fprintf(params, "-----------------------------------------------"
    				"-------------------------------------------\n\n");

    		j++;
    		dGo.push_back(std::vector<double> ());
    	}while(j<n_iso);

    	//Estimate a single value for qmax if it is unknown during the first run through
    	if (dat.qmax == 0 && n_iso > 1)
    	{
    		//Declare some local parameters
    		double temp_norm[n_iso], temp_qmax[n_iso];
    		std::cout << "\nE.Norms and corresponding optimized qmax values\n";
    		std::cout << "--------------------------------------------------\n";
    		for (int j=0; j<n_iso; j++)
    		{
    			if(abs_max_q > dat.opt_qmax[j])
    			{
    				temp_norm[j] = fabs(abs_max_q - dat.opt_qmax[j]);
    				temp_qmax[j] = abs_max_q;
    			}
    			else
    			{
    				temp_norm[j] = fabs(avg_norm[j]);
    				temp_qmax[j] = fabs(dat.opt_qmax[j]);
    			}
    			std::cout << "Norm [" << (j+1) << "] = " << temp_norm[j] << "\t";
    			std::cout << "qmax = " << temp_qmax[j] << std::endl;
    		}
    		dat.qmax = weightedAvg(temp_norm, temp_qmax, n_iso);
    		std::cout << "\nRoutine to be re-run with the following estimate for qmax...";
    		std::cout << "\nWeighted Avg qmax = " << dat.qmax << std::endl << std::endl;
    	}
    	else
    	{
    		dat.qmax = fabs(dat.opt_qmax[0]);
    	}

    	//Close the files
    	fclose(output);
    	fclose(best);
    	fclose(norm);
    	fclose(params);
    	fclose(heats);

    	//Clearing specific memory
    	dat.Fobj.clear();
    	dat.best_par.clear();
    	dat.Kno.clear();
    	dat.all_pars.clear();
    	dat.norms.clear();
    	num_par.clear();
    	best_num_par.clear();

    }while(QMAX == 0);

    //Clear out memory all memory
    m_dat.clear();
    dat.Fobj.clear();
    dat.q.clear();
    dat.P.clear();
    dat.best_par.clear();
    dat.Kno.clear();
    dat.all_pars.clear();
    dat.norms.clear();
    iso_temp.clear();
    num_par.clear();
    best_num_par.clear();
    dat.opt_qmax.clear();
    dHo.clear();
    dSo.clear();
    dGo.clear();

    //Displays the total runtime of the program
    time = clock() - time;
    std::cout << "\nOptimization Runtime: " << (time / CLOCKS_PER_SEC) << " seconds" << std::endl;

    return 0;
}

