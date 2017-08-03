/*

% Version:  v0.9b
% Build:    11050223
% Date:     Mar-02 2011, 6:05 PM EST
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "svm2x.h"

#include "mex.h"

#ifdef MX_API_VER
#if MX_API_VER < 0x07030000
typedef int mwIndex;
#endif
#endif

#define NUM_OF_RETURN_FIELD 10

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

static const char *field_names[] = {
	"Parameters",
	"nr_class",
	"totalSV",
	"rho",
	"Label",
	"ProbA",
	"ProbB",
	"nSV",
	"sv_coef",
	"SVs"
};

const char *model_to_matlab_structure(mxArray *plhs[], int num_of_feature, struct svm_model *model)
{
	int i, j, n;
	double *ptr;
	mxArray *return_model, **rhs;
	int out_id = 0;

	rhs = (mxArray **)mxMalloc(sizeof(mxArray *)*NUM_OF_RETURN_FIELD);

	/* Parameters */
	rhs[out_id] = mxCreateDoubleMatrix(5, 1, mxREAL);
	ptr = mxGetPr(rhs[out_id]);
	ptr[0] = model->param.svm_type;
	ptr[1] = model->param.kernel_type;
	ptr[2] = model->param.degree;
	ptr[3] = model->param.gamma;
	ptr[4] = model->param.coef0;
	out_id++;

	/* nr_class */
	rhs[out_id] = mxCreateDoubleMatrix(1, 1, mxREAL);
	ptr = mxGetPr(rhs[out_id]);
	ptr[0] = model->nr_class;
	out_id++;

	/* total SV */
	rhs[out_id] = mxCreateDoubleMatrix(1, 1, mxREAL);
	ptr = mxGetPr(rhs[out_id]);
	ptr[0] = model->l;
	out_id++;

	/* rho */
	n = model->nr_class*(model->nr_class-1)/2;
	rhs[out_id] = mxCreateDoubleMatrix(n, 1, mxREAL);
	ptr = mxGetPr(rhs[out_id]);
	for(i = 0; i < n; i++)
		ptr[i] = model->rho[i];
	out_id++;

	/* Label */
	if(model->label)
	{
		rhs[out_id] = mxCreateDoubleMatrix(model->nr_class, 1, mxREAL);
		ptr = mxGetPr(rhs[out_id]);
		for(i = 0; i < model->nr_class; i++)
			ptr[i] = model->label[i];
	}
	else
		rhs[out_id] = mxCreateDoubleMatrix(0, 0, mxREAL);
	out_id++;

	/* probA */
	if(model->probA != NULL)
	{
		rhs[out_id] = mxCreateDoubleMatrix(n, 1, mxREAL);
		ptr = mxGetPr(rhs[out_id]);
		for(i = 0; i < n; i++)
			ptr[i] = model->probA[i];
	}
	else
		rhs[out_id] = mxCreateDoubleMatrix(0, 0, mxREAL);
	out_id ++;

	/* probB */
	if(model->probB != NULL)
	{
		rhs[out_id] = mxCreateDoubleMatrix(n, 1, mxREAL);
		ptr = mxGetPr(rhs[out_id]);
		for(i = 0; i < n; i++)
			ptr[i] = model->probB[i];
	}
	else
		rhs[out_id] = mxCreateDoubleMatrix(0, 0, mxREAL);
	out_id++;

	/* nSV */
	if(model->nSV)
	{
		rhs[out_id] = mxCreateDoubleMatrix(model->nr_class, 1, mxREAL);
		ptr = mxGetPr(rhs[out_id]);
		for(i = 0; i < model->nr_class; i++)
			ptr[i] = model->nSV[i];
	}
	else
		rhs[out_id] = mxCreateDoubleMatrix(0, 0, mxREAL);
	out_id++;

	/* sv_coef */
	rhs[out_id] = mxCreateDoubleMatrix(model->l, model->nr_class-1, mxREAL);
	ptr = mxGetPr(rhs[out_id]);
	for(i = 0; i < model->nr_class-1; i++)
		for(j = 0; j < model->l; j++)
			ptr[(i*(model->l))+j] = model->sv_coef[i][j];
	out_id++;

	/* SVs */
	{
		int ir_index, nonzero_element;
		mwIndex *ir, *jc;
		mxArray *pprhs[1], *pplhs[1];

		if(model->param.kernel_type == PRECOMPUTED)
		{
			nonzero_element = model->l;
			num_of_feature = 1;
		}
		else
		{
			nonzero_element = 0;
			for(i = 0; i < model->l; i++) {
				j = 0;
				while(model->SV[i][j].index != -1)
				{
					nonzero_element++;
					j++;
				}
			}
		}

		/* SV in column, easier accessing */
		rhs[out_id] = mxCreateSparse(num_of_feature, model->l, nonzero_element, mxREAL);
		ir = mxGetIr(rhs[out_id]);
		jc = mxGetJc(rhs[out_id]);
		ptr = mxGetPr(rhs[out_id]);
		jc[0] = ir_index = 0;
		for(i = 0;i < model->l; i++)
		{
			if(model->param.kernel_type == PRECOMPUTED)
			{
				/* make a (1 x model->l) matrix */
				ir[ir_index] = 0;
				ptr[ir_index] = model->SV[i][0].value;
				ir_index++;
				jc[i+1] = jc[i] + 1;
			}
			else
			{
				int x_index = 0;
				while (model->SV[i][x_index].index != -1)
				{
					ir[ir_index] = model->SV[i][x_index].index - 1;
					ptr[ir_index] = model->SV[i][x_index].value;
					ir_index++, x_index++;
				}
				jc[i+1] = jc[i] + x_index;
			}
		}
		/* transpose back to SV in row */
		pprhs[0] = rhs[out_id];
		if(mexCallMATLAB(1, pplhs, 1, pprhs, "transpose"))
			return "cannot transpose SV matrix";
		rhs[out_id] = pplhs[0];
		out_id++;
	}

	/* Create a struct matrix contains NUM_OF_RETURN_FIELD fields */
	return_model = mxCreateStructMatrix(1, 1, NUM_OF_RETURN_FIELD, field_names);

	/* Fill struct matrix with input arguments */
	for(i = 0; i < NUM_OF_RETURN_FIELD; i++)
		mxSetField(return_model,0,field_names[i],mxDuplicateArray(rhs[i]));
	/* return */
	plhs[0] = return_model;
	mxFree(rhs);

	return NULL;
}

struct svm_model *matlab_matrix_to_model(const mxArray *matlab_struct, const char **msg)
{
	int i, j, n, num_of_fields;
	double *ptr;
	int id = 0;
	struct svm_node *x_space;
	struct svm_model *model;
	mxArray **rhs;

	num_of_fields = mxGetNumberOfFields(matlab_struct);
	if(num_of_fields != NUM_OF_RETURN_FIELD)
	{
		*msg = "number of return field is not correct";
		return NULL;
	}
	rhs = (mxArray **) mxMalloc(sizeof(mxArray *)*num_of_fields);

	for(i=0;i<num_of_fields;i++)
		rhs[i] = mxGetFieldByNumber(matlab_struct, 0, i);

	model = Malloc(struct svm_model, 1);
	model->rho = NULL;
	model->probA = NULL;
	model->probB = NULL;
	model->label = NULL;
	model->nSV = NULL;
	model->free_sv = 1; /* XXX */

	ptr = mxGetPr(rhs[id]);
	model->param.svm_type = (int)ptr[0];
	model->param.kernel_type  = (int)ptr[1];
	model->param.degree	  = (int)ptr[2];
	model->param.gamma	  = ptr[3];
	model->param.coef0	  = ptr[4];
	id++;

	ptr = mxGetPr(rhs[id]);
	model->nr_class = (int)ptr[0];
	id++;

	ptr = mxGetPr(rhs[id]);
	model->l = (int)ptr[0];
	id++;

	/* rho */
	n = model->nr_class * (model->nr_class-1)/2;
	model->rho = (double*) malloc(n*sizeof(double));
	ptr = mxGetPr(rhs[id]);
	for(i=0;i<n;i++)
		model->rho[i] = ptr[i];
	id++;

	/* label */
	if(mxIsEmpty(rhs[id]) == 0)
	{
		model->label = (int*) malloc(model->nr_class*sizeof(int));
		ptr = mxGetPr(rhs[id]);
		for(i=0;i<model->nr_class;i++)
			model->label[i] = (int)ptr[i];
	}
	id++;

	/* probA */
	if(mxIsEmpty(rhs[id]) == 0)
	{
		model->probA = (double*) malloc(n*sizeof(double));
		ptr = mxGetPr(rhs[id]);
		for(i=0;i<n;i++)
			model->probA[i] = ptr[i];
	}
	id++;

	/* probB */
	if(mxIsEmpty(rhs[id]) == 0)
	{
		model->probB = (double*) malloc(n*sizeof(double));
		ptr = mxGetPr(rhs[id]);
		for(i=0;i<n;i++)
			model->probB[i] = ptr[i];
	}
	id++;

	/* nSV */
	if(mxIsEmpty(rhs[id]) == 0)
	{
		model->nSV = (int*) malloc(model->nr_class*sizeof(int));
		ptr = mxGetPr(rhs[id]);
		for(i=0;i<model->nr_class;i++)
			model->nSV[i] = (int)ptr[i];
	}
	id++;

	/* sv_coef */
	ptr = mxGetPr(rhs[id]);
	model->sv_coef = (double**) malloc((model->nr_class-1)*sizeof(double));
	for( i=0 ; i< model->nr_class -1 ; i++ )
		model->sv_coef[i] = (double*) malloc((model->l)*sizeof(double));
	for(i = 0; i < model->nr_class - 1; i++)
		for(j = 0; j < model->l; j++)
			model->sv_coef[i][j] = ptr[i*(model->l)+j];
	id++;

	/* SV */
	{
		int sr, sc, elements;
		int num_samples;
		mwIndex *ir, *jc;
		mxArray *pprhs[1], *pplhs[1];

		/* transpose SV */
		pprhs[0] = rhs[id];
		if(mexCallMATLAB(1, pplhs, 1, pprhs, "transpose"))
		{
			svm_destroy_model(model);
			*msg = "cannot transpose SV matrix";
			return NULL;
		}
		rhs[id] = pplhs[0];

		sr = (int)mxGetN(rhs[id]);
		sc = (int)mxGetM(rhs[id]);

		ptr = mxGetPr(rhs[id]);
		ir = mxGetIr(rhs[id]);
		jc = mxGetJc(rhs[id]);

		num_samples = (int)mxGetNzmax(rhs[id]);

		elements = num_samples + sr;

		model->SV = (struct svm_node **) malloc(sr * sizeof(struct svm_node *));
		x_space = (struct svm_node *)malloc(elements * sizeof(struct svm_node));

		/* SV is in column */
		for(i=0;i<sr;i++)
		{
			int low = (int)jc[i], high = (int)jc[i+1];
			int x_index = 0;
			model->SV[i] = &x_space[low+i];
			for(j=low;j<high;j++)
			{
				model->SV[i][x_index].index = (int)ir[j] + 1;
				model->SV[i][x_index].value = ptr[j];
				x_index++;
			}
			model->SV[i][x_index].index = -1;
		}

		id++;
	}
	mxFree(rhs);

	return model;
}


#define CMD_LEN 2048
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

void print_null(const char *s) {}

void exit_with_help()
{
	mexPrintf(
	"Usage: model = svmtrain(training_label_vector, training_instance_matrix, 'libsvm_options');\n"
	"libsvm_options:\n"
	"-s svm_type : set type of SVM (default 0)\n"
	"	0 -- C-SVC\n"
	"	1 -- nu-SVC\n"
	"	2 -- one-class SVM\n"
	"	3 -- epsilon-SVR\n"
	"	4 -- nu-SVR\n"
	"-t kernel_type : set type of kernel function (default 2)\n"
	"	0 -- linear: u'*v\n"
	"	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
	"	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
	"	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
	"	4 -- precomputed kernel (kernel values in training_instance_matrix)\n"
	"-d degree : set degree in kernel function (default 3)\n"
	"-g gamma : set gamma in kernel function (default 1/k)\n"
	"-r coef0 : set coef0 in kernel function (default 0)\n"
	"-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"
	"-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
	"-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
	"-m cachesize : set cache memory size in MB (default 100)\n"
	"-e epsilon : set tolerance of termination criterion (default 0.001)\n"
	"-h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
	"-b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
	"-wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
	"-v n : n-fold cross validation mode\n"
	"-q : quiet mode (no outputs)\n"
	);
}

/* svm arguments */
struct svm_parameter param;		/* set by parse_command_line */
struct svm_problem prob;		/* set by read_problem */
struct svm_model *model;
struct svm_node *x_space;
int cross_validation;
int nr_fold;

void (*svm_default_print_string) (const char *) = NULL;

double do_cross_validation()
{
	int i;
	int total_correct = 0;
	double total_error = 0;
	double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
	double *target = Malloc(double,prob.l);
	double retval = 0.0;

	svm_cross_validation(&prob,&param,nr_fold,target);
	if(param.svm_type == EPSILON_SVR ||
	   param.svm_type == NU_SVR)
	{
		for(i=0;i<prob.l;i++)
		{
			double y = prob.y[i];
			double v = target[i];
			total_error += (v-y)*(v-y);
			sumv += v;
			sumy += y;
			sumvv += v*v;
			sumyy += y*y;
			sumvy += v*y;
		}
		mexPrintf("Cross Validation Mean squared error = %g\n",total_error/prob.l);
		mexPrintf("Cross Validation Squared correlation coefficient = %g\n",
			((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/
			((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy))
			);
		retval = total_error/prob.l;
	}
	else
	{
		for(i=0;i<prob.l;i++)
			if(target[i] == prob.y[i])
				++total_correct;
		mexPrintf("Cross Validation Accuracy = %g%%\n",100.0*total_correct/prob.l);
		retval = 100.0*total_correct/prob.l;
	}
	free(target);
	return retval;
}

/* nrhs should be 3 */
int parse_command_line(int nrhs, const mxArray *prhs[], char *model_file_name)
{
	int i, argc = 1;
	char cmd[CMD_LEN];
	char *argv[CMD_LEN/2];

	/* default values */
	param.svm_type = C_SVC;
	param.kernel_type = RBF;
	param.degree = 3;
	param.gamma = 0;	/* 1/k */
	param.coef0 = 0;
	param.nu = 0.5;
	param.cache_size = 100;
	param.C = 1;
	param.eps = 1e-3;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
	cross_validation = 0;
	/* svmtrain loaded only once under matlab */
	if (svm_default_print_string == NULL)
		svm_default_print_string = svm_print_string;
	else
		svm_print_string = svm_default_print_string;


	if(nrhs <= 1)
		return 1;

	if(nrhs > 2)
	{
		/* put options in argv[] */
		mxGetString(prhs[2], cmd, mxGetN(prhs[2]) + 1);
		if((argv[argc] = strtok(cmd, " ")) != NULL)
			while((argv[++argc] = strtok(NULL, " ")) != NULL)
				;
	}

	/* parse options */
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break;
		++i;
		if(i>=argc && argv[i-1][1] != 'q')	/* since option -q has no parameter */
			return 1;
		switch(argv[i-1][1])
		{
			case 's':
				param.svm_type = atoi(argv[i]);
				break;
			case 't':
				param.kernel_type = atoi(argv[i]);
				break;
			case 'd':
				param.degree = atoi(argv[i]);
				break;
			case 'g':
				param.gamma = atof(argv[i]);
				break;
			case 'r':
				param.coef0 = atof(argv[i]);
				break;
			case 'n':
				param.nu = atof(argv[i]);
				break;
			case 'm':
				param.cache_size = atof(argv[i]);
				break;
			case 'c':
				param.C = atof(argv[i]);
				break;
			case 'e':
				param.eps = atof(argv[i]);
				break;
			case 'p':
				param.p = atof(argv[i]);
				break;
			case 'h':
				param.shrinking = atoi(argv[i]);
				break;
			case 'b':
				param.probability = atoi(argv[i]);
				break;
			case 'q':
				svm_print_string = &print_null;
				i--;
				break;
			case 'v':
				cross_validation = 1;
				nr_fold = atoi(argv[i]);
				if(nr_fold < 2)
				{
					mexPrintf("n-fold cross validation: n must >= 2\n");
					return 1;
				}
				break;
			case 'w':
				++param.nr_weight;
				param.weight_label = (int *)realloc(param.weight_label,sizeof(int)*param.nr_weight);
				param.weight = (double *)realloc(param.weight,sizeof(double)*param.nr_weight);
				param.weight_label[param.nr_weight-1] = atoi(&argv[i-1][2]);
				param.weight[param.nr_weight-1] = atof(argv[i]);
				break;
			default:
				mexPrintf("Unknown option -%c\n", argv[i-1][1]);
				return 1;
		}
	}
	return 0;
}

/* read in a problem (in svmlight format) */
int read_problem_dense(const mxArray *label_vec, const mxArray *instance_mat)
{
	int i, j, k;
	int elements, max_index, sc, label_vector_row_num;
	double *samples, *labels;

	prob.x = NULL;
	prob.y = NULL;
	x_space = NULL;

	labels = mxGetPr(label_vec);
	samples = mxGetPr(instance_mat);
	sc = (int)mxGetN(instance_mat);

	elements = 0;
	/* the number of instance */
	prob.l = (int)mxGetM(instance_mat);
	label_vector_row_num = (int)mxGetM(label_vec);

	if(label_vector_row_num!=prob.l)
	{
		mexPrintf("Length of label vector does not match # of instances.\n");
		return -1;
	}

	if(param.kernel_type == PRECOMPUTED)
		elements = prob.l * (sc + 1);
	else
	{
		for(i = 0; i < prob.l; i++)
		{
			for(k = 0; k < sc; k++)
				if(samples[k * prob.l + i] != 0)
					elements++;
			/* count the '-1' element */
			elements++;
		}
	}

	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct svm_node *,prob.l);
	x_space = Malloc(struct svm_node, elements);

	max_index = sc;
	j = 0;
	for(i = 0; i < prob.l; i++)
	{
		prob.x[i] = &x_space[j];
		prob.y[i] = labels[i];

		for(k = 0; k < sc; k++)
		{
			if(param.kernel_type == PRECOMPUTED || samples[k * prob.l + i] != 0)
			{
				x_space[j].index = k + 1;
				x_space[j].value = samples[k * prob.l + i];
				j++;
			}
		}
		x_space[j++].index = -1;
	}

	if(param.gamma == 0 && max_index > 0)
		param.gamma = 1.0/max_index;

	if(param.kernel_type == PRECOMPUTED)
		for(i=0;i<prob.l;i++)
		{
			if((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index)
			{
				mexPrintf("Wrong input format: sample_serial_number out of range\n");
				return -1;
			}
		}

	return 0;
}

int read_problem_sparse(const mxArray *label_vec, const mxArray *instance_mat)
{
	int i, j, k, low, high;
	mwIndex *ir, *jc;
	int elements, max_index, num_samples, label_vector_row_num;
	double *samples, *labels;
	mxArray *instance_mat_col; /* transposed instance sparse matrix */

	prob.x = NULL;
	prob.y = NULL;
	x_space = NULL;

	/* transpose instance matrix */
	{
		mxArray *prhs[1], *plhs[1];
		prhs[0] = mxDuplicateArray(instance_mat);
		if(mexCallMATLAB(1, plhs, 1, prhs, "transpose"))
		{
			mexPrintf("Error: cannot transpose training instance matrix\n");
			return -1;
		}
		instance_mat_col = plhs[0];
		mxDestroyArray(prhs[0]);
	}

	/* each column is one instance */
	labels = mxGetPr(label_vec);
	samples = mxGetPr(instance_mat_col);
	ir = mxGetIr(instance_mat_col);
	jc = mxGetJc(instance_mat_col);

	num_samples = (int)mxGetNzmax(instance_mat_col);

	/* the number of instance */
	prob.l = (int)mxGetN(instance_mat_col);
	label_vector_row_num = (int)mxGetM(label_vec);

	if(label_vector_row_num!=prob.l)
	{
		mexPrintf("Length of label vector does not match # of instances.\n");
		return -1;
	}

	elements = num_samples + prob.l;
	max_index = (int)mxGetM(instance_mat_col);

	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct svm_node *,prob.l);
	x_space = Malloc(struct svm_node, elements);

	j = 0;
	for(i=0;i<prob.l;i++)
	{
		prob.x[i] = &x_space[j];
		prob.y[i] = labels[i];
		low = (int)jc[i], high = (int)jc[i+1];
		for(k=low;k<high;k++)
		{
			x_space[j].index = (int)ir[k] + 1;
			x_space[j].value = samples[k];
			j++;
	 	}
		x_space[j++].index = -1;
	}

	if(param.gamma == 0 && max_index > 0)
		param.gamma = 1.0/max_index;

	return 0;
}

static void fake_answer(mxArray *plhs[])
{
	plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
}

/* Interface function of matlab */
/* now assume prhs[0]: label prhs[1]: features */
void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] )
{
	const char *error_msg;

	/* fix random seed to have same results for each run */
	/* (for cross validation and probability estimation) */
	srand(1);

	/* Transform the input Matrix to libsvm format */
	if(nrhs > 0 && nrhs < 4)
	{
		int err;

		if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])) {
			mexPrintf("Error: label vector and instance matrix must be double\n");
			fake_answer(plhs);
			return;
		}

		if(parse_command_line(nrhs, prhs, NULL))
		{
			exit_with_help();
			svm_destroy_param(&param);
			fake_answer(plhs);
			return;
		}

		if(mxIsSparse(prhs[1]))
		{
			if(param.kernel_type == PRECOMPUTED)
			{
				/* precomputed kernel requires dense matrix, so we make one */
				mxArray *rhs[1], *lhs[1];

				rhs[0] = mxDuplicateArray(prhs[1]);
				if(mexCallMATLAB(1, lhs, 1, rhs, "full"))
				{
					mexPrintf("Error: cannot generate a full training instance matrix\n");
					svm_destroy_param(&param);
					fake_answer(plhs);
					return;
				}
				err = read_problem_dense(prhs[0], lhs[0]);
				mxDestroyArray(lhs[0]);
				mxDestroyArray(rhs[0]);
			}
			else
				err = read_problem_sparse(prhs[0], prhs[1]);
		}
		else
			err = read_problem_dense(prhs[0], prhs[1]);

		/* svmtrain's original code */
		error_msg = svm_check_parameter(&prob, &param);

		if(err || error_msg)
		{
			if (error_msg != NULL)
				mexPrintf("Error: %s\n", error_msg);
			svm_destroy_param(&param);
			free(prob.y);
			free(prob.x);
			free(x_space);
			fake_answer(plhs);
			return;
		}

		if(cross_validation)
		{
			double *ptr;
			plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
			ptr = mxGetPr(plhs[0]);
			ptr[0] = do_cross_validation();
		}
		else
		{
			int nr_feat = (int)mxGetN(prhs[1]);
			const char *error_msg;
			model = svm_train(&prob, &param);
			error_msg = model_to_matlab_structure(plhs, nr_feat, model);
			if(error_msg)
				mexPrintf("Error: can't convert libsvm model to matrix structure: %s\n", error_msg);
			svm_destroy_model(model);
		}
		svm_destroy_param(&param);
		free(prob.y);
		free(prob.x);
		free(x_space);
	}
	else
	{
		exit_with_help();
		fake_answer(plhs);
		return;
	}
}
