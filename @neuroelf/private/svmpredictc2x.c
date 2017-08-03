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

void read_sparse_instance(const mxArray *prhs, int index, struct svm_node *x)
{
	int i, j, low, high;
	mwIndex *ir, *jc;
	double *samples;

	ir = mxGetIr(prhs);
	jc = mxGetJc(prhs);
	samples = mxGetPr(prhs);

	/* each column is one instance */
	j = 0;
	low = (int)jc[index], high = (int)jc[index+1];
	for(i=low;i<high;i++)
	{
		x[j].index = (int)ir[i] + 1;
		x[j].value = samples[i];
		j++;
	}
	x[j].index = -1;
}

static void fake_answer(mxArray *plhs[])
{
	plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
}

void predict(mxArray *plhs[], const mxArray *prhs[], struct svm_model *model, const int predict_probability)
{
	int label_vector_row_num, label_vector_col_num;
	int feature_number, testing_instance_number;
	int instance_index;
	double *ptr_instance, *ptr_label, *ptr_predict_label;
	double *ptr_prob_estimates, *ptr_dec_values, *ptr;
	struct svm_node *x;
	mxArray *pplhs[1]; /* transposed instance sparse matrix */

	int correct = 0;
	int total = 0;
	double error = 0;
	double sump = 0, sumt = 0, sumpp = 0, sumtt = 0, sumpt = 0;

	int svm_type=svm_get_svm_type(model);
	int nr_class=svm_get_nr_class(model);
	double *prob_estimates=NULL;

	/* prhs[1] = testing instance matrix */
	feature_number = (int)mxGetN(prhs[1]);
	testing_instance_number = (int)mxGetM(prhs[1]);
	label_vector_row_num = (int)mxGetM(prhs[0]);
	label_vector_col_num = (int)mxGetN(prhs[0]);

	if(label_vector_row_num!=testing_instance_number)
	{
		mexPrintf("Length of label vector does not match # of instances.\n");
		fake_answer(plhs);
		return;
	}
	if(label_vector_col_num!=1)
	{
		mexPrintf("label (1st argument) should be a vector (# of column is 1).\n");
		fake_answer(plhs);
		return;
	}

	ptr_instance = mxGetPr(prhs[1]);
	ptr_label    = mxGetPr(prhs[0]);

	/* transpose instance matrix */
	if(mxIsSparse(prhs[1]))
	{
		if(model->param.kernel_type == PRECOMPUTED)
		{
			/* precomputed kernel requires dense matrix, so we make one */
			mxArray *rhs[1], *lhs[1];
			rhs[0] = mxDuplicateArray(prhs[1]);
			if(mexCallMATLAB(1, lhs, 1, rhs, "full"))
			{
				mexPrintf("Error: cannot full testing instance matrix\n");
				fake_answer(plhs);
				return;
			}
			ptr_instance = mxGetPr(lhs[0]);
			mxDestroyArray(rhs[0]);
		}
		else
		{
			mxArray *pprhs[1];
			pprhs[0] = mxDuplicateArray(prhs[1]);
			if(mexCallMATLAB(1, pplhs, 1, pprhs, "transpose"))
			{
				mexPrintf("Error: cannot transpose testing instance matrix\n");
				fake_answer(plhs);
				return;
			}
		}
	}

	if(predict_probability)
	{
		if(svm_type==NU_SVR || svm_type==EPSILON_SVR)
			mexPrintf("Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma=%g\n",svm_get_svr_probability(model));
		else
			prob_estimates = (double *) malloc(nr_class*sizeof(double));
	}

	plhs[0] = mxCreateDoubleMatrix(testing_instance_number, 1, mxREAL);
	if(predict_probability)
	{
		/* prob estimates are in plhs[2] */
		if(svm_type==C_SVC || svm_type==NU_SVC)
			plhs[2] = mxCreateDoubleMatrix(testing_instance_number, nr_class, mxREAL);
		else
			plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
	}
	else
	{
		/* decision values are in plhs[2] */
		if(svm_type == ONE_CLASS ||
		   svm_type == EPSILON_SVR ||
		   svm_type == NU_SVR)
			plhs[2] = mxCreateDoubleMatrix(testing_instance_number, 1, mxREAL);
		else
			plhs[2] = mxCreateDoubleMatrix(testing_instance_number, nr_class*(nr_class-1)/2, mxREAL);
	}

	ptr_predict_label = mxGetPr(plhs[0]);
	ptr_prob_estimates = mxGetPr(plhs[2]);
	ptr_dec_values = mxGetPr(plhs[2]);
	x = (struct svm_node*)malloc((feature_number+1)*sizeof(struct svm_node) );
	for(instance_index=0;instance_index<testing_instance_number;instance_index++)
	{
		int i;
		double target_label, predict_label;

		target_label = ptr_label[instance_index];

		if(mxIsSparse(prhs[1]) && model->param.kernel_type != PRECOMPUTED) /* prhs[1]^T is still sparse */
			read_sparse_instance(pplhs[0], instance_index, x);
		else
		{
			for(i=0;i<feature_number;i++)
			{
				x[i].index = i+1;
				x[i].value = ptr_instance[testing_instance_number*i+instance_index];
			}
			x[feature_number].index = -1;
		}

		if(predict_probability)
		{
			if(svm_type==C_SVC || svm_type==NU_SVC)
			{
				predict_label = svm_predict_probability(model, x, prob_estimates);
				ptr_predict_label[instance_index] = predict_label;
				for(i=0;i<nr_class;i++)
					ptr_prob_estimates[instance_index + i * testing_instance_number] = prob_estimates[i];
			} else {
				predict_label = svm_predict(model,x);
				ptr_predict_label[instance_index] = predict_label;
			}
		}
		else
		{
			predict_label = svm_predict(model,x);
			ptr_predict_label[instance_index] = predict_label;

			if(svm_type == ONE_CLASS ||
			   svm_type == EPSILON_SVR ||
			   svm_type == NU_SVR)
			{
				double res;
				svm_predict_values(model, x, &res);
				ptr_dec_values[instance_index] = res;
			}
			else
			{
				double *dec_values = (double *) malloc(sizeof(double) * nr_class*(nr_class-1)/2);
				svm_predict_values(model, x, dec_values);
				for(i=0;i<(nr_class*(nr_class-1))/2;i++)
					ptr_dec_values[instance_index + i * testing_instance_number] = dec_values[i];
				free(dec_values);
			}
		}

		if(predict_label == target_label)
			++correct;
		error += (predict_label-target_label)*(predict_label-target_label);
		sump += predict_label;
		sumt += target_label;
		sumpp += predict_label*predict_label;
		sumtt += target_label*target_label;
		sumpt += predict_label*target_label;
		++total;
	}
#if 0
    if(svm_type==NU_SVR || svm_type==EPSILON_SVR)
	{
		mexPrintf("Mean squared error = %g (regression)\n",error/total);
		mexPrintf("Squared correlation coefficient = %g (regression)\n",
			((total*sumpt-sump*sumt)*(total*sumpt-sump*sumt))/
			((total*sumpp-sump*sump)*(total*sumtt-sumt*sumt))
			);
	}
	else
		mexPrintf("Accuracy = %g%% (%d/%d) (classification)\n",
			(double)correct/total*100,correct,total);
#endif

	/* return accuracy, mean squared error, squared correlation coefficient */
	plhs[1] = mxCreateDoubleMatrix(3, 1, mxREAL);
	ptr = mxGetPr(plhs[1]);
	ptr[0] = (double)correct/total*100;
	ptr[1] = error/total;
	ptr[2] = ((total*sumpt-sump*sumt)*(total*sumpt-sump*sumt))/
				((total*sumpp-sump*sump)*(total*sumtt-sumt*sumt));

	free(x);
	if(prob_estimates != NULL)
		free(prob_estimates);
}

void exit_with_help()
{
	mexPrintf(
		"Usage: [predicted_label, accuracy, decision_values/prob_estimates] = svmpredict(testing_label_vector, testing_instance_matrix, model, 'libsvm_options')\n"
		"Parameters:\n"
		"  model: SVM model structure from svmtrain.\n"
		"  libsvm_options:\n"
		"    -b probability_estimates: whether to predict probability estimates, 0 or 1 (default 0); one-class SVM not supported yet\n"
		"Returns:\n"
		"  predicted_label: SVM prediction output vector.\n"
		"  accuracy: a vector with accuracy, mean squared error, squared correlation coefficient.\n"
		"  prob_estimates: If selected, probability estimate vector.\n"
	);
}

void mexFunction( int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
{
	int prob_estimate_flag = 0;
	struct svm_model *model;

	if(nrhs > 4 || nrhs < 3)
	{
		exit_with_help();
		fake_answer(plhs);
		return;
	}

	if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])) {
		mexPrintf("Error: label vector and instance matrix must be double\n");
		fake_answer(plhs);
		return;
	}

	if(mxIsStruct(prhs[2]))
	{
		const char *error_msg;

		/* parse options */
		if(nrhs==4)
		{
			int i, argc = 1;
			char cmd[CMD_LEN], *argv[CMD_LEN/2];

			/* put options in argv[] */
			mxGetString(prhs[3], cmd,  mxGetN(prhs[3]) + 1);
			if((argv[argc] = strtok(cmd, " ")) != NULL)
				while((argv[++argc] = strtok(NULL, " ")) != NULL)
					;

			for(i=1;i<argc;i++)
			{
				if(argv[i][0] != '-') break;
				if(++i>=argc)
				{
					exit_with_help();
					fake_answer(plhs);
					return;
				}
				switch(argv[i-1][1])
				{
					case 'b':
						prob_estimate_flag = atoi(argv[i]);
						break;
					default:
						mexPrintf("Unknown option: -%c\n", argv[i-1][1]);
						exit_with_help();
						fake_answer(plhs);
						return;
				}
			}
		}

		model = matlab_matrix_to_model(prhs[2], &error_msg);
		if (model == NULL)
		{
			mexPrintf("Error: can't read model: %s\n", error_msg);
			fake_answer(plhs);
			return;
		}

		if(prob_estimate_flag)
		{
			if(svm_check_probability_model(model)==0)
			{
				mexPrintf("Model does not support probabiliy estimates\n");
				fake_answer(plhs);
				svm_destroy_model(model);
				return;
			}
		}
		else
		{
			if(svm_check_probability_model(model)!=0)
				printf("Model supports probability estimates, but disabled in predicton.\n");
		}

		predict(plhs, prhs, model, prob_estimate_flag);
		/* destroy model */
		svm_destroy_model(model);
	}
	else
	{
		mexPrintf("model file should be a struct array\n");
		fake_answer(plhs);
	}

	return;
}
