/*

% Version:  v0.9d
% Build:    14061709
% Date:     Jun-17 2014, 9:30 AM EST
% Author:   Chih-Chung Chang and Chih-Jen Lin, NTU CSIE, Taiwan
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://www.csie.ntu.edu.tw/~cjlin/libsvm/

*/

const char *model_to_matlab_structure(mxArray *plhs[], int num_of_feature, struct svm_model *model);
struct svm_model *matlab_matrix_to_model(const mxArray *matlab_struct, const char **error_message);
