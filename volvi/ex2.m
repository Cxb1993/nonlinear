%% *1D SEISMIC WAVE PROPAGATION IN EQUIVALENT LINEAR SOIL MATERIAL*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _eql_exo0_: function to test EQL-ECP code
%% *N.B.*
% Need for:

%% *SET UP*
ccc;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
%% *PARSE INPUT MOTIONS*
inp = eql_load_records('argost_E83.cyt');

%% *MATERIAL DEFINITION*
mat = eql_mat_ex2;

%% *BATCH ANALYSIS*
eql = eql_run(inp,mat);

%% *PLOT FIGURES*
eql_plot_res_ex2(inp,eql);