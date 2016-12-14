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
inp = eql_ex0_load_records('gabor');

%% *MATERIAL DEFINITION*
mat = eql_ex1_mlp_vel;

%% *BATCH ANALYSIS*
eql = eql_ex0_run_vel(inp,mat);

%% *PLOT FIGURES*
eql_ex0_plot_res(inp,eql);