# sparse_tv_MVAR
The file sparse_time_varying_MVAR_EMpara.m contains implementation of the algorithm usd in WCH's thesis.
Other files contain other algorithms to estimate sparse time varying MVAR (EM & SAEM)  

The method in sparse_time_varying_MVAR_EMpara "estimate_model_VB" is to estimate the model using VB. Here is a demo.

     stv_mvar = sparse_time_varying_MVAR_EMpara(1,[],[],[], data,order); % data is the data to analyzel; order is the model order
     [A, lb]= stv_mvar.estimate_model_VB(); % A is the time varying MVAR model; lb is the variational lower bound

If you have any problem please contact me with xs.wuchaohua@gmail.com 
